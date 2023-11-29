# Copyright (C) 2023 Luigi Pertoldi <gipert@pm.me>
#
# This program is free software: you can redistribute it and/or modify it under
# the terms of the GNU Lesser General Public License as published by the Free
# Software Foundation, either version 3 of the License, or (at your option) any
# later version.
#
# This program is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more
# details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

# ruff: noqa: F821, T201

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
import ROOT
import uproot
from legendmeta import LegendMetadata


def process_mage_id(mage_ids):
    mage_names = {}
    for _mage_id in mage_ids:
        m_id = str(_mage_id)
        is_ged = bool(int(m_id[0]))
        if not is_ged:
            # This should never be the case
            continue

        string = int(m_id[3:5])
        pos = int(m_id[5:7])

        for _name, _meta_dict in chmap.items():
            if _meta_dict["system"] == "geds":
                location = _meta_dict["location"]
                if location["string"] == string and location["position"] == pos:
                    mage_names[_mage_id] = f"ch{_meta_dict['daq']['rawid']}"

    return mage_names


parser = argparse.ArgumentParser(
    prog="build_pdf", description="build LEGEND pdf files from evt tier files"
)
parser.add_argument(
    "--raw-files",
    "-r",
    default=None,
    nargs="+",
    help="path to raw simulation files for number-of-primaries determination",
)
parser.add_argument("--config", "-c", required=True, help="configuration file")
parser.add_argument("--output", "-o", required=True, help="output file name")
parser.add_argument("--metadata", "-m", required=True, help="path to legend-metadata")
parser.add_argument("input_files", nargs="+", help="evt tier files")

args = parser.parse_args()

if not isinstance(args.input_files, list):
    args.input_files = [args.input_files]

# with Path(args.config).open() as f:
with Path(args.config).open() as f:
    rconfig = json.load(f)

meta = LegendMetadata(args.metadata)
chmap = meta.channelmap(rconfig["timestamp"])
geds_mapping = {
    f"ch{_dict['daq']['rawid']}": _name
    for _name, _dict in chmap.items()
    if chmap[_name]["system"] == "geds"
}
n_primaries_total = 0

if args.raw_files:
    for file in args.raw_files:
        with uproot.open(f"{file}:fTree") as fTree:
            n_primaries_total += fTree["fNEvents"].array(entry_stop=1)[0]

# So there are many input files fed into one pdf file
# set up the hists to fill as we go along
# Creat a hist for all dets (even AC ones)
hists = {
    _cut_name: {
        _rawid: ROOT.TH1F(
            f"{_cut_name}_{_rawid}",
            f"{_name} energy deposits",
            rconfig["hist"]["nbins"],
            rconfig["hist"]["emin"],
            rconfig["hist"]["emax"],
        )
        for _rawid, _name in sorted(geds_mapping.items())
    }
    for _cut_name in rconfig["cuts"]
    if rconfig["cuts"][_cut_name]["is_sum"] is False
}

# When we want to start summing the energy of events we have to treat them differently
sum_hists = {
    _cut_name: ROOT.TH1F(
        f"{_cut_name}_all_summed",
        "summed energy deposits",
        rconfig["hist"]["nbins"],
        rconfig["hist"]["emin"],
        rconfig["hist"]["emax"],
    )
    for _cut_name in rconfig["cuts"]
    if rconfig["cuts"][_cut_name]["is_sum"] is True
}

for file_name in args.input_files:
    print(file_name)
    with uproot.open(f"{file_name}:simTree") as pytree:
        if pytree.num_entries == 0:
            msg = f">> Error: MPP evt file {file_name} has 0 events in simTree"
            raise Exception(msg)

        n_primaries = pytree["mage_n_events"].array()[0]
        df_data = pd.DataFrame(
            pytree.arrays(["energy", "npe_tot", "mage_id", "is_good"], library="np")
        )

    # add a column with Poisson(mu=npe_tot) to represent the actual random
    # number of detected photons. This column should be used to determine
    # the LAr classifier
    rng = np.random.default_rng()
    df_data["npe_tot_poisson"] = rng.poisson(df_data.npe_tot)

    df_exploded = df_data.explode(["energy", "mage_id", "is_good"])
    df_ecut = df_exploded.query(f"energy > {rconfig['energy_threshold']}")

    # These give you the multiplicity of events (and events not including AC detectors)
    index_counts = df_ecut.index.value_counts()
    index_counts_is_good = df_ecut.query("is_good == True").index.value_counts()

    df_ecut = df_ecut.copy()
    df_ecut["mul"] = df_ecut.index.map(index_counts)
    df_ecut["mul_is_good"] = df_ecut.index.map(index_counts_is_good)

    if not args.raw_files:
        n_primaries_total += n_primaries

    uniq_mage_ids = df_exploded.dropna(subset=["mage_id"])["mage_id"].unique()
    mage_names = process_mage_id(uniq_mage_ids)
    # out_file["number_of_primaries"] = str(n_primaries)

    for _cut_name, _cut_dict in rconfig["cuts"].items():
        # We want to cut on multiplicity for all detectors >25keV
        # Include them in the dataset then apply cuts - then filter them out
        # Don't store AC detectors
        _cut_string = _cut_dict["cut_string"]
        df_cut = df_ecut.copy() if _cut_string == "" else df_ecut.query(_cut_string)

        df_good = df_cut[df_cut.is_good == True]  # noqa: E712

        if _cut_dict["is_sum"] is False:
            for __mage_id in df_good.mage_id.unique():
                _rawid = mage_names[__mage_id]
                _energy_array = (
                    df_good.query(f"mage_id == {__mage_id}").energy.to_numpy(
                        dtype=float
                    )
                    * 1000
                )  # keV

                hists[_cut_name][_rawid].FillN(
                    len(_energy_array), _energy_array, np.ones(len(_energy_array))
                )
        else:
            _summed_energy_array = (
                df_good.groupby(df_good.index).energy.sum().to_numpy(dtype=float) * 1000
            )  # keV
            sum_hists[_cut_name].FillN(
                len(_summed_energy_array),
                _summed_energy_array,
                np.ones(len(_summed_energy_array)),
            )

# The individual channels have been filled
# now add them together to make the grouped hists
# We don't need to worry about the ac dets as they will have zero entries
for _cut_name in hists:
    hists[_cut_name]["all"] = ROOT.TH1F(
        f"{_cut_name}_all",
        "All energy deposits",
        rconfig["hist"]["nbins"],
        rconfig["hist"]["emin"],
        rconfig["hist"]["emax"],
    )
    for _type in ["bege", "coax", "icpc", "ppc"]:
        hists[_cut_name][_type] = ROOT.TH1F(
            f"{_cut_name}_{_type}",
            f"All {_type} energy deposits",
            rconfig["hist"]["nbins"],
            rconfig["hist"]["emin"],
            rconfig["hist"]["emax"],
        )
    for _rawid, _name in geds_mapping.items():
        hists[_cut_name][chmap[geds_mapping[_rawid]]["type"]].Add(
            hists[_cut_name][_rawid]
        )
        hists[_cut_name]["all"].Add(hists[_cut_name][_rawid])

# write the hists to file (but only if they have none zero entries)
# Changes the names to drop type_ etc
out_file = uproot.recreate(args.output)
for _cut_name, _hist_dict in hists.items():
    dir = out_file.mkdir(_cut_name)
    for key, item in _hist_dict.items():
        if item.GetEntries() > 0:
            dir[key] = item

for _cut_name, _hist in sum_hists.items():
    dir = out_file.mkdir(_cut_name)
    if _hist.GetEntries() > 0:
        dir["all"] = _hist

out_file["number_of_primaries"] = str(int(n_primaries_total))
out_file.close()
