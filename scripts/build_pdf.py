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
import time
from pathlib import Path

import pandas as pd
import ROOT
import uproot
from legendmeta import LegendMetadata

start = time.time()


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
                    mage_names[f"ch{_meta_dict['daq']['rawid']}"] = _mage_id

    return mage_names


parser = argparse.ArgumentParser(
    prog="build_pdf", description="build LEGEND pdf files from evt tier files"
)
parser.add_argument("--config", "-c", help="configuration file")
parser.add_argument("--output", "-o", help="output file name")
parser.add_argument("--metadata", "-m", help="path to legend-metadata")
parser.add_argument("input_files", nargs="+", help="evt tier files")

args = parser.parse_args()

if not isinstance(args.input_files, list):
    args.input_files = [args.input_files]

# with Path(args.config).open() as f:
with Path(args.config).open() as f:
    rconfig = json.load(f)

meta = LegendMetadata(args.metadata)
chmap = meta.channelmap(rconfig["timestamp"])
usable_geds = {
    f"ch{_dict['daq']['rawid']}": _name
    for _name, _dict in chmap.items()
    if chmap[_name]["system"] == "geds"
    and (
        chmap[_name]["analysis"]["usability"] == "on"
        or chmap[_name]["analysis"]["usability"] == "no_psd"
    )
}
n_primaries_total = 0

# So there are many input files fed into one pdf file
# set up the hists to fill as we go along
# Don't store the AC dets
hists = {
    _cut_name: {
        _rawid: ROOT.TH1F(
            f"{_cut_name}_{_rawid}",
            f"{_name} energy deposits",
            rconfig["hist"]["nbins"],
            rconfig["hist"]["emin"],
            rconfig["hist"]["emax"],
        )
        for _rawid, _name in sorted(usable_geds.items())
    }
    for _cut_name in rconfig["cuts"]
}

for file_name in args.input_files:
    with uproot.open(f"{file_name}:simTree") as pytree:
        n_primaries = pytree["mage_n_events"].array()[0]
        df_data = pd.DataFrame(
            pytree.arrays(["energy", "npe_tot", "mage_id"], library="np")
        )
    df_exploded = df_data.explode("energy").explode("mage_id")

    df_ecut = df_exploded[df_exploded["energy"] > rconfig["energy_threshold"]]
    index_counts = df_ecut.index.value_counts()

    n_primaries_total += n_primaries

    uniq_mage_ids = df_exploded.dropna(subset=["mage_id"])["mage_id"].unique()
    mage_names = process_mage_id(uniq_mage_ids)
    # out_file["number_of_primaries"] = str(n_primaries)

    for _cut_name, _cut_string in rconfig["cuts"].items():
        # We want to cut on multiplicity for all detectors >25keV
        # Include them in the dataset then apply cuts - then filter them out
        # Don't store AC detectors
        exec(_cut_string)

        # loop over the geds in the file
        for _rawid, _mage_id in mage_names.items():
            if _rawid not in usable_geds.keys():
                continue
            df_channel = df_cut[df_cut.mage_id == _mage_id]

            for energy in df_channel.energy.to_numpy():
                hists[_cut_name][_rawid].Fill(energy * 1000)  # energy in keV

# The individual channels have been filled
# now add them together to make the grouped hists
# We don't need to worry about the ac dets
for _cut_name in rconfig["cuts"]:
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
    for _rawid, _name in usable_geds.items():
        hists[_cut_name][chmap[usable_geds[_rawid]]["type"]].Add(
            hists[_cut_name][_rawid]
        )
        hists[_cut_name]["all"].Add(hists[_cut_name][_rawid])

# write the hists to file
# Changes the names
out_file = uproot.recreate(args.output)
for _cut_name, _hist_dict in hists.items():
    dir = out_file.mkdir(_cut_name)
    """for _rawid, _name in sorted(usable_geds.items()):
        dir[_rawid] = _hist_dict[_rawid]
    for _type in ["bege", "coax", "icpc", "ppc"]:
        dir[_type] = _hist_dict[_type]
    dir['all'] = _hist_dict['all']"""
    for key, item in _hist_dict.items():
        dir[key] = item
out_file["number_of_primaries"] = str(n_primaries_total)
out_file.close()

print(time.time() - start)
