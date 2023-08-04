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

import ROOT
import uproot
from legendmeta import LegendMetadata


def process_mage_id(mage_id):
    m_id = str(mage_id)
    is_ged = bool(int(m_id[0]))
    if not is_ged:
        return False

    string = int(m_id[3:5])
    pos = int(m_id[5:7])

    for key, value in chmap.items():
        if isinstance(value, dict) and "location" in value:
            location = value["location"]
            usable = value["analysis"]["usability"] == "on"
            if (
                location.get("string") == string
                and location.get("position") == pos
                and usable
            ):
                return {"name": key, "ch": value["daq"]["rawid"], "mage_id": mage_id}

    return False


parser = argparse.ArgumentParser(
    prog="build_pdf", description="build LEGEND pdf files from evt tier files"
)
parser.add_argument("--config", "-c", help="configuration file")
parser.add_argument("--output", "-o", help="output file name")
parser.add_argument("--metadata", "-m", help="path to legend-metadata")
parser.add_argument("input_files", nargs="+", help="evt tier files")

args = parser.parse_args()

meta = LegendMetadata(args.metadata)
chmap = meta.channelmap("20230323T000000Z")

if not isinstance(args.input_files, list):
    args.input_files = [args.input_files]

with Path(args.config).open() as f:
    rconfig = json.load(f)

for file_name in args.input_files:
    # Get the tree from the file
    file = ROOT.TFile(file_name, "READ")
    tree = file.Get("simTree")

    uniq_mage_ids = list(
        {mage_id for n in [list(tree.mage_id) for entry in tree] for mage_id in n}
    )
    mage_names = {
        mage_id: process_mage_id(mage_id)
        for mage_id in uniq_mage_ids
        if process_mage_id(mage_id)
    }

    # Define a dictionary to contain all the histograms to be filled
    # for each cut and each detector
    hists = {
        f"{cut_name}": {
            _mage_id: ROOT.TH1F(
                f"{cut_name}_{_mage_names['ch']}",
                f"{_mage_names['name']} energy deposition",
                rconfig["hist"]["nbins"],
                rconfig["hist"]["emin"],
                rconfig["hist"]["emax"],
            )
            for _mage_id, _mage_names in mage_names.items()
        }
        for cut_name in rconfig["cuts"]
    }
    # Add the grouped dets hists
    for cut_name in rconfig["cuts"]:
        for _type in ["bege", "coax", "icpc", "ppc"]:
            hists[cut_name][f"{_type}"] = ROOT.TH1F(
                f"{cut_name}_type_{_type}",
                f"All {_type} energy deposit",
                rconfig["hist"]["nbins"],
                rconfig["hist"]["emin"],
                rconfig["hist"]["emax"],
            )
        hists[cut_name]["all"] = ROOT.TH1F(
            f"{cut_name}_type_all",
            "All energy deposit",
            rconfig["hist"]["nbins"],
            rconfig["hist"]["emin"],
            rconfig["hist"]["emax"],
        )

    out_file_name = file_name.split("/")[-1].split(".")[0] + "_pdf" + ".root"
    out_file = uproot.recreate(args.output + out_file_name)

    # Loop over the data ONCE and fill as we go along
    for entry in tree:
        # In the real data we have an energy threshold defined in the config
        energy_list = [
            val for val in entry.energy if val >= rconfig["energy_threshold"]
        ]
        mage_id_list = [
            val
            for index, val in enumerate(entry.mage_id)
            if entry.energy[index] >= rconfig["energy_threshold"]
        ]
        npe_tot = entry.npe_tot
        if len(energy_list) == 0:
            continue  # Nothing to fill

        for cut_name, cut_string in rconfig["cuts"].items():
            # In the config is a lambda function string that returns true or false for each cut
            # Maybe this is too much but it allows for more cuts possibly
            exec(cut_string)
            pass_cut = func(energy_list, mage_id_list, npe_tot)
            if not bool(pass_cut):
                continue

            for mage_id, energy in zip(mage_id_list, energy_list):
                energy_kev = energy * 1000
                mage_dict = process_mage_id(mage_id)

                if mage_dict is False:
                    continue  # not an 'on' detector so we don't fill
                hists[cut_name][mage_id].Fill(energy_kev)
                hists[cut_name][chmap[mage_dict["name"]]["type"]].Fill(energy_kev)
                hists[cut_name]["all"].Fill(energy_kev)

    # Create a new root file with the directories
    for cut_name in rconfig["cuts"]:
        dir = out_file.mkdir(cut_name)
        for _key, item in hists[cut_name].items():
            dir[item.GetName()] = item

    out_file.close()
    file.Close()
