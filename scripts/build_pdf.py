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
import ROOT, uproot
from legendmeta import LegendMetadata

parser = argparse.ArgumentParser(
    prog="build_pdf", description="build LEGEND pdf files from evt tier files"
)
parser.add_argument("--output", "-o", help="output file name")
parser.add_argument("--metadata", "-m", help="path to legend-metadata")
parser.add_argument("input_files", nargs="+", help="evt tier files")

args = parser.parse_args()

meta = LegendMetadata(args.metadata)
chmap = meta.channelmap("20230323T000000Z")

if not isinstance(args.input_files, list):
    args.input_files = [args.input_files]

# TODO: move build_pdf script here

def process_mage_id(mage_id, chmap):
    m_id = str(mage_id)
    is_ged = bool(int(m_id[0]))
    if not is_ged:
        return False

    string = int(m_id[3:5])
    pos = int(m_id[5:7])

    for key, value in chmap.items():
        if isinstance(value, dict) and 'location' in value:
            location = value['location']
            if location.get('string') == string and location.get('position') == pos:
                return {'name': key, 'ch': value['daq']['rawid'], 'mage_id': mage_id}

    return False

# Get the tree from the file 
file = ROOT.TFile(args.input_files[0], "READ")
tree = file.Get('simTree')
# Can also access wihtout tree name
# Tobject = file.GetListOfKeys().At(0).ReadObj() # This assumes that there is only one object in the file

uniq_mage_ids = list(set([mage_id for n in [list(tree.mage_id) for entry in tree] for mage_id in n]))
mage_names = {mage_id : process_mage_id(mage_id, chmap) for mage_id in uniq_mage_ids if process_mage_id(mage_id, chmap)}
    
# Set up the root output file
#if args.output_file == './':
#    output_file = args.input_file.split('.')[0] + '_pdf.root'
        
out_file = uproot.create(args.output)

# list of cuts
cuts = {
    'raw': '',
    'mul': 'mage_id@.size()==1',
    'lar': 'npe_tot==0',
    'mul_lar': 'mage_id@.size()==1 && npe_tot==0'
}

for cut_name, cut_string in cuts.items():
    dir = out_file.mkdir(cut_name)

    # Define the grouped histos to be updated 
    grouped = {
        f'{i}': ROOT.TH1D(
                f'hist_{cut_name}_{i}',
                f'{i} energy deposit',
                3000, 0, 3
            ) for i in ['bege', 'coax', 'icpc', 'ppc']
    }

    # loop over the channels
    for mage_id, mg_dict in mage_names.items():
        if cut_string == '': 
            string = f'mage_id=={mage_id}'
        else:
            string = f'mage_id=={mage_id} && {cut_string}'
            
        hist = ROOT.TH1D(f'hist_{cut_name}_{mg_dict["ch"]}', f'{mg_dict["name"]} energy deposit', 3000, 0, 3) # units of MeV and bin width of 1keV
        tree.Project(hist.GetName(), 'energy', string) # add cut
        dir[hist.GetName()] = hist

        tree.Project(grouped[chmap[mg_dict["name"]]['type']].GetName(), 'energy', string)

    for ged_type, hist in grouped.items():
       dir[hist.GetName()] = hist 

    # Then finally ALL
    hist = ROOT.TH1D(f'hist_{cut_name}_all', 'All channels energy deposit', 3000, 0, 3) # units of MeV and bin width of 1keV
    tree.Project(hist.GetName(), 'energy', cut_string) # add cut
    dir[hist.GetName()] = hist

out_file.close()
file.Close()