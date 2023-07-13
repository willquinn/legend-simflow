from __future__ import annotations

import argparse

import ROOT

ROOT.PyConfig.IgnoreCommandLineOptions = False


parser = argparse.ArgumentParser(
    prog="plot_mage_vertices", description="plot MaGe primary event vertices"
)
parser.add_argument("--output", "-o", help="output file name")
parser.add_argument("input_files", nargs="+", help="MaGe output files")

args = parser.parse_args()

if not isinstance(args.input_files, list):
    args.input_files = [args.input_files]

# do not display any of the standard histogram decorations
ROOT.gStyle.SetOptStat(False)
ROOT.gStyle.SetOptFit(0)

tree = ROOT.TChain("fTree")
for f in args.input_files:
    tree.Add(f)

c = ROOT.TCanvas("c", "MaGe event vertices", 1000, 1000)
c.Divide(2, 2)
c.cd(1)
tree.Draw("eventPrimaries.fSteps.fY:eventPrimaries.fSteps.fX")
c.cd(2)
tree.Draw("eventPrimaries.fSteps.fZ:eventPrimaries.fSteps.fY")
c.cd(3)
tree.Draw("eventPrimaries.fSteps.fZ:eventPrimaries.fSteps.fX")
c.cd(4)
tree.Draw("eventPrimaries.fSteps.fZ:eventPrimaries.fSteps.fY:eventPrimaries.fSteps.fX")

if args.output:
    c.SaveAs(args.output)
