# ruff: noqa: F821

import csv
import json
from pathlib import Path
from datetime import timedelta

from utils import patterns, simjobs


def printline(*line):
    print("{:<50} {:>20} {:>6} {:>14} {:>10}".format(*line))


printline("     ", "wall time [s]", "    ", "wall time [s]", "         ")
printline("simid", " (cumulative)", "jobs", "    (per job)", "primaries")
printline("-----", "-------------", "----", "-------------", "---------")

bdir = Path(snakemake.params.setup["paths"]["benchmarks"])

for simd in bdir.glob("*/*"):
    data = {"wall_time": 0}
    for jobd in simd.glob("*.tsv"):
        with jobd.open(newline="") as f:
            this_data = list(csv.DictReader(f, delimiter="\t"))[0]
            data["wall_time"] += float(this_data["s"])

    tdir = patterns.template_macro_dir(snakemake.params.setup, tier=simd.parent.name)
    with (tdir / "simconfig.json").open() as f:
        config = json.load(f)[simd.name]

    nprim = config["number_of_primaries"]
    njobs = simjobs.get_simid_n_macros(snakemake.params.setup, simd.parent.name, simd.name)

    printline(
        simd.parent.name + "." + simd.name,
        str(timedelta(seconds=int(data["wall_time"]))),
        njobs,
        str(timedelta(seconds=int(data["wall_time"])) / njobs),
        "{:.2E}".format(nprim),
    )
