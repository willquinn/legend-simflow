# ruff: noqa: F821, T201

import csv
from pathlib import Path


def printline(*line):
    print("{:<50}{}".format(*line))


printline("simid", "CPU time [s]")
printline("-----", "------------")

bdir = Path(snakemake.params.setup["paths"]["benchmarks"])

for simd in bdir.glob("*/*"):
    data = {"cpu_time": 0}
    for jobd in simd.glob("*.tsv"):
        with jobd.open(newline="") as f:
            this_data = list(csv.DictReader(f, delimiter="\t"))[0]
            data["cpu_time"] += float(this_data["cpu_time"])

    printline(simd.parent.name + "." + simd.name, "{:.5g}".format(data["cpu_time"]))
