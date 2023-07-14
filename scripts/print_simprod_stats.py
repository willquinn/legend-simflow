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
# along with this program.  If not, see <https:#www.gnu.org/licenses/>.

# ruff: noqa: F821, T201

import csv
import json
from datetime import timedelta
from pathlib import Path

from utils import patterns, simjobs


def printline(*line):
    print("{:<50} {:>20} {:>6} {:>14} {:>10}".format(*line))


printline("     ", "wall time [s]", "    ", "wall time [s]", "         ")
printline("simid", " (cumulative)", "jobs", "    (per job)", "primaries")
printline("-----", "-------------", "----", "-------------", "---------")

bdir = Path(snakemake.params.setup["paths"]["benchmarks"])

for simd in sorted(bdir.glob("*/*")):
    data = {"wall_time": 0}
    for jobd in simd.glob("*.tsv"):
        with jobd.open(newline="") as f:
            this_data = list(csv.DictReader(f, delimiter="\t"))[0]
            data["wall_time"] += float(this_data["s"])

    tdir = patterns.template_macro_dir(snakemake.params.setup, tier=simd.parent.name)
    with (tdir / "simconfig.json").open() as f:
        config = json.load(f)[simd.name]

    nprim = config["number_of_primaries"]
    njobs = simjobs.get_simid_n_macros(
        snakemake.params.setup, simd.parent.name, simd.name
    )

    printline(
        simd.parent.name + "." + simd.name,
        str(timedelta(seconds=int(data["wall_time"]))),
        njobs,
        str(timedelta(seconds=int(data["wall_time"])) / njobs),
        f"{nprim:.2E}",
    )
