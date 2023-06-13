# ruff: noqa: F821

import snakemake as smk
import json
from pathlib import Path

from utils import utils, simjobs


with open(snakemake.input.cfgfile) as f:
    config = json.load(f)[snakemake.params.simid]

n_prim = config["number_of_primaries"]
n_macros = None
outver_list = None

if "vertices" in config:
    # get list of ver output files
    outver_list = simjobs.gen_list_of_simid_outputs(
        snakemake.params.setup, "ver", config["vertices"]
    )

    # if number of jobs is not specified, use number of vertices files
    if "number_of_jobs" not in config:
        n_macros = len(outver_list)

if not n_macros:
    n_macros = config["number_of_jobs"]

substitutions = {"NUMBER_OF_PRIMARIES": int(n_prim / n_macros)}

# first substitute global variables
with Path(snakemake.input.template).open() as f:
    text = utils.subst_vars(f.read(), substitutions, ignore_missing=True)

# then substitute macro-specific variables
for i in range(n_macros):
    outname = simjobs.output_simjob_filename(
        snakemake.params.setup,
        snakemake.params.tier,
        snakemake.params.simid,
        i,
    )

    substitutions.update({"OUTPUT_FILE": outname})
    if outver_list:
        substitutions.update({"VERTICES_FILE": outver_list[i]})

    inname = simjobs.input_simjob_filename(
        snakemake.params.setup, snakemake.params.tier, snakemake.params.simid, i
    )

    smk.utils.makedirs(str(inname.parent))

    with inname.open("w") as f:
        f.write(utils.subst_vars(text, substitutions))
