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

# determine whether external vertices are required
if "vertices" in config:
    # get list of ver output files
    outver_list = simjobs.gen_list_of_simid_outputs(
        snakemake.params.setup, "ver", config["vertices"]
    )

    # if number of jobs is not specified, use number of vertices files
    if "number_of_jobs" not in config:
        n_macros = len(outver_list)

# if n_macros could not be determined, there MUST be an explicit reference to
# it in the config
if not n_macros:
    n_macros = config["number_of_jobs"]

# prepare global substitution rules
substitutions = {"NUMBER_OF_PRIMARIES": int(n_prim / n_macros)}

# first substitute global variables
with Path(snakemake.input.template).open() as f:
    text = utils.subst_vars(f.read().strip(), substitutions, ignore_missing=True)

# then substitute macro-specific variables
for i in range(n_macros):
    # determine output file name for this macro
    outname = simjobs.output_simjob_filename(
        snakemake.params.setup,
        snakemake.params.tier,
        snakemake.params.simid,
        i,
    )

    substitutions.update({"OUTPUT_FILE": outname})
    # might also need to set a vertices file name (see above)
    if outver_list:
        substitutions.update({"VERTICES_FILE": outver_list[i]})

    # determine the macro file name for write out
    inname = simjobs.input_simjob_filename(
        snakemake.params.setup, snakemake.params.tier, snakemake.params.simid, i
    )

    text_out = utils.subst_vars(text, substitutions).strip()

    # check if the file exists and open it
    # NOTE: unfortunately, this doesn't produce the desired effect with
    # Snakemake, as it always removes output files before regenerating them:
    # https://stackoverflow.com/questions/49178033/snakemake-avoid-removing-output-files-before-executing-the-shell-command
    if inname.is_file():
        with inname.open() as f:
            # if file content is already the expected one, do not overwrite it
            if f.read().strip() == text_out:
                continue

    # otherwise, prepare for writing and write
    smk.utils.makedirs(str(inname.parent))
    with inname.open("w") as f:
        f.write(text_out)
