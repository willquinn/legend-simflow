"""Prepare pattern strings to be used in Snakemake rules.

Keyword arguments are typically interpreted as variables to be substituted in
the returned (structure of) strings. They are passed to
:func:`snakemake.io.expand`.
"""

import json
from pathlib import Path
from snakemake.io import expand


def simjob_rel_basename():
    """Formats a partial output path for a `simid` and `jobid`."""
    return "{simid}/{simid}_{jobid}"


# TODO: should maybe create a random new directory for each snakemake run?
def genmacro_log_file_path(setup):
    """Formats a log file path for a `simid` and `jobid`."""
    return str(
        Path(setup["paths"]["log"])
        / "macros"
        / "{tier}"
        / (simjob_rel_basename() + ".log")
    )


def template_macro_dir(setup, **kwargs):
    """Returns the directory path to the macro templates for the current `tier`."""
    tier = expand("{tier}", **kwargs, allow_missing=True)[0]
    return Path(setup["paths"]["config"]) / "tier" / tier / setup["experiment"]


def macro_gen_inputs(setup, tier, simid, **kwargs):
    """Return inputs for `generate_macros` Snakemake rule."""
    tdir = template_macro_dir(setup, tier=tier)

    with (tdir / "simconfig.json").open() as f:
        config = json.load(f)[simid]

    if "template" not in config:
        raise RuntimeError("simconfig.json blocks must define a 'template' field.")

    expr = {
        "template": str(tdir / config["template"]),
        "cfgfile": str(tdir / "simconfig.json"),
    }
    for k, v in expr.items():
        expr[k] = expand(v, **kwargs, allow_missing=True)[0]
    return expr


def input_simjob_filename(setup, tier):
    """Returns the full path to the input file for a `simid` and job index."""
    return str(
        Path(setup["paths"]["macros"])
        / f"{tier}"
        / (simjob_rel_basename() + setup["filetypes"]["input"][tier])
    )


def output_simjob_filename(setup, tier, **kwargs):
    """Returns the full path to the output file for a `simid`, `tier` and job index."""
    expr = str(
        Path(setup["paths"][f"tier_{tier}"])
        / (simjob_rel_basename() + setup["filetypes"]["output"][tier])
    )
    return expand(expr, **kwargs, allow_missing=True)[0]


def input_simjob_filenames(setup, n_macros, **kwargs):
    """Returns the full path to all input files for a `simid`."""
    tier = kwargs.get("tier", None)

    if tier is None:
        raise RuntimeError("the 'tier' argument is mandatory")

    pat = str(
        Path(setup["paths"]["macros"])
        / "{tier}"
        / (simjob_rel_basename() + setup["filetypes"]["input"][tier])
    )
    jobids = expand("{id:>04d}", id=list(range(n_macros)))
    return expand(pat, jobid=jobids, **kwargs, allow_missing=True)


def log_file_path(setup, tier):
    """Formats a log file path for a `simid` and `jobid`."""
    return str(
        Path(setup["paths"]["log"]) / f"{tier}" / (simjob_rel_basename() + ".log")
    )


def benchmark_file_path(setup, tier):
    """Formats a benchmark file path for a `simid` and `jobid`."""
    return str(
        Path(setup["paths"]["benchmarks"])
        / f"{tier}"
        / (simjob_rel_basename() + ".tsv")
    )


def run_command(setup, tier):
    """Returns command to build files in tier `tier` prefixed by environment."""
    return "{swenv} " + setup["runcmd"][tier]


def smk_ver_filename_for_raw(setup, wildcards):
    """Returns the vertices file needed for the 'raw' tier job, if needed."""
    tdir = template_macro_dir(setup, tier="raw")

    with (tdir / "simconfig.json").open() as f:
        config = json.load(f)[wildcards.simid]

    if "vertices" in config:
        return output_simjob_filename(setup, "ver", simid=config["vertices"])
    else:
        return []
