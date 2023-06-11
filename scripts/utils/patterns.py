"""Prepare pattern strings to be used in Snakemake rules.

Keyword arguments are typically interpreted as variables to be substituted in
the returned (structure of) strings. They are passed to
:func:`snakemake.io.expand`.
"""

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


def template_macro_dir(setup):
    """Returns the directory path to the macro templates for the current `tier`."""
    return Path(setup["paths"]["config"]) / "tier" / "{tier}" / setup["experiment"]


def macro_gen_inputs(setup, **kwargs):
    """Return inputs for `generate_macros` Snakemake rule."""
    expr = {
        "template": str(template_macro_dir(setup) / "{simid}.template.mage.mac"),
        "cfgfile": str(template_macro_dir(setup) / "simconfig.json"),
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


def output_simjob_filename(setup, tier):
    """Returns the full path to the output file for a `simid`, `tier` and job index."""
    return str(
        Path(setup["paths"][f"tier_{tier}"])
        / (simjob_rel_basename() + setup["filetypes"]["output"][tier])
    )


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


def run_command(setup, tier):
    """"""
    return "{swenv} " + setup["runcmd"][tier]
