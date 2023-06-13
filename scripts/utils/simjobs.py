import json
from pathlib import Path


def simjob_rel_basename(simid, jobid) -> str:
    """Formats a partial output path for a `simid` and `jobid`."""
    return f"{simid}/{simid}_{jobid:>04d}"


def template_macro_dir(setup, tier) -> str:
    """Returns the directory path to the macro templates for `tier`."""
    return Path(setup["paths"]["config"]) / "tier" / tier / setup["experiment"]


def input_simjob_filename(setup, tier, simid, idx) -> str:
    """Returns the full path to the input file for a `tier`,
    `simid` and job index `idx`."""
    destdir = Path(setup["paths"]["macros"]) / tier

    return destdir / (
        simjob_rel_basename(simid, idx) + setup["filetypes"]["input"][tier]
    )


def output_simjob_filename(setup, tier, simid, idx):
    """Returns the full path to the output file for a `tier`, `simid`,
    `tier` and job index `idx`."""
    destdir = Path(setup["paths"][f"tier_{tier}"])

    return destdir / (
        simjob_rel_basename(simid, idx) + setup["filetypes"]["output"][tier]
    )


def get_simid_n_macros(setup, tier, simid):
    """Returns the number of macros that will be generated for a given `tier`
    and `simid`."""
    tdir = template_macro_dir(setup, tier)

    with (Path(tdir) / "simconfig.json").open() as f:
        config = json.load(f)[simid]

    if "vertices" in config and "number_of_jobs" not in config:
        return len(gen_list_of_simid_outputs(setup, "ver", config["vertices"]))
    elif "number_of_jobs" in config:
        return config["number_of_jobs"]
    else:
        raise RuntimeError(
            "simulation config must contain 'vertices' or 'number_of_jobs'"
        )


def gen_list_of_simid_inputs(setup, tier, simid):
    """Generates the full list of input files for a `tier` and `simid`."""
    n_macros = get_simid_n_macros(setup, tier, simid)
    return [str(input_simjob_filename(setup, tier, simid, j)) for j in range(n_macros)]


def gen_list_of_simid_outputs(setup, tier, simid):
    """Generates the full list of output files for a `simid`."""
    n_macros = get_simid_n_macros(setup, tier, simid)
    return [str(output_simjob_filename(setup, tier, simid, j)) for j in range(n_macros)]
