from __future__ import annotations

import json
from pathlib import Path

from . import patterns


def get_simid_n_macros(setup, tier, simid):
    """Returns the number of macros that will be generated for a given `tier`
    and `simid`."""

    if "benchmark" in setup and setup["benchmark"].get("enabled", False):
        return 1

    tdir = patterns.template_macro_dir(setup, tier=tier)

    with (Path(tdir) / "simconfig.json").open() as f:
        config = json.load(f)[simid]

    if "vertices" in config and "number_of_jobs" not in config:
        return len(gen_list_of_simid_outputs(setup, "ver", config["vertices"]))
    elif "number_of_jobs" in config:
        return config["number_of_jobs"]
    else:
        msg = "simulation config must contain 'vertices' or 'number_of_jobs'"
        raise RuntimeError(msg)


def collect_simconfigs(setup, tiers):
    cfgs = []
    for tier in tiers:
        with (
            patterns.template_macro_dir(setup, tier=tier) / "simconfig.json"
        ).open() as f:
            for sid, _val in json.load(f).items():
                cfgs.append((tier, sid, get_simid_n_macros(setup, tier, sid)))

    return cfgs


def gen_list_of_simid_inputs(setup, tier, simid):
    """Generates the full list of input files for a `tier` and `simid`."""
    n_macros = get_simid_n_macros(setup, tier, simid)
    return patterns.input_simjob_filenames(setup, n_macros, tier=tier, simid=simid)


def gen_list_of_simid_outputs(setup, tier, simid, max_files=None):
    """Generates the full list of output files for a `simid`."""
    n_macros = get_simid_n_macros(setup, tier, simid)
    if max_files is not None:
        n_macros = min(n_macros, max_files)
    return patterns.output_simjob_filenames(setup, n_macros, tier=tier, simid=simid)


def gen_list_of_all_simids(setup, tier):
    with (patterns.template_macro_dir(setup, tier=tier) / "simconfig.json").open() as f:
        return json.load(f).keys()


def gen_list_of_all_macros(setup, tier):
    mlist = []
    for sid in gen_list_of_all_simids(setup, tier):
        mlist += gen_list_of_simid_inputs(setup, tier=tier, simid=sid)

    return mlist


def gen_list_of_all_simid_outputs(setup, tier):
    mlist = []
    for sid in gen_list_of_all_simids(setup, tier):
        mlist += gen_list_of_simid_outputs(setup, tier=tier, simid=sid)

    return mlist


def gen_list_of_all_plt_outputs(setup, tier):
    mlist = []
    for sid in gen_list_of_all_simids(setup, tier):
        mlist += [
            patterns.plt_file_path(setup, tier=tier, simid=sid)
            + "/mage-event-vertices.png"
        ]

    return mlist


def gen_list_of_hit_outputs(setup, simid, max_files=None):
    """Generates the full list of output files for a `simid`, 'hit' tier."""
    n_macros = get_simid_n_macros(setup, "raw", simid)
    if max_files is not None:
        n_macros = min(n_macros, max_files)
    return patterns.output_hit_filenames(setup, n_macros, simid=simid)


def gen_list_of_all_hit_outputs(setup):
    mlist = []
    for sid in gen_list_of_all_simids(setup, "raw"):
        mlist += gen_list_of_hit_outputs(setup, simid=sid)

    return mlist
