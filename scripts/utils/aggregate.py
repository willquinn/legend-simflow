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

from __future__ import annotations

import json
from pathlib import Path

from . import patterns


def get_simid_n_macros(setup, tier, simid):
    """Returns the number of macros that will be generated for a given `tier`
    and `simid`."""
    if tier not in ("ver", "raw"):
        tier = "raw"

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
    return patterns.input_simid_filenames(setup, n_macros, tier=tier, simid=simid)


def gen_list_of_simid_outputs(setup, tier, simid, max_files=None):
    """Generates the full list of output files for a `simid`."""
    n_macros = get_simid_n_macros(setup, tier, simid)
    if max_files is not None:
        n_macros = min(n_macros, max_files)
    return patterns.output_simid_filenames(setup, n_macros, tier=tier, simid=simid)


def gen_list_of_all_simids(setup, tier):
    if tier not in ("ver", "raw"):
        tier = "raw"
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


def gen_list_of_all_plots_outputs(setup, tier):
    mlist = []
    for sid in gen_list_of_all_simids(setup, tier):
        mlist += [
            patterns.plots_file_path(setup, tier=tier, simid=sid)
            + "/mage-event-vertices.png"
        ]

    return mlist


def process_simlist_or_all(setup, simlist=None):
    if simlist is None:
        simlist = setup.get("simlist", None)

    if Path(simlist).is_file():
        with Path(simlist).open() as f:
            simlist = [l.rstrip() for l in f.readlines()]
    elif isinstance(simlist, str):
        simlist = [simlist]

    mlist = []
    for line in simlist:
        mlist += gen_list_of_simid_outputs(
            setup, tier=line.split(".")[0], simid=line.split(".")[1].rstrip()
        )

    return mlist
