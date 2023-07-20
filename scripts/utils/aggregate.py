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

from . import patterns, utils


def get_simid_n_macros(config, tier, simid):
    """Returns the number of macros that will be generated for a given `tier`
    and `simid`."""
    if tier not in ("ver", "raw"):
        tier = "raw"

    if "benchmark" in config and config["benchmark"].get("enabled", False):
        return 1

    tdir = patterns.template_macro_dir(config, tier=tier)

    with (Path(tdir) / "simconfig.json").open() as f:
        sconfig = json.load(f)[simid]

    if "vertices" in sconfig and "number_of_jobs" not in sconfig:
        return len(gen_list_of_simid_outputs(config, "ver", sconfig["vertices"]))
    elif "number_of_jobs" in sconfig:
        return sconfig["number_of_jobs"]
    else:
        msg = "simulation config must contain 'vertices' or 'number_of_jobs'"
        raise RuntimeError(msg)


def collect_simconfigs(config, tiers):
    cfgs = []
    for tier in tiers:
        with (
            patterns.template_macro_dir(config, tier=tier) / "simconfig.json"
        ).open() as f:
            for sid, _val in json.load(f).items():
                cfgs.append((tier, sid, get_simid_n_macros(config, tier, sid)))

    return cfgs


def gen_list_of_simid_inputs(config, tier, simid):
    """Generates the full list of input files for a `tier` and `simid`."""
    n_macros = get_simid_n_macros(config, tier, simid)
    return patterns.input_simid_filenames(config, n_macros, tier=tier, simid=simid)


def gen_list_of_simid_outputs(config, tier, simid, max_files=None):
    """Generates the full list of output files for a `simid`."""
    n_macros = get_simid_n_macros(config, tier, simid)
    if max_files is not None:
        n_macros = min(n_macros, max_files)
    return patterns.output_simid_filenames(config, n_macros, tier=tier, simid=simid)


def gen_list_of_all_simids(config, tier):
    if tier not in ("ver", "raw"):
        tier = "raw"
    with (
        patterns.template_macro_dir(config, tier=tier) / "simconfig.json"
    ).open() as f:
        return json.load(f).keys()


def gen_list_of_all_macros(config, tier):
    mlist = []
    for sid in gen_list_of_all_simids(config, tier):
        mlist += gen_list_of_simid_inputs(config, tier=tier, simid=sid)

    return mlist


def gen_list_of_all_simid_outputs(config, tier, only_simid=None):
    mlist = []
    slist = gen_list_of_all_simids(config, tier) if only_simid is None else only_simid
    for sid in slist:
        mlist += gen_list_of_simid_outputs(config, tier=tier, simid=sid)

    return mlist


def gen_list_of_all_plots_outputs(config, tier):
    mlist = []
    for sid in gen_list_of_all_simids(config, tier):
        mlist += [
            patterns.plots_file_path(config, tier=tier, simid=sid)
            + "/mage-event-vertices.png"
        ]

    return mlist


def process_simlist_or_all(config, simlist=None):
    if simlist is None:
        simlist = utils.get_some_list(config["simlist"])

    mlist = []
    for line in simlist:
        mlist += gen_list_of_simid_outputs(
            config, tier=line.split(".")[0], simid=line.split(".")[1].rstrip()
        )

    return mlist
