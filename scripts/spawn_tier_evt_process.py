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

import json
import subprocess
from collections import OrderedDict
from pathlib import Path

import uproot
from snakemake.io import expand
from utils import patterns

# open file with run livetime partitioning
with Path(snakemake.input.run_part_file[0]).open() as f:
    runpart = json.load(f, object_pairs_hook=OrderedDict)

# get total number of mc events
file_evts = uproot.num_entries(
    [f"{file}:simTree" for file in snakemake.input.hit_files]
)
tot_events = 0
for file in file_evts:
    tot_events += file[-1]

runs = list(runpart.keys())
weights = [tot_events * v for v in runpart.values()]

# compute start event and number of events for this run
start_event = sum(weights[: runs.index(snakemake.wildcards.runid)])
n_events = tot_events * runpart[snakemake.wildcards.runid]

# substitute $START_EVENT and $N_EVENTS in the command line
cmd = expand(
    patterns.run_command(snakemake.config, tier="evt"),
    _start_event=start_event,
    _n_events=n_events,
    allow_missing=True,
)[0]

# exec command
subprocess.run(cmd)
