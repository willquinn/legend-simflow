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
from pathlib import Path

# TODO: configure mpp call with snakemake.inputs/outputs

with Path(snakemake.input.run_part_file[0]).open() as f:
    json.load(f)

# TODO: open snakemake.input.hitfiles and calculate total amount of events
# (with PyROOT?), then combine with run_part_file info to determine the mpp
# command line parameters

with Path(snakemake.output[0]).open("w") as f:
    f.write("")
