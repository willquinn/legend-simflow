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
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

# ruff: noqa: F821, T201

from __future__ import annotations

import argparse

parser = argparse.ArgumentParser(
    prog="build_pdf", description="build LEGEND pdf files from evt tier files"
)
parser.add_argument("--output", "-o", help="output file name")
parser.add_argument("--metadata", "-m", help="path to legend-metadata")
parser.add_argument("input_files", nargs="+", help="evt tier files")

args = parser.parse_args()

if not isinstance(args.input_files, list):
    args.input_files = [args.input_files]

# TODO: move build_pdf script here
