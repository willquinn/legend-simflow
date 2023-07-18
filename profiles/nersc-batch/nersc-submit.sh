#!/bin/bash

# module load python
# mamba activate snakemake

cd "$1" || exit 1

logdir=".slurm/$(date +'%Y%m%dT%H%M%SZ')"
mkdir -p "$logdir"

simids=$(python -c '
import json

with open("inputs/simprod/config/tier/raw/l200a/simconfig.json") as f:
    simids = json.load(f).keys()

for s in simids:
    print(f"raw.{s}", end=" ")
')

for s in $simids; do
    sbatch \
        --nodes 1 \
        --ntasks-per-node=1 \
        --account m2676 \
        --constraint cpu \
        --time 12:00:00 \
        --qos regular \
        --licenses scratch,cfs \
        --job-name "$s" \
        --output "$logdir/$s.log" \
        --error "$logdir/$s.log" \
        --wrap "
            srun snakemake \
                --profile workflow/profiles/nersc-interactive \
                --shadow-prefix $PSCRATCH \
                --config simlist=$s
        "
done
