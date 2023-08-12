#!/bin/bash

# module load python
# mamba activate snakemake

version=$(basename "$PWD")

# rm -rf .snakemake

logdir=".slurm/$(date +'%Y%m%dT%H%M%SZ')"
mkdir -p "$logdir"

if [[ "$1" == "parallel" ]]; then
    simids=$(python -c '
    import json

    with open("inputs/simprod/config/tier/raw/l200a/simconfig.json") as f:
        simids = json.load(f).keys()

    for s in simids:
        print(f"evt.{s}", end=" ")
    ')

    for s in $simids; do
        job="${version}_$s"
        echo "INFO: inspecting $job"

        if squeue --me --format '%200j' | grep "$job"; then
            echo "INFO: job already queued, skipping"
            continue
        fi

        snakemake --config simlist="$s" --dry-run | grep 'Nothing to be done' && continue

        echo "INFO: submitting..."
        # https://docs.nersc.gov/development/shifter/faq-troubleshooting/#failed-to-lookup-image
        sbatch \
            --nodes 1 \
            --ntasks-per-node=1 \
            --account m2676 \
            --constraint cpu \
            --time 12:00:00 \
            --qos regular \
            --licenses scratch,cfs \
            --job-name "$job" \
            --output "$logdir/$s.log" \
            --error "$logdir/$s.log" \
            --image "legendexp/legend-base:latest" \
            --wrap "
                srun snakemake \
                    --shadow-prefix $PSCRATCH \
                    --config simlist=$s
            "
    done
else
    job="${version}-legend-pdfs"

    if squeue --me --format '%200j' | grep "$job"; then
        echo "INFO: job already queued"
        exit 1
    fi

    snakemake --dry-run "$@" | grep 'Nothing to be done' && exit 1

    echo "INFO: submitting..."
    # https://docs.nersc.gov/development/shifter/faq-troubleshooting/#failed-to-lookup-image
    sbatch \
        --nodes 1 \
        --ntasks-per-node=1 \
        --account m2676 \
        --constraint cpu \
        --time 12:00:00 \
        --qos regular \
        --licenses scratch,cfs \
        --job-name "$job" \
        --output "$logdir/$job.log" \
        --error "$logdir/$job.log" \
        --image "legendexp/legend-base:latest" \
        --wrap "
            srun snakemake \
                --shadow-prefix $PSCRATCH \
                $*
        "
fi
