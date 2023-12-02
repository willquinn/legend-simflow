#!/bin/bash

doit() {
    echo "INFO: spawning job $1"
    shifter --image=legendexp/legend-software:latest -- \
        ./workflow/scripts/split_K42_vertices.C \
	generated/tier/raw/l200a-lar-cylinder-K42/l200a-lar-cylinder-K42_$1-tier_raw.root \
        generated/tier/raw/l200a-lar-inside-ms-K42/l200a-lar-inside-ms-K42_$1-tier_raw.root \
	generated/tier/raw/l200a-lar-outside-ms-K42/l200a-lar-outside-ms-K42_$1-tier_raw.root \
	&> /dev/null
}
export -f doit

nfiles=$(ls generated/tier/raw/l200a-lar-cylinder-K42/ | wc -l)
seq -f "%04g" 0 $(($nfiles-1)) | parallel doit {}
