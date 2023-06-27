#!/usr/bin/env bash

macro="$1"
logfile="$2"

MaGe --no-colors "$macro" &> "$logfile"

code=0

grep -i "command not found" "$logfile" && code=1
grep -B 10 -i "batch is interrupted" "$logfile" && code=1
grep -B 10 -P "^(Error|Fatal):\w+" "$logfile" && code=1
grep -C 100 -iP "segmentation (fault|violation)" "$logfile" && code=1

if [ "$code" != "0" ]; then
   echo -e "\n"
   echo -e " ┌──────────────────────────────────────────────────────────┐"
   echo -e " │ MaGe returned a non-zero code, see above or the log file │"
   echo -e " └──────────────────────────────────────────────────────────┘"
fi

exit $code
