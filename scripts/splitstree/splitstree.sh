#!/bin/bash

set -euo pipefail

## Software
SPLITSTREE=software/splitstree4/SplitsTree

## Command-line args:
infile=$1
outfile=$2

## Report:
echo "## Starting script splitstree.sh"
date
echo
echo "## Nexus input:      $infile"
echo "## Nexus output:     $outfile"
echo

## Run splitstree
echo "## Running Splitstree..."
$SPLITSTREE -g -i "$infile" -x "UPDATE; SAVE REPLACE=yes FILE=$outfile.tmp; QUIT"

## Remove sequences from nexus output
echo -e "\n## Removing actual sequence from Nexus output..."
START=$(grep -n "BEGIN Characters;" "$outfile".tmp | cut -f1 -d:)
END=$(grep -n "END;.*Characters" "$outfile".tmp | cut -f1 -d:)
sed "$START,${END}d" "$outfile".tmp > "$outfile"

rm -f "$outfile".tmp

echo -e "\n## Done with script splitstree.sh"
date
echo
