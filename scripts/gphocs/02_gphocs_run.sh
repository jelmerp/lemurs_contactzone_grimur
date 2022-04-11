#!/bin/bash

set -euo pipefail

# SET-UP -----------------------------------------------------------------------
## Software
GPHOCS=software/G-PhoCS/bin/G-PhoCS

## Command-line args
CFILE=$1        # Path to control file
NCORES=$2       # Number of cores to use

## Report
echo "## Starting script 02_gphocs_run.sh"
date
echo
echo "## Control file: $CFILE"
echo "## Number of processors: $NCORES"
echo

## Copying old logfile:
LOGFILE=$(grep "trace-file" $CFILE | sed 's/trace-file //')
echo -e "\n## Logfile: $LOGFILE"

if [ -f "$LOGFILE" ]; then
	echo -e "## Logfile exists: moving to $LOGFILE.autocopy...\n"
	mv "$LOGFILE" "$LOGFILE".autocopy.log
fi

## Editing controlfile to include date:
echo "## Editing controlfile..."
DATE=$(date +%Y%m%d-%H%M)
cp "$CFILE" "$CFILE".tmp
sed -e "s/\(date[0-9][0-9][0-9][0-9][0-9][0-9]\).*log/\1.$DATE.log/" "$CFILE".tmp > "$CFILE" 
rm "$CFILE".tmp

echo -e "\n## Showing top of controlfile:"
head "$CFILE"


# RUN GPHOCS -------------------------------------------------------------------
echo -e "\n\n## Starting Gphocs run...\n"

export OMP_NUM_THREADS=$NCORES
$GPHOCS "$CFILE" -n "$NCORES"


## Report:
echo -e "## Done with script 02_gphocs_run.sh"
date
echo
