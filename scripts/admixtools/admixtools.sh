#!/bin/bash

set -euo pipefail

# SET-UP -----------------------------------------------------------------------
## Software
ATOOLS_DSTAT=software/AdmixTools-master/bin/qpDstat
ATOOLS_F4RATIO=software/AdmixTools-master/bin/qpF4ratio

## Command-line arguments:
file_id_full=$1
popfile_line=$2
parfile=$3
output=$4
atools_mode=$5

## Report:
echo
echo "## Starting script admixtools.sh"
date
echo
echo "## File ID:                          $file_id_full"
echo "## Line of popfile to process:       $popfile_line"
echo "## Parfile:                          $parfile"
echo "## Output file:                      $output"
echo "## Admixtools mode:                  $atools_mode"
echo -e "----------------------------\n"


# RUN ADMIXTOOLS ---------------------------------------------------------------
if [ "$atools_mode" == D ]; then

	if [ "$popfile_line" == ALL ]; then

		echo "## Running all popfile popfile_lines at once..."
		echo "## Running qpDstat in D-mode..."
		$ATOOLS_DSTAT -p "$parfile" > "$output"
	
    else
		
        echo "## Running one popfile popfile_line: $popfile_line"
		echo "## Running qpDstat in D-mode..."
		LINEWISE_output=$output.line${popfile_line}
		
		$ATOOLS_DSTAT -l "$popfile_line" -h "$popfile_line" -p "$parfile" > "$LINEWISE_output"
		
		echo "## Output file: $LINEWISE_output"
		cat "$LINEWISE_output"

	fi
fi

if [ "$atools_mode" == "F4RATIO" ]; then
	
    echo "## Running f4-ratio test..."
	
    $ATOOLS_F4RATIO -p "$parfile" > "$output".raw
	grep "result" "$output".raw  > "$output"

fi


# WRAP UP ----------------------------------------------------------------------
echo -e "\n## Done with script admixtools.sh"
date
echo
