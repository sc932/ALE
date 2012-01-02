#!/bin/sh
USAGE="$0 input.ale"
in=$1
if [ ! -f "$in" ]
then
  echo "$USAGE"
  exit 1
fi

awk '!/#/ { sum += $7; count++;} END {print sum / count }' $in
