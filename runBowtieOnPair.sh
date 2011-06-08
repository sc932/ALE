#!/bin/bash

reference=$1
pairedFastq=$2
cores=$3
dir1=$4
ins1=$5
std1=$6
stdWidth=$7
dir2=$8
ins2=$9
std2=${10}

out=${reference##*/}-${2##*/}

USAGE="Proper Usage:
$0 reference.fa paired.fq numCores dir1 ins1 std1 [stdWidth [ dir2 ins2 std2 ]]

"

if [ $# -lt 6 ] || [ ! -f $reference ] || [ ! -f $pairedFastq ]
then
  echo "$USAGE" 1>&2
  exit 1
fi

[ -n "$stdWidth" ] || stdWidth=10

if [ ! -d "$TMPDIR" ]
then
  TMPDIR=/tmp
fi
export TMPDIR

set -e
not1=$TMPDIR/not1.$$
not2=$TMPDIR/not2.$$

CLEAN="$not1 $not2"
cleanup()
{
  state=$?
  trap '' 0 1 2 3 15
  if [ $state -ne 0 ]
  then
    echo "ERROR" 1>&1
  fi
  rm -rf $CLEAN
  exit $?
}
trap cleanup 0 1 2 3 15

set -x 
base=$(dirname $(which $0))

bowtieOpts="-p $cores -a -l 16 -v 3 -S -t -r --sam-nohead --phred64-quals"
$base/illumina2crossbow.pl $pairedFastq \
    | bowtie $bowtieOpts --un $not1 \
             --$dir1 --minins $((ins1-std1*stdWidth)) --maxins $((ins1+std1*stdWidth)) \
             $reference --12 - \
    | awk 'and($2,0x04) != 0x04 {print}' \
    | gzip -c > $out.paired.sam.gz

if [ -n "$dir2" ]
then
    bowtie $bowtieOpts --un $not2 \
             --$dir2  --minins $((ins2-std2*stdWidth)) --maxins $((ins2+std2*stdWidth)) \
             $reference --12 $not1 \
    | awk 'and($2,0x04) != 0x04 {print}' \
    | gzip -c >> $out.paired.sam.gz 
else
  rm $not2
  ln $not1 $not2
fi

if [ -s $not2 ]
then
    bowtie $bowtieOpts --un $out.unmapped $reference --12 $not2 \
        | awk 'and($2,0x04) != 0x04 {print}' \
        | gzip -c > $out.single.sam.gz 
fi

