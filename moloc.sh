#!/bin/bash -x

while read LINE; do
        qsub -v FILE=$LINE moloc.pbs
done < /scratch/msoliai/moloc/r.code/tagc.r.list.txt
