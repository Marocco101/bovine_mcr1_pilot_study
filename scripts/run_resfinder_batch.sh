#!/bin/bash
# ResFinder batch processing for PRJNA1166088 (130 isolates)
# Author: Makiko FUJITA-SUZANNE
# Date: April 6 2026

cd ~/amr_pilot

for dir in ~/Documents/bovine_mcr1_pilot_study/data/GCA_*/; do
    sample=$(basename "$dir")
    echo "Processing: $sample"
    python -m resfinder \
      -ifa "$dir"/*.fna \
      -o ~/amr_pilot/results/"$sample" \
      -s "Escherichia coli" \
      --acquired \
      -t 0.90 \
      -l 0.60 \
      -db_res ~/amr_pilot/databases/resfinder_db
done

echo "=== ALL DONE ==="
