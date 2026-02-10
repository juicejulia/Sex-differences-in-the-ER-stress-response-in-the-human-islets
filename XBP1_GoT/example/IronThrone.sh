#!/bin/bash
#SBATCH -n 25
#SBATCH -p bms_q
#SBATCH --mail-type ALL
#SBATCH --mem=100G

for r1 in *_L001_R1_001.fastq; do
    # strip the R1 suffix to get the sample prefix
    prefix=${r1%%_L001_R1_001.fastq}

    # construct the matching R2 filename
    r2="${prefix}_L001_R2_001.fastq"

    # sanity check R2 exists
    if [[ ! -f "$r2" ]]; then
        echo " Missing R2 for $r1 (expected $r2), skipping."
        continue
    fi

    # set sample name = prefix (everything before _L001_)
    sample="$prefix"

    # make per-sample outdir (optional but cleaner)
    outdir="${sample}_GoT"
    mkdir -p "$outdir"

    echo "Running IronThrone-GoT for sample: $sample"
    ./IronThrone-GoT \
	--umilen 12 \
        --run linear \
        --fastqR1 "$r1" \
        --fastqR2 "$r2" \
        --config config.config \
        --sample "$sample" \
        --whitelist barcodes_trim.tsv \
        --outdir "$outdir" \
        --thread 25
done

