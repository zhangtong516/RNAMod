prefix=$1

WINDOW_SIZE=100

##
sed '1d' ../batch2_config.csv | while IFS=',' read -r s p l f; do 
    if [ -f "${s}_positive_sequences.fa" ]; then
        rm -f ${s}_positive_sequences.fa
    fi

    awk -v window=$WINDOW_SIZE 'BEGIN{OFS="\t"} {
        summit = $2 + $10;
        start = summit - window;
        end = summit + window;
        if(start < 0) start = 0;
        print $1, start, end, $4"_plus", $5, "+";
        print $1, start, end, $4"_minus", $5, "-";
    }' $p > ${s}_positive_regions.bed
    bedtools getfasta -fi ${f} -bed ${s}_positive_regions.bed -fo ${s}_positive_sequences.fa

done

if [ -f "${prefix}_sequences_combined.fa" ]; then
    rm -f ${prefix}_sequences_combined.fa
fi

cat *_positive_sequences.fa > ${prefix}_sequences_combined.fa

## in case the sequences are too many:
seqkit sample -p 0.25 ${prefix}_sequences_combined.fa > ${prefix}_sequences_combined_subsample.fa
#ln -s  ${prefix}_sequences_combined.fa ${prefix}_sequences_combined_subsample.fa
## run streme
module load singularity
memeImage="/home/users/astar/gis/zhangt/scratch/tools/meme_v5.5.8.sif"

singularity run $memeImage streme --p ${prefix}_sequences_combined_subsample.fa \
    --rna \
    --minw 5 \
    --maxw 6 \
    --nmotifs 15 \
    --thresh 0.01 \
    --niter 100 \
    --time 600 \
    --order 2 \
    --nsubsets 5 \
    --seed 12345 \
    --oc streme_${prefix} 

##