set -e
{
name="$1"
GRCh38="$sfbdata/genomes/GENCODE/GRCh38/GRCh38_cr4.0"
GRCm38="$sfbdata/genomes/GENCODE/GRCm38/GRCm38_cr4.0"

# run cellranger count
out_dir_name="cr_count-$name"
cellranger count --id="$out_dir_name" \
                --transcriptome="$GRCh38"\
                --fastqs=./Illumina/fastq/"$name" \

mv "$out_dir_name" ./Illumina/cr_count
mv "__$out_dir_name.mro" ./Illumina/cr_count

# run cellranger aggr
out_dir_name="cr_aggr-$name"
cellranger aggr --id="$out_dir_name"  --csv "./Illumina/aggr_csv/aggr-${name}.csv" --localcores=64
mv "$out_dir_name" ./Illumina/cr_aggr
mv "__$out_dir_name.mro" ./Illumina/cr_aggr
exit;
}