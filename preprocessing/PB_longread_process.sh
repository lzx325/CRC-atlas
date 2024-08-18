{
set -e

# set paths

name="$1"
dataset_root="<root_dir>"
cd "$dataset_root"
references_root="<references_root>"
hg38="$references_root/genomes/GENCODE/GRCh38/GRCh38.primary_assembly.genome.fa"
cupcake_root="<cupcake_root>"

# run PacBio ccs
mkdir -p ccs
ccs "subreads/${name}.subreads.bam" ccs/"$name".ccs.bam --min-rq 0.9

# run PacBio lima, tag, and refine
mkdir -p fltnc

lima --isoseq "ccs/${name}.ccs.bam" "primers/primers.fasta" "fltnc/${name}.fl.bam"
isoseq3 tag "fltnc/${name}.fl.10x_TSO_5p--TruSeq_Read1_3p.bam" "fltnc/${name}.flt.bam" --design T-12U-16B
isoseq3 refine "fltnc/${name}.flt.bam" "primers/primers.fasta" "fltnc/${name}.fltnc.bam" --require-polya

# run PacBio dedup
mkdir -p dedup
isoseq3 dedup "fofn/${name}.fofn" "dedup/${name}.dedup.bam" --log-level INFO

# run minimap2
mkdir -p minimap2
minimap2 -t 30 -ax splice -uf --secondary=no -C5 \
"$hg38" "dedup/${name}.dedup.fasta" \
> "minimap2/${name}.aligned.sam" 2> "minimap2/${name}.aligned.sam.log"

# convert sam to bam, sort, and index
{
    samtools view -b "minimap2/${name}.aligned.sam" > "minimap2/${name}.aligned.bam"
    samtools sort "minimap2/${name}.aligned.bam" -o "minimap2/${name}.aligned.sorted.bam"
    samtools index "minimap2/${name}.aligned.sorted.bam"
    sort -k 3,3 -k 4,4n "minimap2/${name}.aligned.sam" > "minimap2/${name}.aligned.sorted.sam"
}

# run Cupcake to collapse isoforms
{
    cd "$dataset_root"
    mkdir -p downstream
    mkdir -p downstream/"$name"
    cd downstream/"$name"
    if ! [ -L "${name}.aligned.sorted.sam" ]; then 
        ln -s ../../minimap2/"${name}.aligned.sorted.sam" "${name}.aligned.sorted.sam"
    fi

    if ! [ -L "${name}.aligned.sorted.bam" ]; then 
        ln -s ../../minimap2/"${name}.aligned.sorted.bam" "${name}.aligned.sorted.bam"
    fi

    if ! [ -L "${name}.dedup.info.csv" ]; then 
        ln -s ../../dedup/"${name}.dedup.info.csv" "${name}.dedup.info.csv"
    fi

    if ! [ -L "${name}.dedup.fasta" ]; then 
        ln -s ../../dedup/"${name}.dedup.fasta" "${name}.dedup.fasta"
    fi

    collapse_isoforms_by_sam.py --input "$name.dedup.fasta" \
        --bam "$name.aligned.sorted.bam" -c 0.99 -i 0.95 \
        --gen_mol_count \
        -o "$name.aligned.5merge" \
        --cpus 20
}

# run sqanti3_qc.py
SQUANTI3_output_dn="SQANTI3_ref_gtf_merged_gencode_v46_refseq202310"
{
   # activate SQANTI3 environment
   cd "$dataset_root"
   cd downstream/"$name"
   conda deactivate
   conda activate SQANTI3.env
   export PATH="$PATH:<path to SQANTI3>"
   export PYTHONPATH="${cupcake_root}/sequence"
   hg38="$references_root/genomes/GENCODE/GRCh38/GRCh38.primary_assembly.genome.fa"
   gtf="$references_root/genomes/GENCODE/GRCh38/gencode.v37.annotation.filtered.gtf"
   gtf_merged="$references_root/genomes/GENCODE_refseq_combined-GENCODE_v46/gffcmp.combined.namecvt.sorted.gtf"
   
   genePred_fp="$SQUANTI3_output_dn/refAnnotation_${name}.aligned.5merge.collapsed.genePred"
   if [ -f "$genePred_fp" ]; then
      rm "$genePred_fp"
   fi
   sqanti3_qc.py \
      --gtf "$name".aligned.5merge.collapsed.gff \
      "$gtf_merged" "$hg38" \
      --fl_count "$name".aligned.5merge.collapsed.abundance.txt \
      --cage_peak "$references_root/TSS/human.refTSS_v3.1.hg38.bed" \
      --polyA_motif_list "$references_root/polyA/human.polyA.list.txt" \
      --polyA_peak "$references_root/polyA/atlas.clusters.2.0.GRCh38.96.chr_corrected.bed" \
      --dir  "$SQUANTI3_output_dn"
}

# run sqanti3_RulesFilter.py
{
   cd "$dataset_root"
   cd downstream/"$name"
   conda deactivate
   conda activate SQANTI3.env
   export PATH="$PATH:<path to SQANTI3>"
   export PYTHONPATH="${cupcake_root}/sequence"
   # The output is "$name".aligned.5merge.collapsed_classification.filtered_lite_classification.txt
   sqanti3_RulesFilter.py \
       "$SQUANTI3_output_dn"/"$name".aligned.5merge.collapsed_classification.txt \
       "$SQUANTI3_output_dn"/"$name".aligned.5merge.collapsed_corrected.fasta \
       "$name".aligned.5merge.collapsed.gff
}

# run collate_FLNC_gene_info.py
{
   conda deactivate
   conda activate bio
   cd "$dataset_root"
   cd downstream/"$name"
   collate_FLNC_gene_info="$cupcake_root/singlecell/collate_FLNC_gene_info.py"
   # "$name".aligned.annotated.csv is the output filename
   gzip -c "$name".dedup.info.csv > "$name".dedup.info.csv.gz
   annotated_csv_fp="$SQUANTI3_output_dn"/"$name".aligned.annotated.csv

   if [ -f "$annotated_csv_fp" ]; then
      rm "$annotated_csv_fp"
   fi
   
   echo python "$collate_FLNC_gene_info" \
     "$name".aligned.5merge.collapsed.group.txt \
     "$name".dedup.info.csv.gz \
     "$SQUANTI3_output_dn"/"$name".aligned.5merge.collapsed_classification.filtered_lite_classification.txt \
     "$annotated_csv_fp"

   python "$collate_FLNC_gene_info" \
     "$name".aligned.5merge.collapsed.group.txt \
     "$name".dedup.info.csv.gz \
     "$SQUANTI3_output_dn"/"$name".aligned.5merge.collapsed_classification.filtered_lite_classification.txt \
     "$annotated_csv_fp"
}
}