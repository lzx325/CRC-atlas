import pandas as pd
import numpy as np
import re
import os
from copy import copy
from os.path import join
from gtfparse import read_gtf
import csv
import pickle as pkl
from pprint import pprint
from BCBio.GFF import GFFExaminer
from utils.gtf_io import my_gtf_read,my_gtf_write,CategorizedGTF
from utils.isoform_classes import genePredReader, genePredRecord, reference_genePred_TSS_TTS_parser
from utils.transcriptomic_features.shared_classes import \
    SharedTTS, SharedTSS, SharedExon, SharedJunction, SharedSpliceSite, SharedTranscriptStructure
from utils.transcriptomic_features.ref_classes import RefTranscriptIsoform
from utils.transcriptomic_features.parse_transcriptomic_features import reference_genePred_parser
from utils.isoform_analysis import ReferenceSharedFeatureCollection,GffCompareMergedTranscripts,SupportInfo
from collections import defaultdict
from copy import copy
from bx.intervals import Interval, IntervalTree
from matplotlib_venn import venn3_unweighted,venn3
import matplotlib as mpl
import matplotlib.pyplot as plt
from tqdm import tqdm
import time
from utils.file_utils import read_bed,write_bed
import pyfaidx
from utils.transcriptomic_features.ref_classes import TranscriptIsoformWithSeq
from Bio.Seq import Seq
from Bio import SeqIO
import seaborn as sns
from scipy.stats import mannwhitneyu
import scipy.sparse as sp
import anndata as ad
from sklearn.preprocessing import normalize

from utils.isoform_analysis import GroupedIsoformV2

def get_sqanti_sample_fp(sample_name):
    SQUANTI3_output_dn="SQANTI3_ref_gtf_merged_2"
    pacbio_analysis_root="PacBio/downstream"
    sqanti_dir=os.path.join(pacbio_analysis_root,"{}/{}".format(sample_name,SQUANTI3_output_dn))
    sqanti_sample_fp_dict = {
        "filtered_lite_classification":os.path.join(
            sqanti_dir,"{}.aligned.5merge.collapsed_classification.filtered_lite_classification.txt".format(sample_name)
        ),
        "filtered_lite_gtf":os.path.join(
            sqanti_dir,"{}.aligned.5merge.collapsed_classification.filtered_lite.gtf".format(sample_name)
        ),
        "filtered_lite_junctions":os.path.join(
            sqanti_dir,"{}.aligned.5merge.collapsed_classification.filtered_lite_junctions.txt".format(sample_name)
        ),
        "collate_df":os.path.join(
            sqanti_dir,"{}.aligned.annotated.csv".format(sample_name)
        ),
        "genePred":os.path.join(
            sqanti_dir,"{}.aligned.5merge.collapsed_corrected.genePred".format(sample_name)
        )
    }
    for fp in sqanti_sample_fp_dict.values():
        assert os.path.exists(fp)
    return sqanti_sample_fp_dict

UMI_threshold={
    'PS017-T2':100,
    'PS017-T1.1':100,
    'PS018-N':100,
    'PS018-T':100,
    'PS021-N':15,
    'PS025-T':4,
    'PS026-N2':100,
    'PS026-T':4,
    'PS028-T':100
}
sample_name_map={
    1:"PS017-T1.1",
    2:"PS017-T2",
    3:"PS018-N",
    4:"PS018-T",
    5:"PS021-N",
    6:"PS025-T",
    7:"PS026-N2",
    8:"PS026-T",
    9:"PS028-T"
}

gffcmp_fp={
    "genePred_fp":"PacBio/downstream/combined/combined_filter_gtfs-colon/gffcmp.query.all_samples.combined.genePred",
    "gffcmp_gtf_fp":"PacBio/downstream/combined/combined_filter_gtfs-colon/gffcmp.query.all_samples.combined.gtf",
    "gffcmp_tracking_fp":"PacBio/downstream/combined/combined_filter_gtfs-colon/gffcmp.query.all_samples.tracking",
    "filtered_lite_classification":"PacBio/downstream/combined/combined_filter_gtfs-colon/SQANTI3_ref_gtf_merged_gencode_v46_refseq202310/gffcmp.query.all_samples.combined_classification.txt"
}

# path to the Illumina anndata with cluster labels
anndata_fp="expr.clusterfull.enterocyte_split_cl_correct.h5ad"
# path to the cellranger aggregation csv
aggr_fp="aggregation.csv"

# path to the sample name mapping
sample_name_mapping_fp="PacBio/sample_table/PacBio_Illumina_mapping-colon.csv"

illumina_sample_info_args=dict(
    anndata_fp=anndata_fp,
    sample_name_mapping_fp=sample_name_mapping_fp,
    aggr_fp=aggr_fp,
    anndata_field_names={
        "cluster_midway":"ClusterMidway2",
        "cluster_full":"ClusterFull2",
        "condition":"condition"
    },
    aggr_field_names={
        "illumina_sample_name":"library_id",
        "condition":"condition"
    }
)

reference_fp={
    "reference_genePred_fp":"PacBio/downstream/combined/refAnnotation.fromgtf.genePred",
    "TSS_bed_fp":"references/TSS/TSS.merged.formatted.win48.bed",
    "TTS_bed_fp":"references/polyA/TTS.merged.formatted.win48.bed"
}

sample_analysis_fp={
    k:get_sqanti_sample_fp(k) for k in sample_name_map.values()
}

giv2=GroupedIsoformV2(
    sample_name_map=sample_name_map,
    sample_analysis_fp=sample_analysis_fp,
    reference_fp=reference_fp,
    UMI_threshold=UMI_threshold,
    load_ref_collection=True,
    gffcmp_fp=gffcmp_fp,
    illumina_smaple_info_args=illumina_sample_info_args,
    test_mode=False
)

giv2.parse_query_transcript_isoform(parse_exon=False,parse_junction=True,parse_tss_tts=False)

giv2.parse_gffcmp_merged_transcript_isoform()

giv2.differential_analysis(
    chained_transcripts=False,
    ref_collections=False,
    gffcmp_transcripts=True,
    suppa2_events=False
)

def make_simplified_isoform_level_table(self):
    def get_supporting_query_isoform_info(
        feature
    ):
        assert 'cls_record' in feature.associated_data
        cls_record=feature.associated_data.get("cls_record")
        SupportingStructCategory=cls_record["structural_category"]
        gene_id=cls_record["associated_gene"]
        NExons=cls_record["exons"]
        if gene_id in self.gene_name_by_id:
            gene_display_name=self.gene_name_by_id[gene_id]
        else:
            gene_display_name=gene_id
            
        AssociatedGeneSymbols=gene_display_name
        
        if "supporting_query_isoform_id" in feature.associated_data \
        and len(feature.associated_data.get("supporting_query_isoform_id"))>0:
            cage_peak_info=list()
            polya_peak_info=list()
            supporting_qtiso=feature.associated_data.get("supporting_query_isoform_id")
            for qtiso in supporting_qtiso:
                cage_peak_info.append(self.query_iso_by_id[qtiso].cls_record["within_cage_peak"])
                polya_peak_info.append(self.query_iso_by_id[qtiso].cls_record["within_polya_site"])
                WithinCAGEPeak=sum(cage_peak_info)>0
                WithinPolyAPeak=sum(polya_peak_info)>0
        else:
            WithinCAGEPeak=False
            WithinPolyAPeak=False
        return {
            "GeneSymbol":AssociatedGeneSymbols,
            "SupportingStructCategory":SupportingStructCategory,
            "WithinCAGEPeak":WithinCAGEPeak,
            "WithinPolyAPeak":WithinPolyAPeak,
            "NExons":NExons
        }
    
    collection=self.gffcmp_merged_transcripts.shared_iso_by_id
    table_dict=defaultdict(list)
    for feature_id,feature in tqdm(collection.items(),total=len(collection)):
        record = dict()
        record["Isoform"]=feature.id
        record["chrom"]=feature.chrom
        record["start"]=feature.genomic_start
        record["end"]=feature.genomic_end_1based
        result=get_supporting_query_isoform_info(feature)
        record.update(result)
        record.update(result)
        for k,v in record.items():
            table_dict[k].append(v)
    return pd.DataFrame(table_dict)

def isoform_level_count_matrix(self):
    obs_names=list(self.cell_barcodes_all)
    obs_names2index=dict(zip(obs_names,range(len(obs_names))))
    var_names=[
        feature_id 
        for feature_id, feature in self.gffcmp_merged_transcripts.shared_iso_by_id.items() 
        if "differential_analysis" in feature.associated_data
    ]
    var_names2index=dict(zip(var_names,range(len(var_names))))
    raw_count_matrix_data=list()
    raw_count_matrix_row=list()
    raw_count_matrix_col=list()
    
    print("constructing isoform level count matrix:")
    
    for feature_id,feature in tqdm(self.gffcmp_merged_transcripts.shared_iso_by_id.items()):
        if "differential_analysis" in feature.associated_data: 
            support_info=feature.associated_data["differential_analysis"]["support_info"]
            CB_UMI_counts=support_info.get_supporting_CB_UMI_counts()
            CB_UMI_counts={k:v for k,v in CB_UMI_counts.items() if k in obs_names2index}
            
            for cb, umi_count in CB_UMI_counts.items():
                value=umi_count
                row_idx=obs_names2index[cb]
                col_idx=var_names2index[feature_id]
                raw_count_matrix_data.append(value)
                raw_count_matrix_row.append(row_idx)
                raw_count_matrix_col.append(col_idx)
            
    counts_matrix_csr=sp.csr_matrix(
        (raw_count_matrix_data,(raw_count_matrix_row,raw_count_matrix_col)),
        shape=(len(obs_names2index),len(var_names2index))
    )
    
    counts_matrix_ad=ad.AnnData(counts_matrix_csr,dtype=np.int64)
    counts_matrix_ad.obs_names=obs_names
    counts_matrix_ad.var_names=var_names
    
    data_matrix_csr=(normalize(counts_matrix_csr,axis=1,norm="l1")*1e4).log1p()
    data_matrix_ad=ad.AnnData(data_matrix_csr,dtype=np.float64)
    data_matrix_ad.obs_names=obs_names
    data_matrix_ad.var_names=var_names
    return counts_matrix_ad, data_matrix_ad

def gene_level_count_matrix(self,simplified_isoform_level_table,isoform_counts_ad):
    obs_names=list(self.cell_barcodes_all)
    obs_names2index=dict(zip(obs_names,range(len(obs_names))))
    var_names=sorted(simplified_isoform_level_table["GeneSymbol"].unique())
    var_names2index=dict(zip(var_names,range(len(var_names))))
    raw_count_matrix_data=list()
    raw_count_matrix_row=list()
    raw_count_matrix_col=list()
    groupby=simplified_isoform_level_table.groupby("GeneSymbol")
    for gene,group in tqdm(groupby):
        support_info=SupportInfo()
        for iso_id in group["Isoform"]:
            support_info.add_isoform_support(
                self.gffcmp_merged_transcripts.shared_iso_by_id[iso_id]\
                .associated_data.get("differential_analysis")["support_info"]
            )
        CB_UMI_counts=support_info.get_supporting_CB_UMI_counts()
        CB_UMI_counts={k:v for k,v in CB_UMI_counts.items() if k in obs_names2index}

        for cb, umi_count in CB_UMI_counts.items():
            value=umi_count
            row_idx=obs_names2index[cb]
            col_idx=var_names2index[gene]
            raw_count_matrix_data.append(value)
            raw_count_matrix_row.append(row_idx)
            raw_count_matrix_col.append(col_idx)
    counts_matrix_csr=sp.csr_matrix(
        (raw_count_matrix_data,(raw_count_matrix_row,raw_count_matrix_col)),
        shape=(len(obs_names2index),len(var_names2index))
    )
    
    counts_matrix_ad=ad.AnnData(counts_matrix_csr,dtype=np.int64)
    counts_matrix_ad.obs_names=obs_names
    counts_matrix_ad.var_names=var_names
    
    denominator=np.repeat(np.asarray(isoform_counts_ad.X.sum(axis=1))[:,0], counts_matrix_csr.getnnz(axis=1)).astype("float64")
    counts_matrix_csr_float=counts_matrix_csr.astype("float64")
    counts_matrix_csr_float.data/=(denominator+np.finfo("float64").eps)
    
    counts_matrix_csr_normalized=(counts_matrix_csr_float*1e4).log1p()
    data_matrix_csr=counts_matrix_csr_normalized
    data_matrix_ad=ad.AnnData(data_matrix_csr,dtype=np.float64)
    data_matrix_ad.obs_names=obs_names
    data_matrix_ad.var_names=var_names
    return counts_matrix_ad, data_matrix_ad

def collect_anndata(self):
    isoform_counts_matrix_ad, isoform_data_matrix_ad = isoform_level_count_matrix(self)
    simplified_isoform_level_table=make_simplified_isoform_level_table(self)
    simplified_isoform_level_table_indexed=simplified_isoform_level_table.set_index("Isoform")
    assert isoform_counts_matrix_ad.shape[1]==simplified_isoform_level_table.shape[0]
    isoform_metadata=simplified_isoform_level_table_indexed.loc[isoform_counts_matrix_ad.var.index,:]
    
    
    def get_cell_metadata():
        cells_in_illumina=self.cell_barcodes_all.intersection(self.illumina_sample_info.anndata.obs.index)
        cells_not_in_illumina=self.cell_barcodes_all.difference(self.illumina_sample_info.anndata.obs.index)
        metadata=self.illumina_sample_info.anndata.obs.loc[list(cells_in_illumina),:]

        cluster_midway_field_name=self.illumina_sample_info.anndata_field_names["cluster_midway"]
        df_append=pd.DataFrame(
            { cluster_midway_field_name:[None]*len(cells_not_in_illumina) },
            index=list(cells_not_in_illumina)
        )
        metadata = pd.concat([metadata,df_append],axis=0)
        metadata[cluster_midway_field_name]=metadata[cluster_midway_field_name].astype(str)
        metadata.loc[list(self.cluster_midway_sets["_PB_CB_not_found_in_IL"]),cluster_midway_field_name]="_PB_CB_not_found_in_IL"
        metadata.loc[list(self.cluster_midway_sets["_PB_CB_filtered"]),cluster_midway_field_name]="_PB_CB_filtered"
        
        metadata=metadata.loc[isoform_counts_matrix_ad.obs.index,:]
        return metadata
    
    def post_process_metadata(metadata):
        metadata=metadata.rename({"library_id":"library_id_IL","library_id_short":"library_id_short_IL"},axis=1)
        metadata["library_id_PB"]=metadata.index.map(
            lambda x: giv2.illumina_sample_info.pacbio_sample_index_revmap[int(x.split('-')[1])]
        )
        return metadata
    
    
    
    cell_metadata=get_cell_metadata()
    cell_metadata=post_process_metadata(cell_metadata)
    
    isoform_counts_matrix_ad.var=isoform_metadata
    isoform_data_matrix_ad.var=isoform_metadata
    isoform_counts_matrix_ad.obs=cell_metadata
    isoform_data_matrix_ad.obs=cell_metadata
    
    gene_counts_matrix_ad, gene_data_matrix_ad = gene_level_count_matrix(self,simplified_isoform_level_table,isoform_counts_matrix_ad)
    gene_counts_matrix_ad.obs=cell_metadata
    gene_data_matrix_ad.obs=cell_metadata
    
    return dict(
        isoform_counts_matrix_ad=isoform_counts_matrix_ad,
        isoform_data_matrix_ad=isoform_data_matrix_ad,
        gene_counts_matrix_ad=gene_counts_matrix_ad,
        gene_data_matrix_ad=gene_data_matrix_ad
    )

anndata_dict=collect_anndata(giv2)

out_dir="anndata_dir"
os.makedirs(out_dir,exist_ok=True)

anndata_dict["isoform_counts_matrix_ad"].write_h5ad(
    os.path.join(
        out_dir,
        "PacBio-isoform_counts_matrix_ad.h5ad"
    )
)

anndata_dict["isoform_data_matrix_ad"].write_h5ad(
    os.path.join(
        out_dir,
        "PacBio-isoform_data_matrix_ad.h5ad"
    )
)

anndata_dict["gene_counts_matrix_ad"].write_h5ad(
    os.path.join(
        out_dir,
        "PacBio-gene_counts_matrix_ad.h5ad"
    )
)

anndata_dict["gene_data_matrix_ad"].write_h5ad(
    os.path.join(
        out_dir,
        "PacBio-gene_data_matrix_ad.h5ad"
    )
)

