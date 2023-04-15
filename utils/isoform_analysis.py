from concurrent.futures import process
import os
import sys
from turtle import update
from matplotlib import cbook
from scipy.stats import fisher_exact
import numpy as np
from os.path import join,basename,dirname,splitext
from glob import glob
import numpy as np
import scipy.sparse as sp
import pandas as pd
from collections import defaultdict,OrderedDict,deque,Counter
from itertools import chain
from gtfparse import read_gtf
import re
from typing import Dict, DefaultDict, List, Set, Tuple
from torch import negative
from tqdm import tqdm
import anndata
from copy import copy

from utils.transcriptomic_features.parse_transcriptomic_features import reference_genePred_parser, reference_TSS_TTS_parser, \
init_transcriptomic_feature_collections, update_shared_entity_type_id

from utils.transcriptomic_features.ref_classes import RefTranscriptIsoform
from utils.transcriptomic_features.query_classes import QueryTranscriptIsoform
from utils.isoform_classes import genePredReader, genePredRecord

from utils.transcriptomic_features.shared_classes import \
    SharedExon,SharedJunction,SharedSpliceSite,SharedTSS,SharedTTS,SUPPA2Event,SharedTranscriptIsoform
from utils.file_utils import iprint
from utils.gtf_io import CategorizedGTF,my_gtf_read,my_gtf_write
from utils.utils import interval_overlap, canonical_chr_regex
from utils.file_utils import line_count
TEST_CHROM="chr21"
SEP_SIZE=15
class GroupedIsoForm(object):
    def __init__(
        self,
        samples,
        UMI_threshold=dict(),
        pb_gtf=True,
        dataset_files_dir="./PacBio",
        squanti_dir_name="SQANTI3"
    ):
        self.dataset_files_dir=dataset_files_dir
        self.downstream_analysis_dir=join(dataset_files_dir,"downstream")
        self.sqanti_dir_name=squanti_dir_name
        self.sample_names: List[str]=samples
        self.UMI_threshold: Dict =UMI_threshold
        self.__setup_class_df()
        self.__setup_collate_df()
        self.__collate_df_filter_cells()
        self.__get_structural_category_coarse()
        self.class_df_all=self.__class_df_add_supporting_cb_umi_count(self.class_df_all,"isoform","pbid")
        if pb_gtf:
            self.__add_gtf_info()
    
    def __collate_df_filter_cells(self):
        if not self.UMI_threshold:
            return
        assert self.collate_df_all["id"].is_unique
        collate_df_filtered_list=list()
        print("="*SEP_SIZE+"filtering collate df"+"="*SEP_SIZE)

        for sn in self.sample_names:

            collate_df=self.collate_df_all.query(f"sample_name=='{sn}'")
            if sn in self.UMI_threshold:
                thresh=self.UMI_threshold[sn]
                cb_umi_count=collate_df.groupby("BC")["UMI"].nunique().sort_values(ascending=False)
                selected_barcodes=cb_umi_count.index[cb_umi_count>=thresh]
                collate_df_filtered=collate_df[collate_df["BC"].isin(selected_barcodes)]
                collate_df_filtered_list.append(collate_df_filtered)
                print("sample:",sn,"original #CB:",len(cb_umi_count),"after filter #CB:",len(selected_barcodes))
            else:
                collate_df_filtered_list.append(collate_df)
        self.collate_df_all=pd.concat(collate_df_filtered_list,axis=0,ignore_index=True)

    def __get_structural_category_coarse(self):
        def mapper(x):
            return{
                'incomplete-splice_match':"ISM",
                'full-splice_match':"FSM",
                'novel_in_catalog':"NIC",
                'novel_not_in_catalog':"NNC"
            }.get(x,"other")
        self.class_df_all["structural_category_coarse"]=self.class_df_all["structural_category"].map(mapper)

    def __setup_class_df(self):
        df_list=list()
        for name in self.sample_names:
            df=pd.read_csv(join(self.downstream_analysis_dir,name,self.sqanti_dir_name,
                                f"{name}.aligned.5merge.collapsed_classification.filtered_lite_classification.txt"),sep="\t")
            df["sample_name"]=name
            df["isoform"]=df["isoform"].map(lambda x: f"{name}.{x}")
            df_list.append(df)
        df_all=pd.concat(df_list,axis=0,ignore_index=True)
        df_all=df_all.rename({'FL':'FL_count'},axis=1)
        self.class_df_all=df_all
    
    def __class_df_add_supporting_cb_umi_count(self,class_df,class_df_key,collate_df_key):
        collate_df_cb_umi_count=self.collate_df_all.groupby(collate_df_key).agg({"BCrev":"nunique","UMI":"nunique"})
        class_df=pd.merge(class_df,collate_df_cb_umi_count,left_on=class_df_key,right_index=True,how="left")
        class_df=class_df.rename({'BCrev':'n_cells','UMI':'n_umi'},axis=1)
        class_df["n_cells"]=class_df["n_cells"].fillna(0)
        class_df["n_umi"]=class_df["n_umi"].fillna(0)
        return class_df
    def __setup_collate_df(self):
        df_list=list()
        for name in self.sample_names:
            df=pd.read_csv(join(self.downstream_analysis_dir,name,self.sqanti_dir_name,
                                f"{name}.aligned.annotated.csv"),sep="\t")
            df["sample_name"]=name
            df["id"]=df["id"].map(lambda x: f"{name}/{x}")
            df["pbid"]=df["pbid"].map(lambda x: f"{name}.{x}")
            df_list.append(df)
        self.collate_df_all=pd.concat(df_list,axis=0,ignore_index=True)

    def __add_gtf_info(self):
        def gtf_group_mapper(group):
            group_transcript=group[group.feature=="transcript"]
            assert len(group_transcript)==1
            group_exon=group[group.feature=="exon"]
            group_exon.sort_values("start")
            assert group_exon["start"].min()==group_transcript["start"].iloc[0]
            assert group_exon["end"].max()==group_transcript["end"].iloc[0]
            group_transcript["exon_starts"]=[tuple(group_exon["start"])]
            group_transcript["exon_ends"]=[tuple(group_exon["end"])]
#             if len(group_exon["start"])>=2:
#                 group_transcript["exon_starts_inner"]=[tuple(group_exon["start"][1:])]
#                 group_transcript["exon_ends_inner"]=[tuple(group_exon["end"][:-1])]
#             else:
#                 group_transcript["exon_starts_inner"]=group_transcript["exon_starts"]
#                 group_transcript["exon_ends_inner"]=group_transcript["exon_ends"]
            return group_transcript
        def collapse_transcript_info(df):
            transcript_dict=OrderedDict(iter(df.query("feature=='transcript'").groupby("transcript_id")))
            exon_dict=OrderedDict(iter(df.query("feature=='exon'").groupby("transcript_id")))
            
            def loop_fn(transcript_id):
                group_transcript=transcript_dict[transcript_id]
                group_exon=exon_dict[transcript_id]
                assert len(group_transcript)==1
                group_exon=group_exon.sort_values("start")
                assert group_exon["start"].min()==group_transcript["start"].iloc[0]
                assert group_exon["end"].max()==group_transcript["end"].iloc[0]
                group_transcript["exon_starts"]=[tuple(group_exon["start"])]
                group_transcript["exon_ends"]=[tuple(group_exon["end"])]
                return group_transcript
            
#             with joblib.Parallel(n_jobs=30) as parallel: 
#                 res=parallel(joblib.delayed(loop_fn) (tid) for tid in transcript_dict.keys())
            res=(loop_fn(tid) for tid in transcript_dict.keys())
            return pd.concat(res,axis=0,ignore_index=True)
        df_list=list()
        for name in self.sample_names:
            print("processing:",name)
            df=read_gtf(join(self.downstream_analysis_dir,name,self.sqanti_dir_name,
                            f"{name}.aligned.5merge.collapsed_classification.filtered_lite.gtf"))
            print("finish reading gtf:",name)
            df_grouped_result=collapse_transcript_info(df)
            print("finish mapping:",name)
            df_grouped_result["sample_name"]=name
            df_grouped_result["transcript_id"]=df_grouped_result["transcript_id"].map(lambda x: f"{name}.{x}")
            df_grouped_result["gene_id"]=df_grouped_result["gene_id"].map(lambda x: f"{name}.{x}")
            df_list.append(df_grouped_result)
        df_all=pd.concat(df_list,axis=0,ignore_index=True)
        self.gtf_all=df_all

    def dedup_across_samples(self):
        loc2newtid=dict()
        old2newtid=dict()

        if hasattr(self,"class_df_all"):
            print("dedup using class_df_all")
            df_all=self.class_df_all.set_index("isoform")
        else:
            assert False
        
        def get_key(meta,row,isoform_id):
            scc=meta['structural_category_coarse']
            if scc in ["FSM"]:
                return ("FSM",meta["associated_transcript"])
            elif scc in ["ISM"]:
                return ("ISM",meta["associated_transcript"])
            elif scc in ["NIC"]:
                return ("NIC",row["seqname"],row["exon_starts"][1:],row["exon_ends"][0:-1])
            elif scc in ["NNC"]:
                return ("NNC",row["seqname"],row["exon_starts"][1:],row["exon_ends"][0:-1])
            else:
                return ("other",row["seqname"],row["exon_starts"][1:],row["exon_ends"][0:-1])

        for name in self.sample_names:
            gtf_1sample=self.gtf_all[self.gtf_all.sample_name==name]
            for _,row in gtf_1sample.iterrows():
                meta=df_all.loc[row["transcript_id"]]
                key=get_key(meta,row,row["transcript_id"])
                
                if key in loc2newtid:
                    old2newtid[row["transcript_id"]]=loc2newtid[key]
                else:
                    if key is not None:
                        print(key)
                        tid=row["transcript_id"]
                        loc2newtid[key]=f"agg.{tid}"
                        old2newtid[row["transcript_id"]]=loc2newtid[key]
        self.loc2newtid=loc2newtid        
        self.old2newtid=old2newtid

    def collect_dedup_iso_df(self,ENS2SYM=None):
        df_all=self.class_df_all[[
            "isoform","chrom","strand","length","exons","structural_category",'structural_category_coarse',
            "associated_gene","associated_transcript","ref_length","ref_exons","sample_name","FL_count"]].copy()
        self.collate_df_all["newtid"]=self.collate_df_all["pbid"].map(self.old2newtid)
        df_all["newtid"]=df_all["isoform"].map(self.old2newtid)
        self.class_df_all_dedup=df_all.groupby("newtid").agg({
            'chrom':'first',
            'strand':'first',
            'length':'first',
            'exons':'first',
            'structural_category':'first',
            'structural_category_coarse':'first',
            'associated_gene':'first',
            'associated_transcript':'first',
            'ref_length':'first',
            'ref_exons':'first',
            'sample_name':'unique',
            'FL_count':"sum"
        }).reset_index()
        self.class_df_all_dedup=self.__class_df_add_supporting_cb_umi_count(self.class_df_all_dedup,"newtid","newtid")

        if ENS2SYM is not None:
            self.class_df_all_dedup["associated_gene_symbol"]=self.class_df_all_dedup["associated_gene"].str.split('.').map(lambda x: x[0]).map(lambda x: ENS2SYM.get(x,x))
        self.class_df_all_dedup["n_samples"]=self.class_df_all_dedup["sample_name"].map(lambda x:len(x))
        pat=re.compile(r'\.PB\..*')
        pat2=re.compile(r'agg\.')
        self.class_df_all_dedup["first_sample"]=self.class_df_all_dedup["newtid"].map(lambda x: pat2.sub('',pat.sub('',x)))
        self.class_df_all_dedup["sample_name_repr"]=self.class_df_all_dedup["sample_name"].map(lambda x: repr(x.tolist()))
        self.class_df_all_dedup_cell_gt_0=self.class_df_all_dedup.query("n_cells>0").copy()



def get_CDS_intervals(exon_intervals,CDS_genomic_start,CDS_genomic_end_1based,strand):
    assert strand in ("+","-")
    # assert monotonic
    exon_intervals=np.array(exon_intervals)
    exon_intervals_flatten=exon_intervals.flatten()
    assert np.all(exon_intervals_flatten[1:]>exon_intervals_flatten[0:-1])
    cds_list=list()
    if strand=="-":
        exon_intervals=exon_intervals[::-1,:]
    overlap_cumulative=0
    for i in range(exon_intervals.shape[0]):
        (cds_start,cds_end),overlap=interval_overlap((exon_intervals[i,0],exon_intervals[i,1]),(CDS_genomic_start,CDS_genomic_end_1based))
        if overlap>0:
            frame=(3-overlap_cumulative%3)%3
            overlap_cumulative+=overlap
            cds_list.append(((cds_start,cds_end),frame))
    return cds_list
class GroupedIsoformV2(object):
    def __init__(
        self,
        sample_name_map,
        sample_analysis_fp,
        reference_fp,
        sample_names=None,
        load_ref_collection=True,
        UMI_threshold=None,
        suppa2_event_fp=None,
        gffcmp_fp=None,
        chained_fp=None,
        test_mode=False
    ):
        self.test_mode=test_mode
        self.sample_analysis_fp=sample_analysis_fp
        self.reference_fp=reference_fp
        if sample_names is None:
            self.sample_names=list(sample_name_map.values())
        else:
            self.sample_names=sample_names
        self.sample_name_map=sample_name_map
        self.UMI_threshold=UMI_threshold
        self.illumina_sample_info=IlluminaSampleInfo()
        
        self.has_ref_collection=False
        if load_ref_collection:
            self.__load_ref_collection()

        self.load_suppa2_collection=False
        self.suppa2_event_fp=suppa2_event_fp
        self.gffcmp_fp=gffcmp_fp
        self.chained_fp=chained_fp
        self.query_iso_by_id: Dict[str,QueryTranscriptIsoform]=dict()
        self.has_class_df=False
        self.has_collate_info=False
        self.has_gene_id_mapping=False
        self.has_parsed_query_isoforms=False
        self.has_suppa2_event_collection=False
        self.has_differential_analysis=False
        self.has_gffcmp_merged_transcripts=False
        self.has_chained_transcripts=False
        self.__setup_class_df()
        self.__setup_collate_info()
        self.__filter_cells_by_collate_df_UMI()
        self.__filter_cells_by_illumina_barcodes()
        self.__setup_gene_id_mapping()

    def __load_ref_collection(
        self
    ):
        self.ref_collections=ReferenceSharedFeatureCollection(
            self.reference_fp["reference_genePred_fp"],
            self.reference_fp["TSS_bed_fp"],
            self.reference_fp["TTS_bed_fp"],
            self.test_mode
        )
        self.has_ref_collection=True

    def parse_suppa2_events(self):
        assert self.has_class_df and self.has_gffcmp_merged_transcripts
        self.suppa2_event_collection=SUPPA2EventCollection(
            self.suppa2_event_fp,
            self.sample_name_map,
            self.class_df_all,
            self.gffcmp_merged_transcripts.merged_QTISO_mapping,
            test_mode=self.test_mode
        )
        self.has_suppa2_event_collection=True

    def parse_gffcmp_merged_transcript_isoform(self):
        assert self.has_parsed_query_isoforms
        self.gffcmp_merged_transcripts=GffCompareMergedTranscripts(
            gffcmp_fp=self.gffcmp_fp,
            query_iso_by_id=self.query_iso_by_id,
            sample_name_map=self.sample_name_map,
            test_mode=self.test_mode
        )
        self.has_gffcmp_merged_transcripts=True

    def parse_chained_transcript_isoform(self):
        assert self.has_parsed_query_isoforms
        self.chained_transcripts=ChainedTranscripts(
            chained_fp=self.chained_fp,
            query_iso_by_id=self.query_iso_by_id,
            test_mode=self.test_mode
        )
        self.has_chained_transcripts=True

    def __setup_gene_id_mapping(self):
        self.gene_name_by_id=pd.read_csv("/data/liz0f/genomes/GENCODE_refseq_combined/gene_id_to_name.csv",header=None,index_col=0)[1]

    def __setup_class_df(self):
        def process_CDS_start_end(df):
            df=df.copy()
            swapped=df[["strand","CDS_genomic_start","CDS_genomic_end"]].apply(
                lambda x: (x["CDS_genomic_start"],x["CDS_genomic_end"]) 
                if x["strand"]=='+' else (x["CDS_genomic_end"],x["CDS_genomic_start"]), axis=1,result_type="expand"
            )
            df["CDS_genomic_start"]=swapped[0]-1 # make 0-based
            df["CDS_genomic_end"]=swapped[1]
            df["CDS_start"]-=1 # make 0-based
            return df
        def modify_novel_gene_name(df,sample_name):
            df=df.copy()
            df["associated_gene"]=df["associated_gene"].transform(lambda x: f"{sample_name}__{x}" if x.startswith("NovelGene_") else x)
            return df
        iprint("="*SEP_SIZE+"Loading class_df for each sample"+"="*SEP_SIZE)
        df_list=list()
        for name in self.sample_names:
            print(f"Loading class_df for {name}...",end="")
            df=pd.read_csv(self.sample_analysis_fp[name]["filtered_lite_classification"],sep="\t")
            if self.test_mode:
                df=df[df["chrom"]==TEST_CHROM]
            df=process_CDS_start_end(df)
            df=modify_novel_gene_name(df,name)
            df=df.rename({"isoform":"pbid","FL":"FL_count","CDS_genomic_end":"CDS_genomic_end_1based","CDS_end":"CDS_end_1based"},axis=1)
            df["sample_name"]=name
            df["pbid"]=df["pbid"].map(lambda x: f"QTISO__{name}__{x}")
            df["associated_transcript"]=df["associated_transcript"].apply(lambda x: "RTISO__"+x if x!="novel" else "novel")
            
            df_list.append(df)
            print("finished",end="\n")
        df_all=pd.concat(df_list,axis=0,ignore_index=True)
        self.class_df_all=df_all
        self.has_class_df=True

    def __setup_collate_info(self):
        df_list=list()
        isoform_molecule_barcode_map=defaultdict(lambda: {"molecule_id":set(), "BCrev":set(), "BCrev_molecule_id":set()})
        for name in self.sample_names:
            print(f"Loading collate_info for {name}")
            df=pd.read_csv(self.sample_analysis_fp[name]["collate_df"],sep="\t")
            df=df.rename({"id":"molecule_id"},axis=1)
            df["sample_name"]=name
            df["molecule_id"]=df["molecule_id"].map(lambda x: f"{name}/{x}")
            df["pbid"]=df["pbid"].map(lambda x: f"QTISO__{name}__{x}")
            df["BCrev"]=df.apply(lambda row: "{}-{}".format(row["BCrev"],self.illumina_sample_info.pacbio_sample_index_map[row["sample_name"]]),axis=1)
            df_list.append(df)

            for _,row in tqdm(df.iterrows(),total=df.shape[0]):
                isoform_molecule_barcode_map[row["pbid"]]["molecule_id"].add(row["molecule_id"])
                isoform_molecule_barcode_map[row["pbid"]]["BCrev"].add(row["BCrev"])
                isoform_molecule_barcode_map[row["pbid"]]["BCrev_molecule_id"].add((row["BCrev"],row["molecule_id"]))
        self.isoform_molecule_barcode_map=dict(isoform_molecule_barcode_map)
        self.cell_barcodes_all=set()
        for record in self.isoform_molecule_barcode_map.values():
            self.cell_barcodes_all.update(record["BCrev"])
        self.collate_df_all=pd.concat(df_list,axis=0,ignore_index=True)
        self.has_collate_info=True

    def __filter_cells_by_collate_df_UMI(self):
        if self.UMI_threshold:
            assert self.collate_df_all["molecule_id"].is_unique
            cell_barcodes_pass_filter=list()
            print("="*SEP_SIZE+"filter cells by collate df UMI"+"="*SEP_SIZE)

            for sn in self.sample_names:
                collate_df=self.collate_df_all.query(f"sample_name=='{sn}'")
                if sn in self.UMI_threshold:
                    thresh=self.UMI_threshold[sn]
                    cb_umi_count=collate_df.groupby("BCrev")["UMI"].nunique().sort_values(ascending=False)
                    selected_barcodes=cb_umi_count.index[cb_umi_count>=thresh].values
                    cell_barcodes_pass_filter.append(selected_barcodes)
                    print("sample:",sn,"original #CB:",len(cb_umi_count),"after filter #CB:",len(selected_barcodes))
                else:
                    cell_barcodes_pass_filter.append(collate_df["BCrev"].unique().values)

            cell_barcodes_pass_filter=np.concatenate(cell_barcodes_pass_filter)
            self.cell_barcodes_pass_UMI_count_filter=set(cell_barcodes_pass_filter)
            assert len(self.cell_barcodes_pass_UMI_count_filter)==len(cell_barcodes_pass_filter) # assert no duplication

    def __filter_cells_by_illumina_barcodes(self):
        filter_stat=self.illumina_sample_info.get_pacbio_intersection(self.cell_barcodes_pass_UMI_count_filter)

        # add UMI_filtered and PB not found types
        self.cluster_midway_sets=copy(self.illumina_sample_info.cluster_midway_sets_pacbio_intersect)
        self.cell_barcodes_illumina_intersect=copy(self.illumina_sample_info.cell_barcodes_pacbio_intersect)
        self.cluster_midway_sets["_PB_CB_not_found_in_IL"]=self.cell_barcodes_pass_UMI_count_filter.difference(self.cell_barcodes_illumina_intersect)
        self.cluster_midway_sets["_PB_CB_filtered"]=self.cell_barcodes_all.difference(self.cell_barcodes_pass_UMI_count_filter)

        self.sample_condition_sets={
            "sample_normal":set(),"sample_tumor":set()
        }
        for bc in self.cell_barcodes_all:
            sample_index=int(bc.split('-')[1])
            condition=self.illumina_sample_info.illumina_sample_index2condition[sample_index]
            condition="sample_{}".format(condition)
            self.sample_condition_sets[condition].add(bc)

        # sample condition sets
        print("="*SEP_SIZE+"filter cells by Illumina barcodes"+"="*SEP_SIZE)
        # print filter information
        for k,v in filter_stat.items():
            sn=self.illumina_sample_info.pacbio_sample_index_revmap[k]
            print("sample_name: {}, before filter #CB: {}, after filter #CB: {}".format(sn,v[0],v[1]))
        

    def _update_ref_exons_collections(
        self,
        shared_exon_by_id,shared_exon_to_id,transcriptomic_feature_sets_by_gene,
        query_ti,query_cls_rec,
        type_id_map
    ):
        if query_ti.num_exons>=3:
            for exon_i in range(1,query_ti.num_exons-1):
                exon=query_ti.exons[exon_i]
                shared_exon=SharedExon.from_exon(exon)
                if shared_exon in shared_exon_to_id:
                    shared_exon=shared_exon_by_id[shared_exon_to_id[shared_exon]]
                else: 
                    if query_cls_rec["structural_category"] in ("novel_in_catalog","novel_not_in_catalog"):
                        shared_exon=update_shared_entity_type_id(shared_exon,type_id_map)
                        shared_exon_by_id[shared_exon.id]=shared_exon
                        shared_exon_to_id[shared_exon]=shared_exon.id
                        
                # update part
                ## updating "supporting_query_isoform_id"
                shared_exon.update_associated_data("supporting_query_isoform_id",query_ti.id,mode="append")
                if shared_exon.id.startswith("DEXON__"):
                    ## updating "associated_query_id" for novel ones
                    shared_exon.update_associated_data("associated_query_id",exon.id,mode="append")
                    ## updating "associated_gene_id" and feature_sets_by_gene for novel ones
                    shared_exon.update_associated_data("associated_gene_id",(query_ti.id,query_cls_rec["associated_gene"]),mode="append")
                    if query_cls_rec["associated_gene"] in transcriptomic_feature_sets_by_gene["exons"]:
                        transcriptomic_feature_sets_by_gene["exons"][query_cls_rec["associated_gene"]].add(shared_exon)
                    else:
                        transcriptomic_feature_sets_by_gene["exons"][query_cls_rec["associated_gene"]]=set()
                        transcriptomic_feature_sets_by_gene["exons"][query_cls_rec["associated_gene"]].add(shared_exon)

    def _update_ref_junctions_collections(
        self,
        shared_junction_by_id,shared_junction_to_id,shared_donors_by_id, shared_donors_to_id,shared_acceptors_by_id, shared_acceptors_to_id, transcriptomic_feature_sets_by_gene,
        query_ti,query_cls_rec,
        type_id_map
    ):
        def annotate_junction_coding_type(junction: SharedJunction, query_cls_rec: pd.Series):
            if query_cls_rec["coding"]=="coding":
                if junction.strand=="+":
                    if junction.genomic_start<junction.genomic_end_1based<=query_cls_rec.CDS_genomic_start:
                        coding_type=("5_UTR_splice_junction")
                    elif query_cls_rec.CDS_genomic_start <= junction.genomic_start  < junction.genomic_end_1based <= query_cls_rec.CDS_genomic_end_1based:
                        coding_type=("coding_region_splice_junction")
                    elif query_cls_rec.CDS_genomic_end_1based <= junction.genomic_start < junction.genomic_end_1based:
                        coding_type=("3_UTR_splice_junction")
                    else:
                        assert False
                elif junction.strand=="-":
                    if junction.genomic_start<junction.genomic_end_1based<=query_cls_rec.CDS_genomic_start:
                        coding_type=("3_UTR_splice_junction")
                    elif query_cls_rec.CDS_genomic_start <= junction.genomic_start  < junction.genomic_end_1based <= query_cls_rec.CDS_genomic_end_1based:
                        coding_type=("coding_region_splice_junction")
                    elif query_cls_rec.CDS_genomic_end_1based <= junction.genomic_start < junction.genomic_end_1based:
                        coding_type=("5_UTR_splice_junction")
                    else:
                        assert False

            elif query_cls_rec["coding"]=="non_coding":
                coding_type=("non_coding_transcript_splice_junction")
            else:
                assert False
            return coding_type

        for junction_i in range(0,len(query_ti.junctions)):
            junction=query_ti.junctions[junction_i]
            donor=junction.splice_donor
            acceptor=junction.splice_acceptor
            shared_junction=SharedJunction.from_junction(junction)
            shared_donor=SharedSpliceSite.from_splicesite(donor)
            shared_acceptor=SharedSpliceSite.from_splicesite(acceptor)

            coding_type=annotate_junction_coding_type(shared_junction,query_cls_rec)

            if shared_junction in shared_junction_to_id:
                shared_junction=shared_junction_by_id[shared_junction_to_id[shared_junction]]

            else: # novel junction
                if query_cls_rec["structural_category"] in ("novel_in_catalog","novel_not_in_catalog"):
                    shared_junction=update_shared_entity_type_id(shared_junction,type_id_map)
                    shared_junction_by_id[shared_junction.id]=shared_junction
                    shared_junction_to_id[shared_junction]=shared_junction.id

            # update junctions
            shared_junction.update_associated_data("supporting_query_isoform_id",query_ti.id,mode="append")
            shared_junction.update_associated_data("supporting_query_isoform_coding_type",(query_ti.id,coding_type),mode="append")
            if shared_junction.id.startswith("DJUNC__"):
                ## updating "associated_query_id" for novel ones
                shared_junction.update_associated_data("associated_query_id",junction.id,mode="append")
                ## updating "associated_gene_id" and feature_sets_by_gene for novel ones
                shared_junction.update_associated_data("associated_gene_id",(query_ti.id,query_cls_rec["associated_gene"]),mode="append")
                if query_cls_rec["associated_gene"] in transcriptomic_feature_sets_by_gene["junctions"]:
                    transcriptomic_feature_sets_by_gene["junctions"][query_cls_rec["associated_gene"]].add(shared_junction)
                else:
                    transcriptomic_feature_sets_by_gene["junctions"][query_cls_rec["associated_gene"]]=set()
                    transcriptomic_feature_sets_by_gene["junctions"][query_cls_rec["associated_gene"]].add(shared_junction)

            if shared_donor in shared_donors_to_id:
                shared_donor=shared_donors_by_id[shared_donors_to_id[shared_donor]]
                
            else:
                if query_cls_rec["structural_category"] in ("novel_not_in_catalog",):
                    shared_donor=update_shared_entity_type_id(shared_donor,type_id_map)
                    shared_donors_by_id[shared_donor.id]=shared_donor
                    shared_donors_to_id[shared_donor]=shared_donor.id

            # update donors
            shared_donor.update_associated_data("supporting_query_isoform_id",query_ti.id,mode="append")
            shared_donor.update_associated_data("supporting_query_isoform_coding_type",(query_ti.id,coding_type),mode="append")
            if shared_donor.id.startswith("DSD__"):
                shared_donor.update_associated_data("associated_query_id",donor.id,mode="append")
                shared_donor.update_associated_data("associated_gene_id",(query_ti.id,query_cls_rec["associated_gene"]),mode="append")
                if query_cls_rec["associated_gene"] in transcriptomic_feature_sets_by_gene["donors"]:
                    transcriptomic_feature_sets_by_gene["donors"][query_cls_rec["associated_gene"]].add(shared_donor)
                else:
                    transcriptomic_feature_sets_by_gene["donors"][query_cls_rec["associated_gene"]]=set()
                    transcriptomic_feature_sets_by_gene["donors"][query_cls_rec["associated_gene"]].add(shared_donor)

            if shared_acceptor in shared_acceptors_to_id:
                shared_acceptor=shared_acceptors_by_id[shared_acceptors_to_id[shared_acceptor]]
            else:
                if query_cls_rec["structural_category"] in ("novel_not_in_catalog",):
                    shared_acceptor=update_shared_entity_type_id(shared_acceptor,type_id_map)
                    shared_acceptors_by_id[shared_acceptor.id]=shared_acceptor
                    shared_acceptors_to_id[shared_acceptor]=shared_acceptor.id

            # update acceptors
            shared_acceptor.update_associated_data("supporting_query_isoform_id",query_ti.id,mode="append")
            shared_acceptor.update_associated_data("supporting_query_isoform_coding_type",(query_ti.id,coding_type),mode="append")

            if shared_acceptor.id.startswith("DSA__"):
                shared_acceptor.update_associated_data("associated_query_id",acceptor.id,mode="append")
                shared_acceptor.update_associated_data("associated_gene_id",(query_ti.id,query_cls_rec["associated_gene"]),mode="append")
                if query_cls_rec["associated_gene"] in transcriptomic_feature_sets_by_gene["acceptors"]:
                    transcriptomic_feature_sets_by_gene["acceptors"][query_cls_rec["associated_gene"]].add(shared_acceptor)
                else:
                    transcriptomic_feature_sets_by_gene["acceptors"][query_cls_rec["associated_gene"]]=set()
                    transcriptomic_feature_sets_by_gene["acceptors"][query_cls_rec["associated_gene"]].add(shared_acceptor)

    def _update_ref_tss_tts_collections(
        self,
        shared_tss_by_id,shared_tss_to_id,shared_tts_by_id, shared_tts_to_id, tss_tts_itree_by_chr, tss_tts_by_gene,
        query_ti,query_cls_rec,
        type_id_map,
        accept_new=True
    ):
        tss=query_ti.TSS
        tts=query_ti.TTS
        # update TSS record
        if query_cls_rec.structural_category=="full-splice_match" or \
            (query_cls_rec.structural_category=="incomplete-splice_match" and query_cls_rec.subcategory not in ("3prime_fragment","internal_fragment")) or \
            (query_cls_rec.structural_category=="novel_not_in_catalog") or \
            (query_cls_rec.structural_category=="novel_in_catalog"):
            shared_tss_to_update=tss_tts_itree_by_chr["tss"][(tss.chrom,tss.strand)].find(tss.genomic_center,tss.genomic_center+1)
            if len(shared_tss_to_update)>0:
                for shared_tss in shared_tss_to_update:
                    shared_tss.update_associated_data("supporting_query_isoform_id",query_ti.id,mode="append")
            else:
                assert accept_new
                shared_tss=SharedTSS.from_TSS(tss)
                shared_tss.genomic_window_start=shared_tss.genomic_center+1-24
                shared_tss.genomic_window_end_1based=shared_tss.genomic_center+1+24

                shared_tss=update_shared_entity_type_id(shared_tss,type_id_map)
                shared_tss_by_id[shared_tss.id]=shared_tss
                shared_tss_to_id[shared_tss]=shared_tss.id
                tss_tts_itree_by_chr["tss"][(shared_tss.chrom,shared_tss.strand)].insert(
                    shared_tss.genomic_window_start,
                    shared_tss.genomic_window_end_1based,
                    shared_tss
                )
                shared_tss_to_update=[shared_tss]

            # update part
            for shared_tss in shared_tss_to_update:
                ## updating "supporting_query_isoform_id"
                shared_tss.update_associated_data("supporting_query_isoform_id",query_ti.id,mode="append")

                ## updating "associated_gene_id" and feature_sets_by_gene
                shared_tss.update_associated_data("associated_gene_id",(query_ti.id,query_cls_rec["associated_gene"]),mode="append")
                if query_cls_rec["associated_gene"] in tss_tts_by_gene["tss"]:
                    tss_tts_by_gene["tss"][query_cls_rec["associated_gene"]].add(shared_tss)
                else:
                    tss_tts_by_gene["tss"][query_cls_rec["associated_gene"]]=set()
                    tss_tts_by_gene["tss"][query_cls_rec["associated_gene"]].add(shared_tss)

                if shared_tss.id.startswith("DTSS__"):
                    ## updating "associated_query_id" for novel ones
                    shared_tss.update_associated_data("associated_query_id",tss.id,mode="append")

                    

        if query_cls_rec.structural_category=="full-splice_match" or \
            (query_cls_rec.structural_category=="incomplete-splice_match" and query_cls_rec.subcategory not in ("5prime_fragment","internal_fragment")) or \
            (query_cls_rec.structural_category=="novel_not_in_catalog") or \
            (query_cls_rec.structural_category=="novel_in_catalog"):

            shared_tts_to_update=tss_tts_itree_by_chr["tts"][(tts.chrom,tts.strand)].find(tts.genomic_center,tts.genomic_center+1)
            if len(shared_tts_to_update)>0:
                for ref_tts in shared_tts_to_update:
                    ref_tts.update_associated_data("supporting_query_isoform_id",query_ti.id,mode="append")
            else:
                assert accept_new
                shared_tts=SharedTTS.from_TTS(tts)
                shared_tts.genomic_window_start=shared_tts.genomic_center+1-24
                shared_tts.genomic_window_end_1based=shared_tts.genomic_center+1+24

                shared_tts=update_shared_entity_type_id(shared_tts,type_id_map)
                shared_tts_by_id[shared_tts.id]=shared_tts
                shared_tts_to_id[shared_tts]=shared_tts.id

                tss_tts_itree_by_chr["tts"][(shared_tts.chrom,shared_tts.strand)].insert(
                    shared_tts.genomic_window_start,
                    shared_tts.genomic_window_end_1based,
                    shared_tts
                )

                shared_tts_to_update=[shared_tts]

            for shared_tts in shared_tts_to_update:
                ## updating "supporting_query_isoform_id"
                shared_tts.update_associated_data("supporting_query_isoform_id",query_ti.id,mode="append")
                ## updating "associated_gene_id" and feature_sets_by_gene
                shared_tts.update_associated_data("associated_gene_id",(query_ti.id,query_cls_rec["associated_gene"]),mode="append")
                if query_cls_rec["associated_gene"] in tss_tts_by_gene["tts"]:
                    tss_tts_by_gene["tts"][query_cls_rec["associated_gene"]].add(shared_tts)
                else:
                    tss_tts_by_gene["tts"][query_cls_rec["associated_gene"]]=set()
                    tss_tts_by_gene["tts"][query_cls_rec["associated_gene"]].add(shared_tts)
                if shared_tts.id.startswith("DTTS__"):
                    ## updating "associated_query_id" for novel ones
                    shared_tts.update_associated_data("associated_query_id",tts.id,mode="append")

    def reset_supporting_query_isoform_id(
        self
    ):
        self.has_parsed_query_isoforms=False
        if self.has_ref_collection:
            self.ref_collections.reset_supporting_query_isoform_id()

    def reset_differential_analysis(self,ref_collections,gffcmp_transcripts,suppa2_events,chained_transcripts):
        if ref_collections:
            self.ref_collections.reset_differential_analysis()
        if suppa2_events:
            self.suppa2_event_collection.reset_differential_analysis()
        if gffcmp_transcripts:
            self.gffcmp_merged_transcripts.reset_differential_analysis()
        if chained_transcripts:
            self.chained_transcripts.reset_differential_analysis()
        self.has_differential_analysis=False

    def parse_query_transcript_isoform(
        self,
        parse_exon=True,
        parse_junction=True,
        parse_tss_tts=True,
    ):
        assert self.has_class_df and self.has_collate_info
        def get_supporting_molecule_count_by_cell_type(barcode_umi_pair):
            mol_count_by_cell_type=defaultdict(int)
            mol_count_by_cell_type_by_sample=defaultdict(int)
            for barcode,molecule_id in barcode_umi_pair:
                sample_name=molecule_id.split('/')[0]
                if barcode in self.cell_barcodes_illumina_intersect:
                    cell_type=self.illumina_sample_info.illumina_barcodes_indexed_table_pacbio_intersect.loc[barcode,"ClusterMidway2"]
                elif barcode in self.cell_barcodes_pass_UMI_count_filter:
                    cell_type="_PB_CB_not_found_in_IL"
                else:
                    cell_type="_PB_CB_filtered"
                mol_count_by_cell_type[cell_type]+=1
                mol_count_by_cell_type_by_sample[(sample_name,cell_type)]+=1
                
            mol_count_by_cell_type=dict(mol_count_by_cell_type)
            mol_count_by_cell_type_by_sample=dict(mol_count_by_cell_type_by_sample)
            return mol_count_by_cell_type,mol_count_by_cell_type_by_sample

        def map_cell_type(cb,with_sample_name=False):
            if cb in self.cell_barcodes_illumina_intersect:
                cell_type=self.illumina_sample_info.illumina_barcodes_indexed_table_pacbio_intersect.loc[cb,"ClusterMidway2"]
            elif cb in self.cell_barcodes_pass_UMI_count_filter:
                cell_type="_PB_CB_not_found_in_IL"
            else:
                cell_type="_PB_CB_filtered"
            if with_sample_name:
                sample_index=int(cb.split('-')[1])
                sample_name=self.illumina_sample_info.pacbio_sample_index_revmap[sample_index]
                return (sample_name,cell_type)
            else:
                return cell_type

        def construct_query_transcript_isoform(genePred_record: genePredRecord, query_cls_record: pd.Series):
            query_ti=QueryTranscriptIsoform.from_genePred_record(genePred_record, query_cls_record)
            
            barcodes_all=self.isoform_molecule_barcode_map[genePred_record.id]["BCrev"]

            supporting_molecule_id=self.isoform_molecule_barcode_map[genePred_record.id]["molecule_id"]
            supporting_molecule_count=len(self.isoform_molecule_barcode_map[genePred_record.id]["molecule_id"])

            support_info=SupportInfo()
            support_info.add_CB_UMI_pair(
                self.isoform_molecule_barcode_map[genePred_record.id]["BCrev_molecule_id"]
            )
            
            counts_by_cell_type=support_info.cell_identity_stats(
                lambda cb: map_cell_type(cb, with_sample_name=False)
            )
            counts_by_cell_type_by_sample=support_info.cell_identity_stats(
                lambda cb: map_cell_type(cb, with_sample_name=True)
            )

            counts1=sum(val[1] for val in counts_by_cell_type.values())
            counts2=sum(val[1] for val in counts_by_cell_type_by_sample.values())
            assert counts1==supporting_molecule_count, "{} != {}".format(counts1,supporting_molecule_count)
            assert counts2==supporting_molecule_count, "{} != {}".format(counts2,supporting_molecule_count)

            query_ti.update_associated_data(
                "support_info",
                support_info,
                mode="assign"
            )

            query_ti.update_associated_data(
                "supporting_counts_by_cell_type",
                counts_by_cell_type,
                mode="assign"
            )

            query_ti.update_associated_data(
                "supporting_counts_by_cell_type_by_sample",
                counts_by_cell_type_by_sample,
                mode="assign"
            )

            return query_ti
        self.reset_supporting_query_isoform_id()
        from utils.file_utils import line_count
        type_id_map={
            "TS":"DTS",
            "EXON":"DEXON",
            "JUNC":"DJUNC",
            "TSS":"DTSS",
            "TTS":"DTTS",
            "SD":"DSD",
            "SA":"DSA"
        }
        if self.has_ref_collection:
            self.ref_collections.reset_supporting_query_isoform_id()
        filtered_lite_classification=self.class_df_all.set_index("pbid")
        for sample_name in self.sample_names:
            print(f"Parsing Isoforms in {sample_name}")
            genePred_fp=self.sample_analysis_fp[sample_name]["genePred"]
            
            for i,r in enumerate(tqdm(genePredReader(genePred_fp),total=line_count(genePred_fp))):
                r.id=f"QTISO__{sample_name}__{r.id}"
                if r.id in filtered_lite_classification.index and (not self.test_mode or r.chrom==TEST_CHROM):

                    query_cls_rec=filtered_lite_classification.loc[r.id]
                    query_ti=construct_query_transcript_isoform(r,query_cls_rec)
                    
                    self.query_iso_by_id[query_ti.id]=query_ti

                    if self.has_ref_collection and query_cls_rec["structural_category"] in ("full-splice_match","incomplete-splice_match","novel_in_catalog","novel_not_in_catalog"):
                        if parse_exon:
                            self._update_ref_exons_collections(
                                self.ref_collections.transcriptomic_feature_by_id["exons"],self.ref_collections.transcriptomic_feature_to_id["exons"],self.ref_collections.transcriptomic_feature_sets_by_gene,
                                query_ti,query_cls_rec,
                                type_id_map
                            )
                        if parse_junction:
                            self._update_ref_junctions_collections(
                                self.ref_collections.transcriptomic_feature_by_id["junctions"],self.ref_collections.transcriptomic_feature_to_id["junctions"],
                                self.ref_collections.transcriptomic_feature_by_id["donors"],self.ref_collections.transcriptomic_feature_to_id["donors"],
                                self.ref_collections.transcriptomic_feature_by_id["acceptors"],self.ref_collections.transcriptomic_feature_to_id["acceptors"],
                                self.ref_collections.transcriptomic_feature_sets_by_gene,
                                query_ti,query_cls_rec,
                                type_id_map
                            )
                        if parse_tss_tts:
                            self._update_ref_tss_tts_collections(
                                self.ref_collections.TSS_TTS_feature_by_id["tss"], self.ref_collections.TSS_TTS_feature_to_id["tss"],
                                self.ref_collections.TSS_TTS_feature_by_id["tts"], self.ref_collections.TSS_TTS_feature_to_id["tts"],
                                self.ref_collections.TSS_TTS_feature_itree_by_chr,self.ref_collections.TSS_TTS_by_gene,
                                query_ti,query_cls_rec,
                                type_id_map
                            )
        if self.has_ref_collection and parse_tss_tts:
            for sample_name in self.sample_names:
                print(f"Second Round Parsing Isoforms in {sample_name} (For TSS/TTS)")
                for query_ti_id, query_ti in tqdm(self.query_iso_by_id.items(),total=len(self.query_iso_by_id)):
                    if self.test_mode and query_ti.chrom!=TEST_CHROM:
                        continue
                    query_cls_rec=query_ti.cls_record
                    if query_cls_rec["structural_category"] in ("full-splice_match","incomplete-splice_match","novel_in_catalog","novel_not_in_catalog"):
                        self._update_ref_tss_tts_collections(
                            self.ref_collections.TSS_TTS_feature_by_id["tss"], self.ref_collections.TSS_TTS_feature_to_id["tss"],
                            self.ref_collections.TSS_TTS_feature_by_id["tts"], self.ref_collections.TSS_TTS_feature_to_id["tts"],
                            self.ref_collections.TSS_TTS_feature_itree_by_chr,self.ref_collections.TSS_TTS_by_gene,
                            query_ti,query_cls_rec,
                            type_id_map,
                            accept_new=False
                        )
        self.has_parsed_query_isoforms=True

    def differential_analysis(
        self,
        ref_collections=True,
        gffcmp_transcripts=True,
        chained_transcripts=True,
        suppa2_events=True,
    ):
        self.reset_differential_analysis(
            ref_collections=ref_collections,
            gffcmp_transcripts=gffcmp_transcripts,
            suppa2_events=suppa2_events,
            chained_transcripts=chained_transcripts
        )
        if ref_collections and self.has_ref_collection:
            for feature_name in ("exons","junctions","donors","acceptors","tss","tts"):
                print(f"Processing {feature_name} features")
                for feature_id,feature in tqdm(self.ref_collections._get_feature_collection(feature_name).items()):
                    if "supporting_query_isoform_id" in feature.associated_data:
                        support_info=SupportInfo()
                        for iso_id in feature.associated_data.get("supporting_query_isoform_id"):
                            support_info.add_isoform_support(self.query_iso_by_id[iso_id].associated_data.get("support_info",None))
                        supporting_pacbio_cell_barcodes_illumina_intersect=support_info.get_supporting_CBs(
                            legal_CB_set=self.cell_barcodes_illumina_intersect
                        )
                        ctable,oddsr,pvalue=self.fisher_test(supporting_pacbio_cell_barcodes_illumina_intersect)
                        cell_type_stats=self.cell_type_stats(supporting_pacbio_cell_barcodes_illumina_intersect)
                        feature.update_associated_data(
                            "differential_analysis",
                            {
                                "support_info":support_info,
                                "cell_type_stats":cell_type_stats,
                                "supporting_illumina_cell_barcodes":supporting_pacbio_cell_barcodes_illumina_intersect,
                                "fisher_exact":{
                                    "ctable":ctable,
                                    "oddsr":oddsr,
                                    "pvalue":pvalue,
                                }
                            },
                            mode="assign"
                        )
        if gffcmp_transcripts and self.has_gffcmp_merged_transcripts:
            print("Processing gffcmp transcripts")
            for feature_id,feature in tqdm(self.gffcmp_merged_transcripts.shared_iso_by_id.items()):
                if "supporting_query_isoform_id" in feature.associated_data:
                    # supporting_pacbio_cell_barcodes_pass_UMI_count_filter=set()
                    # supporting_illumina_cell_barcodes=set()
                    # supporting_molecule_count_by_cell_type=defaultdict(int)
                    # supporting_molecule_count_by_cell_type_by_sample=defaultdict(int)
                    support_info=SupportInfo()
                    supporting_query_isoform_ids=feature.associated_data.get("supporting_query_isoform_id")
                    for iso_id in supporting_query_isoform_ids:
                        support_info.add_isoform_support(self.query_iso_by_id[iso_id].associated_data.get("support_info",None))
                    # supporting_pacbio_cell_barcodes_pass_UMI_count_filter=support_info.get_supporting_CBs(self.cell_barcodes_pass_UMI_count_filter)
                    
                    # supporting_pacbio_cell_barcodes_illumina_intersect=support_info.get_supporting_CBs(self.cell_barcodes_illumina_intersect)

                    # # aggregated support counts
                    # supporting_molecule_count_by_cell_type=defaultdict(int)
                    # supporting_molecule_count_by_cell_type_by_sample=defaultdict(int)

                    # for iso_id in supporting_query_isoform_ids:
                    #     molc_ct=self.query_iso_by_id[iso_id].associated_data.get("supporting_counts_by_cell_type",None)
                    #     for cell_type, count in molc_ct.items():
                    #         supporting_molecule_count_by_cell_type[cell_type]+=count[1]
                        
                    #     molc_ct_sample=self.query_iso_by_id[iso_id].associated_data.get("supporting_counts_by_cell_type_by_sample",None)
                    #     for tp, count in molc_ct_sample.items():
                    #         supporting_molecule_count_by_cell_type_by_sample[tp]+=count[1]

                    # ctable,oddsr,pvalue=self.fisher_test(supporting_pacbio_cell_barcodes_illumina_intersect)
                    # cell_type_stats=self.cell_type_stats(supporting_pacbio_cell_barcodes_illumina_intersect)
                    feature.update_associated_data(
                        "differential_analysis",
                        {
                            "support_info":support_info,
                        },
                        mode="assign"
                    )
        if chained_transcripts and self.has_chained_transcripts:
            print("Processing gffcmp transcripts")
            for feature_id,feature in tqdm(self.chained_transcripts.shared_iso_by_id.items()):
                if "supporting_query_isoform_id" in feature.associated_data:
                    support_info=SupportInfo()
                    supporting_query_isoform_ids=feature.associated_data.get("supporting_query_isoform_id")
                    for iso_id in supporting_query_isoform_ids:
                        support_info.add_isoform_support(self.query_iso_by_id[iso_id].associated_data.get("support_info",None))
                        feature.update_associated_data(
                        "differential_analysis",
                        {
                            "support_info":support_info,
                        },
                        mode="assign"
                    )
        if suppa2_events and self.has_suppa2_event_collection:
            for feature_name in ("A3","A5","AF","AL","MX","RI","SE"):
                print(f"Processing {feature_name} features")
                for feature_id,feature in tqdm(self.suppa2_event_collection._get_feature_collection(feature_name).items()):
                    if "supporting_query_isoform_id_pos" in feature.associated_data and \
                        "supporting_query_isoform_id_neg" in feature.associated_data:
                        support_info_pos=SupportInfo()
                        support_info_neg=SupportInfo()
                        supporting_query_isoform_ids_pos=feature.associated_data.get("supporting_query_isoform_id_pos",None)
                        supporting_query_isoform_ids_neg=feature.associated_data.get("supporting_query_isoform_id_neg",None)
                        for iso_id in supporting_query_isoform_ids_pos:
                            support_info_pos.add_isoform_support(self.query_iso_by_id[iso_id].associated_data.get("support_info",None))
                        for iso_id in supporting_query_isoform_ids_neg:
                            support_info_neg.add_isoform_support(self.query_iso_by_id[iso_id].associated_data.get("support_info",None))
                        supporting_cell_barcodes_pos=support_info_pos.get_supporting_CBs(legal_CB_set=self.cell_barcodes_illumina_intersect)
                        supporting_cell_barcodes_neg=support_info_neg.get_supporting_CBs(legal_CB_set=self.cell_barcodes_illumina_intersect)
                        if False:
                            supporting_cell_barcodes_pos=set()
                            for iso_id in feature.associated_data.get("supporting_query_isoform_id_pos",None):
                                supporting_cell_barcodes_pos.update(self.query_iso_by_id[iso_id].associated_data.get("supporting_cell_barcodes_illumina_intersect",None))
                            supporting_cell_barcodes_neg=set() # TODO: may have overlapping barcodes with the positive set, need to modify
                            for iso_id in feature.associated_data.get("supporting_query_isoform_id_neg",None):
                                supporting_cell_barcodes_neg.update(self.query_iso_by_id[iso_id].associated_data.get("supporting_cell_barcodes_illumina_intersect",None))
                        ctable,oddsr,pvalue=self.fisher_test(supporting_cell_barcodes_pos,supporting_cell_barcodes_neg)
                        cell_type_stats=self.cell_type_stats(supporting_cell_barcodes_pos,supporting_cell_barcodes_neg)

                        feature.update_associated_data(
                            "differential_analysis",
                            {
                                "cell_type_stats":cell_type_stats,
                                "support_info_pos":support_info_pos,
                                "support_info_neg":support_info_neg,
                                "supporting_illumina_cell_barcodes_pos":supporting_cell_barcodes_pos,
                                "supporting_cell_barcodes_neg":supporting_cell_barcodes_neg,
                                "fisher_exact":{
                                    "ctable":ctable,
                                    "oddsr":oddsr,
                                    "pvalue":pvalue,
                                }
                            },
                            mode="assign"
                        )
        self.has_differential_analysis=True

    def fisher_test(self,positive_cell_barcodes, negative_cell_barcodes=None, illumina_intersect=False):
        if illumina_intersect:
            positive_cell_barcodes=self.cell_barcodes_illumina_intersect.intersection(positive_cell_barcodes)
            if negative_cell_barcodes:
                negative_cell_barcodes=self.cell_barcodes_illumina_intersect.intersection(negative_cell_barcodes)
        epi_normal_set=self.illumina_sample_info.cluster_midway_sets_pacbio_intersect["Epi_normal"]
        epi_tumor_set=self.illumina_sample_info.cluster_midway_sets_pacbio_intersect["Epi_tumor"]
        ctable=np.zeros((2,2),dtype=np.int64)
        n_tumor_in_set=0
        n_normal_in_set=0
        for bc in positive_cell_barcodes:
            if bc in epi_normal_set:
                n_normal_in_set+=1
            if bc in epi_tumor_set:
                n_tumor_in_set+=1
        if negative_cell_barcodes==None:
            n_normal_not_in_set=len(epi_normal_set)-n_normal_in_set
            n_tumor_not_in_set=len(epi_tumor_set)-n_tumor_in_set
        else:
            n_normal_not_in_set=0
            n_tumor_not_in_set=0
            for bc in negative_cell_barcodes:
                if bc in epi_normal_set:
                    n_normal_not_in_set+=1
                if bc in epi_tumor_set:
                    n_tumor_not_in_set+=1
        ctable[0,0]=n_tumor_in_set
        ctable[0,1]=n_tumor_not_in_set
        ctable[1,0]=n_normal_in_set
        ctable[1,1]=n_normal_not_in_set
        oddsr, p = fisher_exact(ctable, alternative='two-sided')
        return ctable.flatten(), oddsr, p

    def cell_type_stats(self,positive_cell_barcodes, negative_cell_barcodes=None):
        positive_cell_barcodes=set(positive_cell_barcodes)
        cell_type_stats=dict()
        for cell_type, all_barcodes in self.illumina_sample_info.cluster_midway_sets_pacbio_intersect.items():
            n_in_set=len(positive_cell_barcodes.intersection(all_barcodes))
            if negative_cell_barcodes is None:
                n_total=len(all_barcodes)
                n_not_in_set=n_total-n_in_set
            else:
                n_not_in_set=len(negative_cell_barcodes.intersection(all_barcodes))
            cell_type_stats[cell_type]=(n_in_set,n_not_in_set)
        return cell_type_stats

    def print_shared_features(self,*args,**kwargs):
        if self.has_ref_collection:
            self.ref_collections.print_shared_features(*args,**kwargs)
        if self.has_suppa2_event_collection:
            self.suppa2_event_collection.print_shared_features(*args,**kwargs)
        if self.has_gffcmp_merged_transcripts:
            self.gffcmp_merged_transcripts.print_shared_features(*args,**kwargs)
        if self.has_chained_transcripts:
            self.chained_transcripts.print_shared_features(*args,**kwargs)

class SupportInfo(object):
    def __init__(self):
        self.__init_collections()
    def __init_collections(self):
        self._supporting_CB_UMI_counts: DefaultDict[str,int]=defaultdict(int)

    def add_CB_UMI_pair(self, CB_UMI_pairs):
        CB_UMI_pairs=set(CB_UMI_pairs)
        for pair in CB_UMI_pairs:
            cb, umi = pair
            self._supporting_CB_UMI_counts[cb]+=1

    def add_isoform_support(self,isoform_support: "SupportInfo"):
        for cb, umi_count in isoform_support._supporting_CB_UMI_counts.items():
            self._supporting_CB_UMI_counts[cb]+=umi_count

    def get_supporting_CB_UMI_counts(self,log1p=False) -> Dict[str,int]:
        if log1p:
            return {k:np.log1p(v) for k,v in self._supporting_CB_UMI_counts.items()}
        else:
            return dict(self._supporting_CB_UMI_counts)

    def get_supporting_CB_counts(self,legal_CB_set=None):
        supporting_CBs=self._supporting_CB_UMI_counts.keys()
        if legal_CB_set is not None:
            supporting_CBs=legal_CB_set.intersection(supporting_CBs)
        return len(supporting_CBs)

    def get_supporting_CBs(self,legal_CB_set=None):
        supporting_CBs=self._supporting_CB_UMI_counts.keys()
        if legal_CB_set is not None:
            supporting_CBs=legal_CB_set.intersection(supporting_CBs)
        return supporting_CBs

    def get_supporting_UMI_counts(self,legal_CB_set=None):
        supporting_CB_UMI_counts=self._supporting_CB_UMI_counts
        if legal_CB_set is not None:
            supporting_CB_UMI_counts=sum(supporting_CB_UMI_counts[k] for k in legal_CB_set)
            return supporting_CB_UMI_counts
        else:
            return sum(supporting_CB_UMI_counts.values())

    def cell_identity_stats(self,cell_identity_map):
        cell_identity_stats=defaultdict(lambda: [0,0])
        for cb, umi_count in self._supporting_CB_UMI_counts.items():
            c_identity=cell_identity_map(cb)
            cell_identity_stats[c_identity][0]+=1
            cell_identity_stats[c_identity][1]+=umi_count
        return dict(cell_identity_stats)

class IlluminaSampleInfo(object):
    def __init__(self):
        illumina_cr_aggr_run_root="/datawaha_sfb/liz0f/Long-read-RNA-seq/Illumina/cr_aggr/cr_aggr-Colon_PS017-PS033"
        # (Illumina sample level info) Load sample name mapping between PacBio and Illumina
        sample_name_mapping_fp="PacBio/sample_table/PacBio_Illumina_mapping.csv"
        sample_name_mapping=pd.read_csv(sample_name_mapping_fp)
        sample_name_mapping.columns=["pacbio_sample_name","illumina_sample_name"]
        self.pacbio_illumina_sample_name_mapping=sample_name_mapping

        # (Illumina sample level info) Load Illumina sample aggregation table
        aggr_fp=os.path.join(illumina_cr_aggr_run_root,"outs/aggregation.csv")
        aggr_table=pd.read_csv(aggr_fp)
        aggr_table.columns=["illumina_sample_name","molecule_h5_fp","condition"]
        aggr_table=pd.merge(aggr_table,sample_name_mapping,on="illumina_sample_name",how="left")
        aggr_table["illumina_sample_index_0based"]=aggr_table.index
        aggr_table["illumina_sample_index_1based"]=aggr_table.index+1
        self.illumina_aggr_table=aggr_table

        illumina_aggr_table_pb_sub=self.illumina_aggr_table.dropna(subset=["pacbio_sample_name"])
        self.pacbio_sample_index_map=dict(zip(illumina_aggr_table_pb_sub.pacbio_sample_name,illumina_aggr_table_pb_sub.illumina_sample_index_1based))
        self.pacbio_sample_index_revmap=dict(zip(illumina_aggr_table_pb_sub.illumina_sample_index_1based,illumina_aggr_table_pb_sub.pacbio_sample_name))
        self.illumina_sample_index_map=dict(zip(self.illumina_aggr_table.illumina_sample_name,self.illumina_aggr_table.illumina_sample_index_1based))
        self.illumina_sample_index_revmap=dict(zip(self.illumina_aggr_table.illumina_sample_index_1based,self.illumina_aggr_table.illumina_sample_name))
        self.illumina_sample_index2condition=dict(zip(self.illumina_aggr_table.illumina_sample_index_1based,self.illumina_aggr_table.condition))
        # anndata_fp="/datawaha_sfb/liz0f/Long-read-RNA-seq/downstream/Colon_PS017-PS033-bc295-combined_visualization-ClusterFull-AllCellTypes/hdf5/expr.clusterfull.h5ad"
        # anndata_fp = "/data/liz0f/Long-read-RNA-seq/Illumina/expr.clusterfull.h5ad"
        anndata_fp = "/datawaha_sfb/liz0f/Long-read-RNA-seq/downstream/Colon_PS017-PS033-bc295-17976_features-ClusterFull-AllCellTypes/hdf5/expr.clusterfull.enterocyte_split.h5ad"
        self.anndata=anndata.read_h5ad(anndata_fp)
        illumina_barcodes_table=self.anndata.obs.copy()

        illumina_barcodes_table["ClusterFull"]=illumina_barcodes_table["ClusterFull"].map(lambda x: {
            'In-house-N cE01 (Stem/TA-like)':'Epi_normal_cE01_stem_TA',
            'In-house-N cE02 (Stem/TA-like/Immature Goblet)': 'Epi_normal_cE02_stem_TA_im_goblet',
            'In-house-N cE03 (Stem/TA-like prolif)': 'Epi_normal_cE03_stem_TA_prolif',
            'In-house-N cE04 (Enterocyte 1)': 'Epi_normal_cE04_enterocyte1',
            'In-house-N cE05 (Enterocyte 2)': 'Epi_normal_cE05_enterocyte2',
            'In-house-N cE06 (Immature Goblet)': 'Epi_normal_cE06_immature_goblet',
            'In-house-N cE07 (Goblet/Enterocyte)': 'Epi_normal_cE07_goblet_enterocyte',
            'In-house-N cE08 (Goblet)': 'Epi_normal_cE08_goblet',
            'In-house-N cE09 (Best4)': 'Epi_normal_cE09_best4',
            'In-house-N cE10 (Tuft)': 'Epi_normal_cE10_tuft',
            'In-house-N cE11 (Enteroendocrine)': 'Epi_normal_cE11_enteroendocrine',
            'In-house-T cE01 (Stem/TA-like)':'Epi_tumor_cE01_stem_TA',
            'In-house-T cE02 (Stem/TA-like/Immature Goblet)': 'Epi_tumor_cE02_stem_TA_im_goblet',
            'In-house-T cE03 (Stem/TA-like prolif)': 'Epi_tumor_cE03_stem_TA_prolif',
            'In-house-T cE04 (Enterocyte 1)': 'Epi_tumor_cE04_enterocyte1',
            'In-house-T cE05 (Enterocyte 2)': 'Epi_tumor_cE05_enterocyte2',
            'In-house-T cE06 (Immature Goblet)': 'Epi_tumor_cE06_im_goblet',
            'In-house-T cE07 (Goblet/Enterocyte)': 'Epi_tumor_cE07_goblet_enterocyte',
            'In-house-T cE08 (Goblet)': 'Epi_tumor_cE08_goblet',
            'In-house-T cE09 (Best4)': 'Epi_tumor_cE09_best4',
            'In-house-T cE10 (Tuft)': 'Epi_tumor_cE10_tuft',
            'In-house-T cE11 (Enteroendocrine)': 'Epi_tumor_cE11_enteroendocrine'
        }.get(x,x))

        illumina_barcodes_table["ClusterFull2"]=illumina_barcodes_table["ClusterFull2"].map(lambda x: {
            'In-house-N cE01 (Stem/TA-like)':'Epi_normal_cE01_stem_TA',
            'In-house-N cE02 (Stem/TA-like/Immature Goblet)': 'Epi_normal_cE02_stem_TA_im_goblet',
            'In-house-N cE03 (Stem/TA-like prolif)': 'Epi_normal_cE03_stem_TA_prolif',
            'In-house-N cE04 (Enterocyte 1)': 'Epi_normal_cE04_enterocyte1',
            'In-house-N cE05a (Enterocyte 2)': 'Epi_normal_cE05a_enterocyte2',
            'In-house-N cE05b (Stem/TA-like/Enterocyte)':'Epi_normal_cE05b_stem_TA_enterocyte',
            'In-house-N cE06 (Immature Goblet)': 'Epi_normal_cE06_im_goblet',
            'In-house-N cE07 (Goblet/Enterocyte)': 'Epi_normal_cE07_goblet_enterocyte',
            'In-house-N cE08 (Goblet)': 'Epi_normal_cE08_goblet',
            'In-house-N cE09 (Best4)': 'Epi_normal_cE09_best4',
            'In-house-N cE10 (Tuft)': 'Epi_normal_cE10_tuft',
            'In-house-N cE11 (Enteroendocrine)': 'Epi_normal_cE11_enteroendocrine',
            'In-house-T cE01 (Stem/TA-like)':'Epi_tumor_cE01_stem_TA',
            'In-house-T cE02 (Stem/TA-like/Immature Goblet)': 'Epi_tumor_cE02_stem_TA_im_goblet',
            'In-house-T cE03 (Stem/TA-like prolif)': 'Epi_tumor_cE03_stem_TA_prolif',
            'In-house-T cE04 (Enterocyte 1)': 'Epi_tumor_cE04_enterocyte1',
            'In-house-T cE05 (Enterocyte 2)': 'Epi_tumor_cE05_enterocyte2',
            'In-house-T cE06 (Immature Goblet)': 'Epi_tumor_cE06_im_goblet',
            'In-house-T cE07 (Goblet/Enterocyte)': 'Epi_tumor_cE07_goblet_enterocyte',
            'In-house-T cE08 (Goblet)': 'Epi_tumor_cE08_goblet',
            'In-house-T cE09 (Best4)': 'Epi_tumor_cE09_best4',
            'In-house-T cE10 (Tuft)': 'Epi_tumor_cE10_tuft',
            'In-house-T cE11 (Enteroendocrine)': 'Epi_tumor_cE11_enteroendocrine'
        }.get(x,x))

        illumina_barcodes_table.index.name="illumina_indexed_barcode"
        self.illumina_barcodes_indexed_table=illumina_barcodes_table
        self.illumina_barcodes_table=illumina_barcodes_table.reset_index()
        self.illumina_barcodes=set(self.illumina_barcodes_table["illumina_indexed_barcode"].values)
        self.get_cell_type_sets()
        self.get_all_per_gene_stats()

    @property
    def gene_symbols(self):
        return self.anndata.var.index

    def get_cell_type_sets(self):
        cell_barcodes_classified=self.illumina_barcodes_table.groupby("ClusterMidway2")["illumina_indexed_barcode"].unique()
        cluster_midway_sets={k:set(cell_barcodes_classified[k]) for k in cell_barcodes_classified.index}
        self.cluster_midway_sets=cluster_midway_sets
        self.cluster_midway_stats={k:len(self.cluster_midway_sets[k]) for k in self.cluster_midway_sets}

        cell_barcodes_classified=self.illumina_barcodes_table.groupby("ClusterFull2")["illumina_indexed_barcode"].unique()
        cluster_full_sets={k:set(cell_barcodes_classified[k]) for k in cell_barcodes_classified.index if k not in ["UNK"]}
        self.cluster_full_sets=cluster_full_sets
        self.cluster_full_stats={k:len(self.cluster_full_sets[k]) for k in self.cluster_full_sets}

        cell_barcodes_classified=self.illumina_barcodes_table.groupby("condition")["illumina_indexed_barcode"].unique()
        sample_condition_sets={k:set(cell_barcodes_classified[k]) for k in cell_barcodes_classified.index}
        self.sample_condition_sets=sample_condition_sets
        self.sample_condition_stats={k:len(self.sample_condition_sets[k]) for k in self.sample_condition_sets}

    def get_pacbio_intersection(self,pb_barcodes):
        # compute barcode intersection
        filter_stat=defaultdict(lambda:[0,0])
        cell_barcodes_pass_filter=list()
        if self.illumina_barcodes:
            for barcode in pb_barcodes:
                bc,sample_index=barcode.split("-")
                sample_index=int(sample_index)
                filter_stat[sample_index][0]+=1
                if barcode in self.illumina_barcodes:
                    filter_stat[sample_index][1]+=1
                    cell_barcodes_pass_filter.append(barcode)
        self.cell_barcodes_pacbio_intersect=set(cell_barcodes_pass_filter)

        # compute cell type sets intersection
        self.cluster_midway_sets_pacbio_intersect={k:v.intersection(self.cell_barcodes_pacbio_intersect) for k,v in self.cluster_midway_sets.items()}
        self.cluster_midway_stats_pacbio_intersect={k:len(self.cluster_midway_sets_pacbio_intersect[k]) for k in self.cluster_midway_sets_pacbio_intersect}

        self.cluster_full_sets_pacbio_intersect={k:v.intersection(self.cell_barcodes_pacbio_intersect) for k,v in self.cluster_full_sets.items()}
        self.cluster_full_stats_pacbio_intersect={k:len(self.cluster_full_sets_pacbio_intersect[k]) for k in self.cluster_full_sets_pacbio_intersect}
        
        self.illumina_barcodes_indexed_table_pacbio_intersect=self.illumina_barcodes_indexed_table[
            self.illumina_barcodes_indexed_table.index.isin(self.cell_barcodes_pacbio_intersect)
        ]

        self.illumina_barcodes_table_pacbio_intersect=self.illumina_barcodes_table[
            self.illumina_barcodes_table["illumina_indexed_barcode"].isin(self.cell_barcodes_pacbio_intersect)
        ].copy()

        self.illumina_barcodes_table_pacbio_intersect["PB_sample_name"]=self.illumina_barcodes_table_pacbio_intersect\
            ["illumina_indexed_barcode"].map(lambda x: self.pacbio_sample_index_revmap[int(x.split("-")[1])])

        self.illumina_barcodes_cluster_midway_sn_map_pacbio_intersect=dict(
            zip(
                self.illumina_barcodes_table_pacbio_intersect["illumina_indexed_barcode"],
                zip(self.illumina_barcodes_table_pacbio_intersect["PB_sample_name"],self.illumina_barcodes_table_pacbio_intersect["ClusterMidway2"])
            )
        )

        self.illumina_barcodes_cluster_full_sn_map_pacbio_intersect=dict(
            zip(
                self.illumina_barcodes_table_pacbio_intersect["illumina_indexed_barcode"],
                zip(self.illumina_barcodes_table_pacbio_intersect["PB_sample_name"],self.illumina_barcodes_table_pacbio_intersect["ClusterFull2"])
            )
        )
        return filter_stat

    def get_per_gene_stats2(self,gene_symbol):
        assert gene_symbol in self.gene_symbols
        data=defaultdict(list)
        for cell_type in chain(self.cluster_midway_sets.keys(),["ALL"]):
            if cell_type=="ALL":
                subset_index=self.anndata.obs.index
            else:
                subset_index=self.anndata.obs.index[self.anndata.obs["ClusterMidway2"]==cell_type]
            X_subset=self.anndata[subset_index,gene_symbol].X
            X_subset_nonzero=X_subset[X_subset>0]
            n_barcodes_nonzero=X_subset_nonzero.shape[1]
            
            if X_subset.shape[1]>0:
                UMI_mean=X_subset.mean()
            else:
                UMI_mean=0

            if n_barcodes_nonzero>0:
                UMI_mean_nonzero=X_subset_nonzero.mean()
            else:
                UMI_mean_nonzero=0
            
            data["cell_type"].append(cell_type)
            data["n_barcodes"].append(len(subset_index))
            data["n_barcodes_nonzero"].append(n_barcodes_nonzero)
            data["UMI_mean"].append(UMI_mean)
            data["UMI_mean_nonzero"].append(UMI_mean_nonzero)

        data_df=pd.DataFrame(data)
        return data_df
    def get_per_gene_stats(self,gene_symbol):
        return self.all_per_gene_stats[gene_symbol]
    def get_all_per_gene_stats(self,log1p=True):
        subset_mean_list=list()
        subset_mean_nonzero_list=list()
        n_barcodes_list=list()
        n_barcodes_nonzero_list=list()
        cell_types_list=list(self.cluster_midway_sets.keys())+["ALL"]
        index2int=dict(zip(self.anndata.obs.index,range(len(self.anndata.obs.index))))
        for cell_type in cell_types_list:
            if cell_type=="ALL":
                subset_index=list(range(len(self.anndata.obs.index)))
            else:
                subset_index_str=self.anndata.obs.index[self.anndata.obs["ClusterMidway2"]==cell_type]
                subset_index=[index2int[s] for s in subset_index_str]
            if log1p:
                X_subset=self.anndata[subset_index,:].X
                X_subset=sp.csc_matrix(X_subset)
                X_subset=X_subset.expm1()
                X_subset_sum=np.asarray(X_subset.sum(axis=0)).squeeze()
                counts = np.diff(X_subset.indptr)
                X_subset_mean = np.log1p(X_subset_sum/X_subset.shape[0])
                X_subset_mean_nonzero = np.log1p(X_subset_sum/(counts+np.finfo(np.float64).eps))
            else:
                X_subset=self.anndata[subset_index,:].raw.X
                X_subset=sp.csc_matrix(X_subset)
                X_subset_sum=np.asarray(X_subset.sum(axis=0)).squeeze()
                counts = np.diff(X_subset.indptr)
                X_subset_mean = X_subset_sum/X_subset.shape[0]
                X_subset_mean_nonzero = (X_subset_sum/(counts+np.finfo(np.float64).eps))
            n_barcodes_list.append(X_subset.shape[0])
            n_barcodes_nonzero_list.append(counts)
            subset_mean_list.append(X_subset_mean)
            subset_mean_nonzero_list.append(X_subset_mean_nonzero)
        subset_mean_all=np.vstack(subset_mean_list)
        subset_mean_nonzero_all=np.vstack(subset_mean_nonzero_list)
        n_barcodes_all=np.array(n_barcodes_list)
        n_barcodes_nonzero_all=np.vstack(n_barcodes_nonzero_list)

        all_per_gene_stats=dict()
        for i,gene_symbol in enumerate(self.gene_symbols):
            all_per_gene_stats[gene_symbol]=pd.DataFrame({
                "cell_type":cell_types_list,
                "n_barcodes":n_barcodes_list,
                "n_barcodes_nonzero":n_barcodes_nonzero_all[:,i],
                "UMI_mean":subset_mean_all[:,i],
                "UMI_mean_nonzero":subset_mean_nonzero_all[:,i]
            })
        self.all_per_gene_stats=all_per_gene_stats
    
        

class ReferenceSharedFeatureCollection(object):
    def __init__(
        self,
        reference_genePred_fp,
        TSS_bed_fp,
        TTS_bed_fp,
        test_mode=False,
        filter_length=True
    ):
        self.test_mode=test_mode
        # TODO: update
        print("="*SEP_SIZE+"Loading reference genePred"+"="*SEP_SIZE)
        iso_searchable_data_struct, \
        ref_iso_by_id, \
        transcriptomic_feature_by_id, transcriptomic_feature_to_id, \
        transcriptomic_feature_itree_by_chr, transcriptomic_feature_sets_by_gene= \
        reference_genePred_parser(
            reference_genePred_fp,
            test_mode_chr=TEST_CHROM if self.test_mode else None,
            filter_length=filter_length
        )

        self.ref_iso_by_id=ref_iso_by_id
        self.iso_searchable_data_struct=iso_searchable_data_struct
        self.transcriptomic_feature_by_id=transcriptomic_feature_by_id
        self.transcriptomic_feature_to_id=transcriptomic_feature_to_id
        self.transcriptomic_feature_itree_by_chr=transcriptomic_feature_itree_by_chr
        self.transcriptomic_feature_sets_by_gene=transcriptomic_feature_sets_by_gene
        
        print("="*SEP_SIZE+"Loading reference TSS&TTS"+"="*SEP_SIZE)
        TSS_TTS_feature_by_id, TSS_TTS_feature_to_id, TSS_TTS_feature_itree_by_chr, TSS_TTS_by_gene = \
        reference_TSS_TTS_parser(
            TSS_bed_fp,
            TTS_bed_fp,
            ref_iso_by_id,
            test_mode_chr=TEST_CHROM if self.test_mode else None
        )

        self.TSS_TTS_feature_by_id = TSS_TTS_feature_by_id
        self.TSS_TTS_feature_to_id = TSS_TTS_feature_to_id
        self.TSS_TTS_feature_itree_by_chr = TSS_TTS_feature_itree_by_chr
        self.TSS_TTS_by_gene = TSS_TTS_by_gene
        self.feature_type_id_map={
            "exons":("REXON","DEXON"),
            "junctions":("RJUNC","DJUNC"),
            "tss":("RTSS","DTSS"),
            "tts":("RTTS","DTTS"),
            "donors":("RSD","DSD"),
            "acceptors":("RSA","DSA")
        }
    def reset_supporting_query_isoform_id(self):
        for feature_type, feature_dict in self.transcriptomic_feature_by_id.items():
            for fid,feature in feature_dict.items():
                feature.reset_associated_data("supporting_query_isoform_id")

    def reset_differential_analysis(self):
        for feature_type, feature_dict in self.transcriptomic_feature_by_id.items():
            for fid,feature in feature_dict.items():
                feature.reset_associated_data("differential_analysis")

    def _get_feature_collection(self,feature_name):
        if feature_name in ("exons","junctions","donors","acceptors"):
            feature_collection=self.transcriptomic_feature_by_id[feature_name]
        elif feature_name in ("tss","tts"):
            feature_collection=self.TSS_TTS_feature_by_id[feature_name]
        else:
            assert False
        return feature_collection


    def print_shared_features(
        self,
        print_limit=5,
        seed=123,
        print_example=True
    ):
        np.random.seed(seed)
        
        for feature_name in ("exons","junctions","donors","acceptors","tss","tts"):
            iprint(f"*feature_name: {feature_name}")
            feature_collection=self._get_feature_collection(feature_name)
            ref_feature_ids = [k for k in feature_collection.keys() if k.split("__")[0]==self.feature_type_id_map[feature_name][0]]
            novel_feature_ids = [k for k in feature_collection.keys() if k.split("__")[0]==self.feature_type_id_map[feature_name][1]]
            iprint(f"**ALL: # reference features {len(ref_feature_ids)}, # novel features {len(novel_feature_ids)}",ident=1)
            ref_feature_ids_supported= [k for k in ref_feature_ids if "supporting_query_isoform_id" in feature_collection[k].associated_data]
            novel_feature_ids_supported= [k for k in novel_feature_ids if "supporting_query_isoform_id" in feature_collection[k].associated_data]
            iprint(f"**Supported: # reference features {len(ref_feature_ids_supported)}, # novel features {len(novel_feature_ids_supported)}",ident=1)
            if print_example:
                for feature_id in chain(
                    np.random.choice(ref_feature_ids_supported,min(print_limit,len(ref_feature_ids_supported))),
                    np.random.choice(novel_feature_ids_supported,min(print_limit,len(novel_feature_ids_supported)))
                ):
                    iprint("->",feature_id,sorted(feature_collection[feature_id].associated_data.get("supporting_query_isoform_id")),ident=1)

            ref_feature_ids_differential_analysis= [k for k in ref_feature_ids if "differential_analysis" in feature_collection[k].associated_data]
            novel_feature_ids_differential_analysis= [k for k in novel_feature_ids if "differential_analysis" in feature_collection[k].associated_data]
            iprint(f"**Differential Analysis: # reference features {len(ref_feature_ids_differential_analysis)}, # novel features {len(novel_feature_ids_differential_analysis)}",ident=1)
            if print_example:
                for feature_id in chain(
                    np.random.choice(ref_feature_ids_differential_analysis,min(print_limit,len(ref_feature_ids_differential_analysis))),
                    np.random.choice(novel_feature_ids_differential_analysis,min(print_limit,len(novel_feature_ids_differential_analysis)))
                ):
                    diff_analysis=feature_collection[feature_id].associated_data.get("differential_analysis")
                    if len(diff_analysis["supporting_illumina_cell_barcodes"])>0:
                        iprint("->",feature_id,ident=1)
                        iprint("# supporting_illumina_cell_barcodes",len(diff_analysis["supporting_illumina_cell_barcodes"]),ident=2)
                        iprint("ctable",repr(diff_analysis["fisher_exact"]["ctable"]),ident=3)
                        iprint("oddsr",diff_analysis["fisher_exact"]["oddsr"],ident=3)
                        iprint("pvalue",diff_analysis["fisher_exact"]["pvalue"],ident=3)
                        iprint("cell_type_stats",diff_analysis["cell_type_stats"],ident=3)
                        iprint()

class GffCompareMergedTranscripts(object):
    def __init__(
        self,
        gffcmp_fp,
        query_iso_by_id,
        sample_name_map,
        test_mode=False
    ):
        self.test_mode=test_mode
        self.gffcmp_fp=gffcmp_fp
        self.query_iso_by_id=query_iso_by_id
        self.sample_name_map=sample_name_map
        self.__load_gffcmp_tracking_file()
        self.__parse_transcript_isoform()

    def __load_gffcmp_tracking_file(self):
        merged_QTISO_mapping=defaultdict(set)
        tracking_table=pd.read_csv(self.gffcmp_fp["gffcmp_tracking_fp"],sep="\t",header=None)
        for _, row in tracking_table.iterrows():
            gffcmp_merged_isoid=row[0]
            original_isoform_correspondance=row[4:]
            for o_iso in original_isoform_correspondance:
                if o_iso!="-":
                    sample_number=int(o_iso.split(":")[0].replace("q",""))
                    if sample_number in self.sample_name_map:
                        o_iso_info_str=o_iso.split(":")[1]
                        for o_iso_info_str_sub in o_iso_info_str.split(','):
                            sample_name=self.sample_name_map[sample_number]
                            o_isoid="QTISO__{}__{}".format(sample_name,o_iso_info_str_sub.split("|")[1])
                            merged_QTISO_mapping[gffcmp_merged_isoid].add(o_isoid)
        self.merged_QTISO_mapping=dict(merged_QTISO_mapping) # gffcompare merged isoform id -> unmerged isoform id
    
    def __parse_transcript_isoform(self):
        def construct_shared_transcript_isoform(
            genePred_record: genePredRecord,
            cls_record: pd.Series,
            query_iso_by_id: Dict[str,QueryTranscriptIsoform],
            merged_QTISO_mapping: Dict[str,set]
        ):
            shared_ti=SharedTranscriptIsoform.from_genePred_record(genePred_record)

            shared_ti.update_associated_data(
                "cls_record",
                cls_record,
                mode="assign"
            )

            qtiso_collection=merged_QTISO_mapping[shared_ti.id]
            shared_ti.update_associated_data(
                "supporting_query_isoform_id",
                qtiso_collection,
                mode="update"
            )

            if query_iso_by_id is not None:
                shared_ti.update_associated_data(
                    "supporting_query_isoform_cls_rec",
                    {qtiso:query_iso_by_id[qtiso].cls_record 
                    for qtiso in qtiso_collection 
                    if qtiso.split("__")[1] in self.sample_name_map.values()},
                    mode="assign"
                )
                shared_ti.update_associated_data(
                    "associated_gene_id",
                    [(qtiso,query_iso_by_id[qtiso].cls_record["associated_gene"])
                        for qtiso in qtiso_collection 
                        if qtiso.split("__")[1] in self.sample_name_map.values()],
                    mode="update"
                )


                shared_ti.update_associated_data(
                    "associated_rtiso_id",
                    [
                        (qtiso,query_iso_by_id[qtiso].cls_record["associated_transcript"])
                        for qtiso in qtiso_collection
                        if qtiso.split("__")[1] in self.sample_name_map.values()
                    ],
                    mode="update"
                )
            return shared_ti
        
        def process_CDS_start_end(df):
            df=df.copy()
            swapped=df[["strand","CDS_genomic_start","CDS_genomic_end"]].apply(
                lambda x: (x["CDS_genomic_start"],x["CDS_genomic_end"]) 
                if x["strand"]=='+' else (x["CDS_genomic_end"],x["CDS_genomic_start"]), axis=1,result_type="expand"
            )
            df["CDS_genomic_start"]=swapped[0]-1 # make 0-based
            df["CDS_genomic_end"]=swapped[1]
            df["CDS_start"]-=1 # make 0-based
            return df

        class_df=pd.read_csv(self.gffcmp_fp["filtered_lite_classification"],sep="\t")
        class_df=process_CDS_start_end(class_df)
        class_df=class_df.rename({"isoform":"pbid","FL":"FL_count","CDS_genomic_end":"CDS_genomic_end_1based","CDS_end":"CDS_end_1based"},axis=1)
        class_df["pbid"]=class_df["pbid"].map(lambda x: x)
        class_df["associated_transcript"]=class_df["associated_transcript"].apply(lambda x: "RTISO__"+x if x!="novel" else "novel")
        class_df=class_df.set_index("pbid")
        shared_iso_by_id=dict()
        for i,r in enumerate(tqdm(genePredReader(self.gffcmp_fp["genePred_fp"]),total=line_count(self.gffcmp_fp["genePred_fp"]))):
            r.id = f"{r.id}"
            if r.id in self.merged_QTISO_mapping and (not self.test_mode or r.chrom==TEST_CHROM):
                cls_record=class_df.loc[r.id]
                shared_ti=construct_shared_transcript_isoform(r,cls_record,self.query_iso_by_id,self.merged_QTISO_mapping)
                shared_iso_by_id[shared_ti.id]=shared_ti
        self.shared_iso_by_id=shared_iso_by_id
        self.class_df_all=class_df

    def print_shared_features(self,print_limit=5,seed=123,print_example=True):
        feature_collection=self.shared_iso_by_id
        iprint("*feature_name: gffcmp_merged_transcripts")
        iprint(f"**ALL: # events {len(feature_collection)}",ident=1)
        feature_collection_keys=list(feature_collection.keys())
        if print_example:
            for feature_id in np.random.choice(feature_collection_keys,min(print_limit,len(feature_collection_keys))):
                iprint("->",feature_id,ident=1)
                iprint("supporting_query_isoform_id",sorted(feature_collection[feature_id].associated_data.get("supporting_query_isoform_id")),ident=2)

        feature_ids_differential_analysis= [k for k in feature_collection_keys if "differential_analysis" in feature_collection[k].associated_data]
        iprint(f"**Differential Analysis: # features {len(feature_ids_differential_analysis)}",ident=1)
        if print_example:
            for feature_id in np.random.choice(feature_ids_differential_analysis,min(print_limit,len(feature_ids_differential_analysis))):
                diff_analysis=feature_collection[feature_id].associated_data.get("differential_analysis")
                if len(diff_analysis["supporting_illumina_cell_barcodes"])>0:
                    iprint("->",feature_id,ident=1)
                    iprint("# supporting_illumina_cell_barcodes",len(diff_analysis["supporting_illumina_cell_barcodes"]),ident=2)
                    iprint("ctable",repr(diff_analysis["fisher_exact"]["ctable"]),ident=3)
                    iprint("oddsr",diff_analysis["fisher_exact"]["oddsr"],ident=3)
                    iprint("pvalue",diff_analysis["fisher_exact"]["pvalue"],ident=3)
                    iprint("cell_type_stats",diff_analysis["cell_type_stats"],ident=3)
                    iprint()

    def reset_differential_analysis(self):
        for fid,feature in self.shared_iso_by_id.items():
            feature.reset_associated_data("differential_analysis")

class ChainedTranscripts(object):
    def __init__(
        self,
        chained_fp,
        query_iso_by_id=None,
        test_mode=False
    ):
        self.test_mode=test_mode
        self.chained_fp=chained_fp
        self.query_iso_by_id=query_iso_by_id
        self.__load_chained_tracking_file()
        self.__parse_transcript_isoform()
        
    def __load_chained_tracking_file(self):
        merged_QTISO_mapping=defaultdict(set)
        tracking_table=pd.read_csv(self.chained_fp["chained_tracking_fp"],sep="\t")
        self.sample_names=tracking_table.columns[1:]
        for _, row in tracking_table.iterrows():
            chained_pbid="CHAINED__{}".format(row["superPBID"])
            for sample_name in tracking_table.columns[1:]:
                pbid=row[sample_name]
                if not pd.isna(pbid):
                    o_isoid="QTISO__{}__{}".format(sample_name,pbid)
                    merged_QTISO_mapping[chained_pbid].add(o_isoid)
        self.merged_QTISO_mapping=dict(merged_QTISO_mapping) # gffcompare merged isoform id -> unmerged isoform id
    
    def __parse_transcript_isoform(self):
        def construct_shared_transcript_isoform(
            genePred_record: genePredRecord,
            cls_record: pd.Series,
            query_iso_by_id: Dict[str,QueryTranscriptIsoform],
            merged_QTISO_mapping: Dict[str,set]
        ):
            shared_ti=SharedTranscriptIsoform.from_genePred_record(genePred_record,cls_record)
            shared_ti.update_associated_data(
                "cls_record",
                cls_record,
                mode="assign"
            )

            # add supporting QTISO information
            qtiso_collection=merged_QTISO_mapping[shared_ti.id]
            shared_ti.update_associated_data(
                "supporting_query_isoform_id",
                qtiso_collection,
                mode="update"
            )
            shared_ti.update_associated_data(
                "supporting_query_isoform_cls_rec",
                {
                    qtiso:query_iso_by_id[qtiso].cls_record 
                    for qtiso in qtiso_collection 
                    if qtiso.split("__")[1] in self.sample_names
                },
                mode="assign"
            )
            shared_ti.update_associated_data(
                "associated_gene_id",
                [
                    (qtiso,query_iso_by_id[qtiso].cls_record["associated_gene"])
                    for qtiso in qtiso_collection 
                    if qtiso.split("__")[1] in self.sample_names
                ],
                mode="update"
            )

            shared_ti.update_associated_data(
                "associated_rtiso_id",
                [
                    (qtiso,query_iso_by_id[qtiso].cls_record["associated_transcript"])
                    for qtiso in qtiso_collection
                    if qtiso.split("__")[1] in self.sample_names
                ],
                mode="update"
            )
            
            return shared_ti
        
        def process_CDS_start_end(df):
            df=df.copy()
            swapped=df[["strand","CDS_genomic_start","CDS_genomic_end"]].apply(
                lambda x: (x["CDS_genomic_start"],x["CDS_genomic_end"]) 
                if x["strand"]=='+' else (x["CDS_genomic_end"],x["CDS_genomic_start"]), axis=1,result_type="expand"
            )
            df["CDS_genomic_start"]=swapped[0]-1 # make 0-based
            df["CDS_genomic_end"]=swapped[1]
            df["CDS_start"]-=1 # make 0-based
            return df
        
        class_df=pd.read_csv(self.chained_fp["filtered_lite_classification"],sep="\t")
        class_df=process_CDS_start_end(class_df)
        class_df=class_df.rename({"isoform":"pbid","FL":"FL_count","CDS_genomic_end":"CDS_genomic_end_1based","CDS_end":"CDS_end_1based"},axis=1)
        class_df["pbid"]=class_df["pbid"].map(lambda x: f"CHAINED__{x}")
        class_df["associated_transcript"]=class_df["associated_transcript"].apply(lambda x: "RTISO__"+x if x!="novel" else "novel")
        class_df=class_df.set_index("pbid")
        shared_iso_by_id=dict()
        for i,r in enumerate(tqdm(genePredReader(self.chained_fp["genePred_fp"]),total=line_count(self.chained_fp["genePred_fp"]))):
            r.id = f"CHAINED__{r.id}"
            if r.id in self.merged_QTISO_mapping and (not self.test_mode or r.chrom==TEST_CHROM):
                cls_record=class_df.loc[r.id]
                shared_ti=construct_shared_transcript_isoform(r,cls_record,self.query_iso_by_id,self.merged_QTISO_mapping)
                shared_iso_by_id[shared_ti.id]=shared_ti
        self.shared_iso_by_id=shared_iso_by_id
        self.class_df_all=class_df

    def print_shared_features(self,print_limit=5,seed=123,print_example=False):
        feature_collection=self.shared_iso_by_id
        iprint("*feature_name: chained transcripts")
        iprint(f"**ALL: # events {len(feature_collection)}",ident=1)
        feature_collection_keys=list(feature_collection.keys())
        if print_example:
            for feature_id in np.random.choice(feature_collection_keys,min(print_limit,len(feature_collection_keys))):
                iprint("->",feature_id,ident=1)
                iprint("supporting_query_isoform_id",sorted(feature_collection[feature_id].associated_data.get("supporting_query_isoform_id")),ident=2)

        feature_ids_differential_analysis= [k for k in feature_collection_keys if "differential_analysis" in feature_collection[k].associated_data]
        iprint(f"**Differential Analysis: # features {len(feature_ids_differential_analysis)}",ident=1)

    def reset_differential_analysis(self):
        for fid,feature in self.shared_iso_by_id.items():
            feature.reset_associated_data("differential_analysis")
class PacBioIlluminaCombinedInfo(object):
    def __init__(self):
        self.pacbio_illumina_combined_dir="PacBio/downstream/combined"
        self.load_illumina_info()

    def load_illumina_info(self):
        illumina_cr_aggr_run_root="/datawaha_sfb/liz0f/Long-read-RNA-seq/Illumina/cr_aggr/cr_aggr-Colon_PS017-PS033"
        self.illumina_sample_info=IlluminaSampleInfo()

        # (Illumina cell level info) Load Illumina cell barcode and metadata
        barcodes_fp=join(self.pacbio_illumina_combined_dir,"illumina_barcodes.xgboost.infercnv.midwaycorrect.csv")
        illumina_barcodes_table=pd.read_csv(barcodes_fp)
        illumina_barcodes_table.columns=["illumina_indexed_barcode","sample_source","library_id","ClusterTop","ClusterMidway","ClusterMidway2","condition","cnv_predicted_tumor","xgboost_predicted_tumor"]

        self.illumina_barcodes_table=illumina_barcodes_table

    def load_pacbio_info(self):
        # (PacBio isoform level info) Load selected PacBio isoforms for analysis
        selected_isoforms_fp=join(self.pacbio_illumina_combined_dir,"selected_isoforms_filtered.csv")
        selected_isoforms=pd.read_csv(selected_isoforms_fp)
        selected_isoforms=selected_isoforms[[
            "pacbio_tid","chrom","strand","length","exons","structural_category","structural_category_coarse",
            "associated_transcript","associated_gene","associated_gene_symbol",
            'ref_length', 'ref_exons', 'sample_name_repr', 
            'n_umi', 'n_cells', 'n_samples', 'first_sample'
        ]]
        selected_isoforms["display_name"]=selected_isoforms.apply(
            lambda row: row["pacbio_tid"] + "." + row["associated_gene_symbol"]+"."+row["structural_category_coarse"],
            axis=1
        )
        self.selected_isoforms=selected_isoforms

        # (PacBio molecule level info) Load correspondence between PacBio isoforms and PB cell barcodes
        molecule_level_collate_df_fp=join(self.pacbio_illumina_combined_dir,"molecule_level_collate_df_filtered.csv")
        molecule_level_collate_df=pd.read_csv(molecule_level_collate_df_fp)
        molecule_level_collate_df=molecule_level_collate_df[[
            "molecule_id",'pacbio_tid_old', 'length', 'transcript', 'gene',
            'category', 'ontarget', 'ORFgroup', 'UMI', 'UMIrev', 'BC', 'BCrev',
            'sample_name', 'pacbio_tid'
        ]]
        sample_index_map=self.illumina_aggr_table.set_index("pacbio_sample_name")
        sample_index_map=dict(zip(sample_index_map.index,sample_index_map.illumina_sample_index_1based))

        molecule_level_collate_df["pacbio_indexed_barcode"] = molecule_level_collate_df.apply(
            lambda row: row["BCrev"]+'-'+str(sample_index_map[row["sample_name"]]),
            axis=1
        )
        self.molecule_level_collate_df=molecule_level_collate_df

        print(
            "collate_df_all (before filter with Illumina):",
            "shape",molecule_level_collate_df.shape,
            "n_barcodes",molecule_level_collate_df["BCrev"].nunique(),
            "n_pacbio_tid",molecule_level_collate_df["pacbio_tid"].nunique()
        )
        molecule_level_collate_df_il_itersect=molecule_level_collate_df[
            molecule_level_collate_df["pacbio_indexed_barcode"].isin(self.illumina_barcodes_table["illumina_indexed_barcode"])
        ].reset_index(drop=True)
       
        self.molecule_level_collate_df_il_itersect=molecule_level_collate_df_il_itersect
        print(
            "molecule_level_collate_df_il_itersect (after filter with Illumina):",
            "shape",molecule_level_collate_df_il_itersect.shape,
            "n_barcodes",molecule_level_collate_df_il_itersect["BCrev"].nunique(),
            "n_pacbio_tid",molecule_level_collate_df_il_itersect["pacbio_tid"].nunique()
        )
        self.selected_isoforms_il_intersect=self.selected_isoforms[
            self.selected_isoforms["pacbio_tid"].isin(self.molecule_level_collate_df_il_itersect["pacbio_tid"])
        ]


    def get_cell_level_isoform_molecule_count(self):
        cell_level_umi_count_table_fp=join(self.pacbio_illumina_combined_dir,"cell_level_umi_count_table.csv")
        if not os.path.isfile(cell_level_umi_count_table_fp):
            umi_series_list=list()
            display_name_list=list()
            for i,(index,row) in enumerate(self.selected_isoforms_il_intersect.iterrows()):
                print(f"{i+1}/{len(self.selected_isoforms)}")
                newtid=row["pacbio_tid"]
                c_df=self.molecule_level_collate_df_il_itersect[
                    self.molecule_level_collate_df_il_itersect["pacbio_tid"]==newtid
                ]
                umi_series_list.append(c_df.groupby("pacbio_indexed_barcode")["UMI"].nunique())
                display_name_list.append(row["display_name"])
            umi_count_table=pd.concat(umi_series_list,axis=1)
            umi_count_table.columns=display_name_list
            umi_count_table.to_csv(cell_level_umi_count_table_fp)
        else:
            umi_count_table=pd.read_csv(cell_level_umi_count_table_fp,index_col=0)
        return umi_count_table

    def load_fisher_table(self):
        fisher_table_fp=os.path.join(self.pacbio_illumina_combined_dir,"fisher_table.csv")
        fisher_table=pd.read_csv(fisher_table_fp,index_col=0)
        fisher_table["gene_symbol"]=fisher_table["newtid"].map(lambda x: x.split('.')[-2])
        fisher_table["pacbio_tid"]=fisher_table["newtid"].map(lambda x:'.'.join(x.split('.')[:-2]))
        
        return fisher_table

class SUPPA2EventCollection(object):
    def __init__(
        self,
        suppa2_event_fp,
        sample_name_map,
        class_df_all,
        merged_QTISO_mapping,
        sample_names=None,
        test_mode=False
    ):
        self.test_mode=test_mode
        self.suppa2_event_fp=suppa2_event_fp
        self.sample_name_map=sample_name_map
        self.class_df_all=class_df_all
        self.merged_QTISO_mapping=merged_QTISO_mapping
        if sample_names is None:
            self.sample_names=list(self.sample_name_map.values())
        else:
            self.sample_names=sample_names
        self.parse_suppa2_events()
    
    def parse_suppa2_events(self):
        class_df_all=self.class_df_all.set_index("pbid")
        print("="*SEP_SIZE+"Parsing SUPPA2 events"+"="*SEP_SIZE)
        suppa2_events_by_id={
            "A3":dict(),
            "A5":dict(),
            "AF":dict(),
            "AL":dict(),
            "MX":dict(),
            "RI":dict(),
            "SE":dict()
        }
        for event_type in suppa2_events_by_id.keys():
            if event_type in self.suppa2_event_fp["suppa2_local_event_fp"]:
                print("Parsing Event Type: {}".format(event_type))
                current_event_fp=self.suppa2_event_fp["suppa2_local_event_fp"][event_type]
                with open(current_event_fp) as f:
                    f.readline()
                    for line in tqdm(f,total=line_count(current_event_fp)):
                        line=line.rstrip()
                        suppa2event=SUPPA2Event.from_line(line)
                        if self.test_mode and suppa2event.chrom!=TEST_CHROM:
                            continue
                        assert suppa2event.id not in suppa2_events_by_id

                        supporting_query_isoform_id_pos=set()
                        for tid in suppa2event.positive_transcripts:
                            qtiso=self.merged_QTISO_mapping[tid]
                            supporting_query_isoform_id_pos.update(qtiso)

                        supporting_query_isoform_id_total=set()
                        associated_gene_id_set=set()
                        for tid in suppa2event.total_transcripts:
                            qtiso=self.merged_QTISO_mapping[tid]
                            supporting_query_isoform_id_total.update(qtiso)
                            for iso in qtiso:
                                if iso in class_df_all.index:
                                    gene = class_df_all.loc[iso,"associated_gene"]
                                    associated_gene_id_set.add((iso,gene))

                        supporting_query_isoform_id_neg=supporting_query_isoform_id_total.difference(supporting_query_isoform_id_pos)

                        suppa2event.update_associated_data(
                            "supporting_query_isoform_id_pos",
                            supporting_query_isoform_id_pos,
                            mode="assign"
                        )

                        suppa2event.update_associated_data(
                            "supporting_query_isoform_id_neg",
                            supporting_query_isoform_id_neg,
                            mode="assign"
                        )

                        suppa2event.update_associated_data(
                            "supporting_query_isoform_id_total",
                            supporting_query_isoform_id_total,
                            mode="assign"
                        )

                        suppa2event.update_associated_data(
                            "associated_gene_id",
                            associated_gene_id_set,
                            mode="assign"
                        )

                        suppa2_events_by_id[event_type][suppa2event.id]=suppa2event

        self.suppa2_events_by_id=suppa2_events_by_id
    
    def reset_differential_analysis(self):
        for feature_type, feature_dict in self.suppa2_events_by_id.items():
            for fid,feature in feature_dict.items():
                feature.reset_associated_data("differential_analysis")

    def _get_feature_collection(self,feature_name):
        if feature_name in ("A3","A5","AF","AL","MX","RI","SE"):
            feature_collection=self.suppa2_events_by_id[feature_name]
        else:
            assert False
        return feature_collection
    
    def print_shared_features(
        self,
        print_limit=5,
        seed=123,
        print_example=True
    ):
        np.random.seed(seed)
        
        for feature_name in (("A3","A5","AF","AL","MX","RI","SE")):
            iprint(f"*feature_name: {feature_name}")
            feature_collection=self._get_feature_collection(feature_name)
            
            iprint(f"**ALL: # events {len(feature_collection)}",ident=1)
            
            feature_collection_keys=list(feature_collection.keys())
            if print_example:
                for feature_id in np.random.choice(feature_collection_keys,min(print_limit,len(feature_collection_keys))):
                    iprint("->",feature_id,ident=1)
                    iprint("supporting_query_isoform_id pos",sorted(feature_collection[feature_id].associated_data.get("supporting_query_isoform_id_pos")),ident=2)
                    iprint("supporting_query_isoform_id neg",sorted(feature_collection[feature_id].associated_data.get("supporting_query_isoform_id_neg")),ident=2)

            feature_ids_differential_analysis= [k for k in feature_collection_keys if "differential_analysis" in feature_collection[k].associated_data]
            iprint(f"**Differential Analysis: # features {len(feature_ids_differential_analysis)}",ident=1)
            if print_example:
                for feature_id in np.random.choice(feature_ids_differential_analysis,min(print_limit,len(feature_ids_differential_analysis))):
                    diff_analysis=feature_collection[feature_id].associated_data.get("differential_analysis")
                    if len(diff_analysis["supporting_illumina_cell_barcodes_pos"])>0 and len(diff_analysis["supporting_illumina_cell_barcodes_neg"])>0:
                        iprint("->",feature_id,ident=1)
                        iprint("# supporting_illumina_cell_barcodes pos",len(diff_analysis["supporting_illumina_cell_barcodes_pos"]),ident=2)
                        iprint("# supporting_illumina_cell_barcodes neg",len(diff_analysis["supporting_illumina_cell_barcodes_neg"]),ident=2)

                        iprint("ctable",repr(diff_analysis["fisher_exact"]["ctable"]),ident=3)
                        iprint("oddsr",diff_analysis["fisher_exact"]["oddsr"],ident=3)
                        iprint("pvalue",diff_analysis["fisher_exact"]["pvalue"],ident=3)
                        iprint("cell_type_stats",diff_analysis["cell_type_stats"],ident=3)
                        iprint()

if __name__=="__main__":
    ils=IlluminaSampleInfo()
