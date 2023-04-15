from copy import copy
from collections import defaultdict
from bx.intervals import Interval, IntervalTree
from tqdm import tqdm
import pandas as pd
import re

from utils.isoform_classes import genePredReader, genePredRecord
from utils.file_utils import line_count
from .shared_classes import \
    SharedTTS, SharedTSS, SharedExon, SharedJunction, SharedSpliceSite, SharedTranscriptStructure
from .ref_classes import RefTranscriptIsoform
from utils.utils import canonical_chr_regex

def update_shared_entity_type_id(shared_entity,id_map):
    shared_entity=copy(shared_entity)
    id_split=shared_entity.id.split("__")
    id_split[0]=id_map[id_split[0]]
    new_id="__".join(id_split)
    shared_entity.id=new_id
    return shared_entity

def update_ref_iso_collections(
    ref_iso_by_id,ref_iso_itree_by_chr,ref_iso_itree_by_chr2,
    ref_ti,
    type_id_map
):
    ref_ti=copy(ref_ti)
    ref_ti.id="RTISO__"+ref_ti.id
    assert ref_ti.id not in ref_iso_by_id
    ref_iso_by_id[ref_ti.id]=ref_ti
    ref_iso_itree_by_chr[(ref_ti.chrom,ref_ti.strand)].insert(ref_ti.genomic_start,ref_ti.genomic_end_1based,ref_ti.id)
    ref_iso_itree_by_chr2[(ref_ti.chrom,ref_ti.strand)].insert(ref_ti.genomic_start,ref_ti.genomic_end_1based,ref_ti.id)

def update_ref_ts_collections(
    shared_ts_by_id,shared_ts_to_id,
    ref_ti,
    type_id_map
):
    shared_ts = SharedTranscriptStructure.from_transcript_structure(ref_ti.transcript_structure,unique_by_id=True)
    shared_ts = update_shared_entity_type_id(shared_ts,type_id_map)
    
    if shared_ts not in shared_ts_to_id:
        assert shared_ts.id not in shared_ts_by_id
        shared_ts_by_id[shared_ts.id]=shared_ts
        shared_ts_to_id[shared_ts]=shared_ts.id
        existing_id=shared_ts.id
    else:
        existing_id = shared_ts_to_id[shared_ts]
    shared_ts = shared_ts_by_id[existing_id]
    shared_ts.update_associated_data("associated_ref_id",ref_ti.transcript_structure.id,mode="append")

def update_ref_exons_collections(
    shared_exon_by_id,shared_exon_to_id,exon_itree_by_chr,exon_sets_by_gene,
    ref_ti,
    type_id_map
):
    for exon in ref_ti.exons.values():
        shared_exon=SharedExon.from_exon(exon)
        shared_exon = update_shared_entity_type_id(shared_exon,type_id_map)
        if shared_exon not in shared_exon_to_id:
            assert shared_exon.id not in shared_exon_by_id
            shared_exon_by_id[shared_exon.id]=shared_exon
            shared_exon_to_id[shared_exon]=shared_exon.id
            
            exon_itree_by_chr[(ref_ti.chrom,ref_ti.strand)].insert(shared_exon.genomic_start,shared_exon.genomic_end_1based,shared_exon)
            existing_id = shared_exon.id
        else:
            existing_id = shared_exon_to_id[shared_exon]
        shared_exon=shared_exon_by_id[existing_id]
        shared_exon.update_associated_data("associated_ref_id",exon.id,mode="append")
        for gene_id in ref_ti.genes:
            shared_exon.update_associated_data("associated_gene_id",(exon.associated_transcript_id,gene_id),mode="append")
            exon_sets_by_gene[gene_id].add(shared_exon)
        
def update_ref_junctions_collections(
    shared_junction_by_id,shared_junction_to_id,junction_itree_by_chr,junction_sets_by_gene,
    ref_ti,
    type_id_map
):
    for junc in ref_ti.junctions.values():
        shared_junc=SharedJunction.from_junction(junc)
        shared_junc = update_shared_entity_type_id(shared_junc,type_id_map)
        if shared_junc not in shared_junction_to_id:
            assert junc.id not in shared_junction_by_id
            shared_junction_by_id[shared_junc.id]=shared_junc
            shared_junction_to_id[shared_junc]=shared_junc.id
            
            junction_itree_by_chr[(ref_ti.chrom,ref_ti.strand)].insert(shared_junc.genomic_start,shared_junc.genomic_end_1based,shared_junc)
            existing_id = shared_junc.id
        else:
            existing_id = shared_junction_to_id[shared_junc]
        shared_junction=shared_junction_by_id[existing_id]
        shared_junction.update_associated_data("associated_ref_id",junc.id,mode="append")

        for gene_id in ref_ti.genes:
            shared_junction.update_associated_data("associated_gene_id",(junc.associated_transcript_id,gene_id),mode="append")
            junction_sets_by_gene[gene_id].add(shared_junction)
        
def update_donor_collections(
    shared_donors_by_id, shared_donors_to_id, donor_itree_by_chr, donor_sets_by_gene,
    ref_ti,
    type_id_map
):
    donor_sites_list=[junc.splice_donor for junc in ref_ti.junctions.values()]
    
    for donor in donor_sites_list:
        shared_donor=SharedSpliceSite.from_splicesite(donor)
        shared_donor = update_shared_entity_type_id(shared_donor,type_id_map)
        
        if shared_donor not in shared_donors_to_id:
            assert shared_donor.id not in shared_donors_by_id
            shared_donors_by_id[shared_donor.id]=shared_donor
            shared_donors_to_id[shared_donor]=shared_donor.id
            
            donor_itree_by_chr[(ref_ti.chrom,ref_ti.strand)].insert(shared_donor.genomic_loc,shared_donor.genomic_loc+1,shared_donor)
            
            existing_id = shared_donor.id
        else:
            existing_id = shared_donors_to_id[shared_donor]
        shared_donor=shared_donors_by_id[existing_id]
        shared_donor.update_associated_data("associated_ref_id",donor.id,mode="append")
        for gene_id in ref_ti.genes:
            shared_donor.update_associated_data("associated_gene_id",(ref_ti.id,gene_id),mode="append")
            donor_sets_by_gene[gene_id].add(shared_donor)
def update_acceptor_collections(
    shared_acceptors_by_id, shared_acceptors_to_id, acceptor_itree_by_chr, acceptor_sets_by_gene,
    ref_ti,
    type_id_map
):
    
    acceptor_sites_list=[junc.splice_acceptor for junc in ref_ti.junctions.values()]
        
    for acceptor in acceptor_sites_list:
        shared_acceptor = SharedSpliceSite.from_splicesite(acceptor)
        shared_acceptor = update_shared_entity_type_id(shared_acceptor,type_id_map)
        
        if shared_acceptor not in shared_acceptors_to_id:
            assert shared_acceptor.id not in shared_acceptors_by_id
            shared_acceptors_by_id[shared_acceptor.id]=shared_acceptor
            shared_acceptors_to_id[shared_acceptor]=shared_acceptor.id
            
            acceptor_itree_by_chr[(ref_ti.chrom,ref_ti.strand)].insert(shared_acceptor.genomic_loc,shared_acceptor.genomic_loc+1,shared_acceptor)
            existing_id = shared_acceptor.id
        else:
            existing_id = shared_acceptors_to_id[shared_acceptor]
            
        shared_acceptor=shared_acceptors_by_id[existing_id]
        shared_acceptor.update_associated_data("associated_ref_id",acceptor.id,mode="append")
        for gene_id in ref_ti.genes:
            shared_acceptor.update_associated_data("associated_gene_id",(ref_ti.id,gene_id),mode="append")
            acceptor_sets_by_gene[gene_id].add(shared_acceptor)

def update_tss_collections(
    shared_tss_by_id,shared_tss_to_id,tss_itree_by_chr,tss_sets_by_gene,
    bed_row,
    type_id_map,
    reference_isoforms
):
    shared_tss=SharedTSS(
        chrom=bed_row["chrom"],
        genomic_center=(bed_row["start"]+bed_row["end"])//2,
        strand=bed_row["strand"],
        genomic_window_start=bed_row["start"],
        genomic_window_end_1based=bed_row["end"]
    )
    shared_tss = update_shared_entity_type_id(shared_tss,type_id_map)
    if shared_tss not in shared_tss_to_id:
        assert shared_tss.id not in shared_tss_by_id
        shared_tss_by_id[shared_tss.id]=shared_tss
        shared_tss_to_id[shared_tss]=shared_tss.id
        tss_itree_by_chr[(shared_tss.chrom,shared_tss.strand)].insert(
            shared_tss.genomic_window_start,
            shared_tss.genomic_window_end_1based,
            shared_tss
        )
        
        existing_id = shared_tss.id
    else:
        existing_id = shared_tss_to_id[shared_tss]

    shared_tss = shared_tss_by_id[existing_id]
    for associated_ref_id in bed_row["name"].split(","):
        shared_tss.update_associated_data("associated_ref_id",associated_ref_id,mode="append")
        if associated_ref_id.startswith("GV37RS:"):
            associated_ref_id=re.sub("GV37RS:TSS__","RTISO__",associated_ref_id)
            if associated_ref_id in reference_isoforms:
                for gene_id in reference_isoforms[associated_ref_id].genes:
                    shared_tss.update_associated_data("associated_gene_id",(associated_ref_id,gene_id),mode="append")
                    tss_sets_by_gene[gene_id].add(shared_tss)

def update_tts_collections(
    shared_tts_by_id,shared_tts_to_id,tts_itree_by_chr,tts_sets_by_gene,
    bed_row,
    type_id_map,
    reference_isoforms
):
    shared_tts=SharedTTS(
        chrom=bed_row["chrom"],
        genomic_center=(bed_row["start"]+bed_row["end"])//2,
        strand=bed_row["strand"],
        genomic_window_start=bed_row["start"],
        genomic_window_end_1based=bed_row["end"]
    )
    shared_tts = update_shared_entity_type_id(shared_tts,type_id_map)
    if shared_tts not in shared_tts_to_id:
        assert shared_tts.id not in shared_tts_by_id
        shared_tts_by_id[shared_tts.id]=shared_tts
        shared_tts_to_id[shared_tts]=shared_tts.id
        tts_itree_by_chr[(shared_tts.chrom,shared_tts.strand)].insert(
            shared_tts.genomic_window_start,
            shared_tts.genomic_window_end_1based,
            shared_tts
        )
        
        existing_id = shared_tts.id
    else:
        existing_id = shared_tts_to_id[shared_tts]
    shared_tts = shared_tts_by_id[existing_id]
    for associated_ref_id in bed_row["name"].split(","):
        shared_tts.update_associated_data("associated_ref_id",associated_ref_id,mode="append")
        if associated_ref_id.startswith("GV37RS:"):
            associated_ref_id=re.sub("GV37RS:TTS__","RTISO__",associated_ref_id)
            if associated_ref_id in reference_isoforms:
                for gene_id in reference_isoforms[associated_ref_id].genes:
                    shared_tts.update_associated_data("associated_gene_id",(associated_ref_id,gene_id),mode="append")
                    tts_sets_by_gene[gene_id].add(shared_tts)

def init_reference_isoform_collections():
    ref_iso_by_id = dict()
    ref_iso_itree_by_chr = defaultdict(lambda: IntervalTree())
    ref_iso_1exon_itree_by_chr = defaultdict(lambda: IntervalTree())
    ref_iso_multiexon_itree_by_chr = defaultdict(lambda: IntervalTree())
    return ref_iso_by_id, ref_iso_itree_by_chr, ref_iso_1exon_itree_by_chr, ref_iso_multiexon_itree_by_chr

def init_transcriptomic_feature_collections():
    transcriptomic_feature_by_id = {'ts':dict(), 'exons':dict(), 'junctions': dict(), 'donors': dict(), 'acceptors': dict()}
    transcriptomic_feature_to_id = {'ts':dict(), 'exons':dict(), 'junctions': dict(), 'donors': dict(), 'acceptors': dict()}
    transcriptomic_feature_itree_by_chr = {
        "exons": defaultdict(lambda: IntervalTree()),
        "junctions": defaultdict(lambda: IntervalTree()),
        "donors": defaultdict(lambda: IntervalTree()),
        "acceptors": defaultdict(lambda: IntervalTree())
    }

    transcriptomic_feature_sets_by_gene = {
        "exons": defaultdict(lambda: set()),
        "junctions": defaultdict(lambda: set()),
        "donors": defaultdict(lambda: set()),
        "acceptors": defaultdict(lambda: set())
    }

    return transcriptomic_feature_by_id, \
        transcriptomic_feature_to_id, \
        transcriptomic_feature_itree_by_chr, \
        transcriptomic_feature_sets_by_gene

def init_TSS_TTS_feature_collections():
    TSS_TTS_feature_by_id = {"tss": dict(), "tts": dict()}
    TSS_TTS_feature_to_id = {"tss": dict(), "tts": dict()}

    TSS_TTS_feature_itree_by_chr = {
        "tss": defaultdict(lambda: IntervalTree()),
        "tts": defaultdict(lambda: IntervalTree())
    }

    TSS_TTS_feature_sets_by_gene = {
        "tss": defaultdict(lambda: set()),
        "tts": defaultdict(lambda: set())
    }

    return TSS_TTS_feature_by_id, \
        TSS_TTS_feature_to_id, \
        TSS_TTS_feature_itree_by_chr, \
        TSS_TTS_feature_sets_by_gene
        
def reference_genePred_parser(ref_genePred_fp,test_mode_chr=None, filter_length=True):
    # parse reference annotation

    ## Principles
    ### 1. ignore all miRNAs (< 200 bp)
    ### 2. separately store single exon and multi-exon references

    ## Parse reference isoforms (RISO)
    ref_iso_by_id, ref_iso_itree_by_chr, ref_iso_1exon_itree_by_chr, ref_iso_multiexon_itree_by_chr = \
        init_reference_isoform_collections()

    ## Parse reference internal junction structures, exons and junctions (RIJS, REXON, RJUNC, deduplicated)
    transcriptomic_feature_by_id, \
    transcriptomic_feature_to_id, \
    transcriptomic_feature_itree_by_chr, \
    transcriptomic_feature_sets_by_gene = init_transcriptomic_feature_collections()

    ## count number of dropped transcripts
    num_dropped_transcript_isoforms=0
    
    ## type_id_map
    type_id_map={
        "TS":"RTS",
        "EXON":"REXON",
        "JUNC":"RJUNC",
        "SD":"RSD",
        "SA":"RSA"
    }
    for i,r in tqdm(enumerate(genePredReader(ref_genePred_fp)),total=line_count(ref_genePred_fp)):
        if (filter_length and (r.length < 200 or not re.match(canonical_chr_regex(),r.chrom) ) ) or (test_mode_chr is not None and r.chrom!=test_mode_chr):
            num_dropped_transcript_isoforms+=1
            continue
        ref_ti=RefTranscriptIsoform.from_genePred_record(r)
        
        if ref_ti.num_exons == 1:
            update_ref_iso_collections(
                ref_iso_by_id,ref_iso_itree_by_chr,ref_iso_1exon_itree_by_chr,
                ref_ti,
                type_id_map
            )
            
            update_ref_ts_collections(
                transcriptomic_feature_by_id["ts"],transcriptomic_feature_to_id["ts"],
                ref_ti,
                type_id_map
            )
            
            update_ref_exons_collections(
                transcriptomic_feature_by_id["exons"],transcriptomic_feature_to_id["exons"],transcriptomic_feature_itree_by_chr["exons"],transcriptomic_feature_sets_by_gene["exons"],
                ref_ti,
                type_id_map
            )
        else:
            update_ref_iso_collections(
                ref_iso_by_id,ref_iso_itree_by_chr,ref_iso_multiexon_itree_by_chr,
                ref_ti,
                type_id_map
            )
            
            update_ref_ts_collections(
                transcriptomic_feature_by_id["ts"],transcriptomic_feature_to_id["ts"],
                ref_ti,
                type_id_map
            )
            
            update_ref_exons_collections(
                transcriptomic_feature_by_id["exons"],transcriptomic_feature_to_id["exons"],transcriptomic_feature_itree_by_chr["exons"],transcriptomic_feature_sets_by_gene["exons"],
                ref_ti,
                type_id_map
            )
            
            update_ref_junctions_collections(
                transcriptomic_feature_by_id["junctions"],transcriptomic_feature_to_id["junctions"],transcriptomic_feature_itree_by_chr["junctions"],transcriptomic_feature_sets_by_gene["junctions"],
                ref_ti,
                type_id_map
            )
            
            update_donor_collections(
                transcriptomic_feature_by_id["donors"], transcriptomic_feature_to_id["donors"],transcriptomic_feature_itree_by_chr["donors"],transcriptomic_feature_sets_by_gene["donors"],
                ref_ti,
                type_id_map
            )
            
            update_acceptor_collections(
                transcriptomic_feature_by_id["acceptors"], transcriptomic_feature_to_id["acceptors"],transcriptomic_feature_itree_by_chr["acceptors"],transcriptomic_feature_sets_by_gene["acceptors"],
                ref_ti,
                type_id_map
            )
            
    # convert the content of junctions_by_chr to sorted list
    ref_iso_itree_by_chr = dict(ref_iso_itree_by_chr)
    ref_iso_1exon_itree_by_chr = dict(ref_iso_1exon_itree_by_chr)
    ref_iso_multiexon_itree_by_chr = dict(ref_iso_multiexon_itree_by_chr)
    transcriptomic_feature_itree_by_chr = {k: dict(v) for k,v in transcriptomic_feature_itree_by_chr.items()}
    transcriptomic_feature_sets_by_gene = {k: dict(v) for k,v in transcriptomic_feature_sets_by_gene.items()}

    ref_iso_searchable_data_struct={
        "ref_iso_itree_by_chr":ref_iso_itree_by_chr,
        "ref_iso_1exon_itree_by_chr":ref_iso_1exon_itree_by_chr,
        "ref_iso_multiexon_itree_by_chr":ref_iso_multiexon_itree_by_chr,
    }
    print("# dropped transcript isoforms: ",num_dropped_transcript_isoforms)

    return ref_iso_searchable_data_struct, \
            ref_iso_by_id, \
            transcriptomic_feature_by_id, transcriptomic_feature_to_id, \
            transcriptomic_feature_itree_by_chr, transcriptomic_feature_sets_by_gene

def reference_TSS_TTS_parser(
        ref_tss_bed_fp,ref_tts_bed_fp,
        reference_isoforms,
        test_mode_chr=None
    ):
    ## Parse reference TSS TTS
    TSS_TTS_feature_by_id, \
    TSS_TTS_feature_to_id, \
    TSS_TTS_feature_itree_by_chr, \
    TSS_TTS_by_gene = init_TSS_TTS_feature_collections()

    type_id_map={
        "TSS":"RTSS",
        "TTS":"RTTS"
    }

    tss_bed=pd.read_csv(ref_tss_bed_fp,sep="\t",header=None)
    tss_bed.columns=["chrom","start","end","name","score","strand","thick_start", "thick_end","color"]
    if test_mode_chr is not None:
        tss_bed=tss_bed[tss_bed["chrom"]==test_mode_chr]

    print("Parsing TSS...")
    for _, bed_row in tqdm(tss_bed.iterrows(),total=tss_bed.shape[0]):
        update_tss_collections(
            TSS_TTS_feature_by_id["tss"],TSS_TTS_feature_to_id["tss"],TSS_TTS_feature_itree_by_chr["tss"],TSS_TTS_by_gene["tss"],
            bed_row,
            type_id_map,
            reference_isoforms
        )

    tts_bed=pd.read_csv(ref_tts_bed_fp,sep="\t",header=None)
    tts_bed.columns=["chrom","start","end","name","score","strand","thick_start", "thick_end","color"]

    if test_mode_chr is not None:
        tts_bed=tts_bed[tts_bed["chrom"]==test_mode_chr]

    print("Parsing TTS...")
    for _, bed_row in tqdm(tts_bed.iterrows(),total=tts_bed.shape[0]):
        update_tts_collections(
            TSS_TTS_feature_by_id["tts"],TSS_TTS_feature_to_id["tts"],TSS_TTS_feature_itree_by_chr["tts"],TSS_TTS_by_gene["tts"],
            bed_row,
            type_id_map,
            reference_isoforms
        )

    TSS_TTS_feature_itree_by_chr = {k:dict(v) for k,v in TSS_TTS_feature_itree_by_chr.items()}
    TSS_TTS_by_gene = {k:dict(v) for k,v in TSS_TTS_by_gene.items()}

    return TSS_TTS_feature_by_id, TSS_TTS_feature_to_id, TSS_TTS_feature_itree_by_chr, TSS_TTS_by_gene