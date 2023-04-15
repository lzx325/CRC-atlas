from ast import Str

from utils.transcriptomic_features.ref_classes import RefTranscriptIsoform
from .base_classes import TranscriptomicFeature,Exon,Junction,TSS,TTS,SpliceSite,TranscriptStructure
from typing import DefaultDict, Dict, List, Set
from collections import defaultdict
import pandas as pd
from utils.isoform_classes import genePredRecord

def _update_associated_data(associated_data, key, value, mode):
    if mode=="append":
        associated_data[key].add(value)
    elif mode=="assign":
        associated_data[key]=value
    elif mode=="update":
        associated_data[key].update(value)

def _reset_associated_data(associated_data: DefaultDict[Str,Set[Str]],key: Str):
    if key in associated_data:
        del associated_data[key]

class SharedExon(Exon):
    def __init__(
        self,
        chrom: str,
        genomic_start: int,
        genomic_end_1based: int,
        strand: str,
        id: str = None,
    ):
        super().__init__(
            chrom=chrom,
            genomic_start=genomic_start,
            genomic_end_1based=genomic_end_1based,
            strand=strand,
            id=id
        )
        self.associated_data: DefaultDict[Str,Set[Str]] =defaultdict(set)
    @staticmethod
    def from_exon(exon:"Exon"):
        return SharedExon(
            chrom=exon.chrom,
            genomic_start=exon.genomic_start,
            genomic_end_1based=exon.genomic_end_1based,
            strand=exon.strand
        )
    
    def update_associated_data(self,key,value,mode):
        _update_associated_data(self.associated_data,key,value,mode)

    def reset_associated_data(self,key):
        _reset_associated_data(self.associated_data,key)

    def __eq__(self, rhs: "SharedExon"):
        assert isinstance(rhs,SharedExon)
        return (self.chrom,self.genomic_start,self.genomic_end_1based,self.strand)==(rhs.chrom,rhs.genomic_start,rhs.genomic_end_1based,rhs.strand)

    def __hash__(self):
        return hash((self.chrom,self.genomic_start,self.genomic_end_1based,self.strand))

class SharedJunction(Junction):
    def __init__(
        self,
        chrom: str,
        genomic_start: int,
        genomic_end_1based: int,
        strand: str,
        id: str = None,
    ):
        super().__init__(
            chrom=chrom,
            genomic_start=genomic_start,
            genomic_end_1based=genomic_end_1based,
            strand=strand,
            id=id
        )
        self.associated_data: DefaultDict[Str,Set[Str]] =defaultdict(set)
    @staticmethod
    def from_junction(junction:"Junction"):
        return SharedJunction(
            chrom=junction.chrom,
            genomic_start=junction.genomic_start,
            genomic_end_1based=junction.genomic_end_1based,
            strand=junction.strand
        )

    def __eq__(self, rhs: "SharedJunction"):
        assert isinstance(rhs,SharedJunction)
        return (self.chrom,self.genomic_start,self.genomic_end_1based,self.strand)==(rhs.chrom,rhs.genomic_start,rhs.genomic_end_1based,rhs.strand)

    def __hash__(self):
        return hash((self.chrom,self.genomic_start,self.genomic_end_1based,self.strand))

    def update_associated_data(self,key,value,mode):
        _update_associated_data(self.associated_data,key,value,mode)

    def reset_associated_data(self,key):
        _reset_associated_data(self.associated_data,key)

class SharedTSS(TSS):
    def __init__(
        self,
        chrom: str,
        genomic_center: int,
        strand: str,
        genomic_window_start: int = None,
        genomic_window_end_1based: int = None,
        id: str = None,

    ):
        super().__init__(
            chrom=chrom,
            genomic_center=genomic_center,
            strand=strand,
            genomic_window_start=genomic_window_start,
            genomic_window_end_1based=genomic_window_end_1based,
            id=id
        )
        self.associated_data: DefaultDict[Str,Set[Str]] =defaultdict(set)

    def __eq__(self,rhs:"SharedTSS"):
        assert isinstance(rhs,SharedTSS)
        return (self.chrom, self.genomic_center, self.strand, self.genomic_window_start, self.genomic_window_end_1based)\
                == (rhs.chrom, rhs.genomic_center, rhs.strand, rhs.genomic_window_start, rhs.genomic_window_end_1based)

    
    def __hash__(self):
        return hash((
            self.chrom,
            self.genomic_center,
            self.strand,
            self.genomic_window_start,
            self.genomic_window_end_1based,
        ))
    

    @staticmethod
    def from_TSS(tss: "TSS"):
        return SharedTSS(
            chrom=tss.chrom,
            genomic_center=tss.genomic_center,
            strand=tss.strand,
            genomic_window_start=tss.genomic_window_start,
            genomic_window_end_1based=tss.genomic_window_end_1based
        )
    
    def within_region(self,tss: "TSS"):
        return (self.chrom,self.strand)==(tss.chrom,tss.strand) \
            and self.genomic_window_start<=tss.genomic_center<self.genomic_window_end_1based

    def update_associated_data(self,key,value,mode):
        _update_associated_data(self.associated_data,key,value,mode)
    
    def reset_associated_data(self,key):
        _reset_associated_data(self.associated_data,key)

class SharedTTS(TTS):
    def __init__(
        self,
        chrom: str,
        genomic_center: int,
        strand: str,
        genomic_window_start: int = None,
        genomic_window_end_1based: int = None,
        id: str = None,

    ):
        super().__init__(
            chrom=chrom,
            genomic_center=genomic_center,
            strand=strand,
            genomic_window_start=genomic_window_start,
            genomic_window_end_1based=genomic_window_end_1based,
            id=id
        )
        self.associated_data: DefaultDict[Str,Set[Str]] =defaultdict(set)
    @staticmethod
    def from_TTS(tts: "TTS"):
        return SharedTTS(
            chrom=tts.chrom,
            genomic_center=tts.genomic_center,
            strand=tts.strand,
            genomic_window_start=tts.genomic_window_start,
            genomic_window_end_1based=tts.genomic_window_end_1based
        )

    def __eq__(self,rhs:"SharedTTS"):
        assert isinstance(rhs,SharedTTS)
        return (self.chrom,self.genomic_center,self.strand,self.genomic_window_start,self.genomic_window_end_1based)==(rhs.chrom,rhs.genomic_center,rhs.strand,rhs.genomic_window_start,rhs.genomic_window_end_1based)

    
    def __hash__(self):
        return hash((self.chrom,self.genomic_center,self.strand,self.genomic_window_start,self.genomic_window_end_1based))

    def update_associated_data(self,key,value,mode):
        _update_associated_data(self.associated_data,key,value,mode)
    
    def reset_associated_data(self,key):
        _reset_associated_data(self.associated_data,key)
        
class SharedSpliceSite(SpliceSite):
    def __init__(
        self,
        ss_type: str,
        chrom: str,
        genomic_loc: str,
        strand: str,
        id: str = None,

    ):
        super().__init__(
            ss_type=ss_type,
            chrom=chrom,
            genomic_loc=genomic_loc,
            strand=strand,
            id=id
        )
        self.associated_data: DefaultDict[Str,Set[Str]] =defaultdict(set)
    @staticmethod
    def from_splicesite(ss:"SpliceSite"):
        return SharedSpliceSite(
            ss_type=ss.ss_type,
            chrom=ss.chrom,
            genomic_loc=ss.genomic_loc,
            strand=ss.strand,
        )
    
    def __eq__(self,rhs:"SharedSpliceSite"):
        assert isinstance(rhs,SharedSpliceSite)
        return (self.ss_type,self.chrom,self.genomic_loc,self.strand)==(rhs.ss_type,rhs.chrom,rhs.genomic_loc,rhs.strand)
    
    def __hash__(self):
        return hash((self.ss_type,self.chrom,self.genomic_loc,self.strand))
    
    def update_associated_data(self,key,value,mode):
        _update_associated_data(self.associated_data,key,value,mode)
    
    def reset_associated_data(self,key):
        _reset_associated_data(self.associated_data,key)

class SharedTranscriptStructure(TranscriptStructure):
    def __init__(
        self,
        chrom: str,
        strand: str,
        junction_intervals: tuple,
        id: str = None,
        unique_by_id: bool = False,

    ):
        super().__init__(
            chrom=chrom,
            strand=strand,
            junction_intervals=junction_intervals,
            id = id
        )
        self.unique_by_id = unique_by_id
        self.associated_data: DefaultDict[Str,Set[Str]] =defaultdict(set)
    @staticmethod
    def from_transcript_structure(ijs: "TranscriptStructure", unique_by_id=False):
        return SharedTranscriptStructure(
            chrom=ijs.chrom,
            strand=ijs.strand,
            junction_intervals=ijs.junction_intervals,
            id=ijs.id if unique_by_id else None,
            unique_by_id=unique_by_id
        )
    
    def __eq__(self,rhs:"SharedTranscriptStructure"):
        assert isinstance(rhs,SharedTranscriptStructure)
        if not self.unique_by_id:
            return (self.unique_by_id, self.chrom,self.strand, self.junction_intervals)==(rhs.unique_by_id, rhs.chrom,rhs.strand, self.junction_intervals)
        else:
            return (self.unique_by_id,self.id) == (rhs.unique_by_id, rhs.id)
    
    def __hash__(self):
        if not self.unique_by_id:
            return hash((self.unique_by_id,self.chrom,self.strand, self.junction_intervals))
        else:
            return hash((self.unique_by_id,self.id))
    
    def update_associated_data(self,key,value,mode):
        _update_associated_data(self.associated_data,key,value,mode)
    
    def reset_associated_data(self,key):
        _reset_associated_data(self.associated_data,key)

class SharedTranscriptIsoform(RefTranscriptIsoform):
    def __init__(
        self,
        id: str,
        chrom: str,
        genomic_start: int, 
        genomic_end_1based: int,
        strand: str ,
        txStruct,
        CDS_genomic_start: int = -1,
        CDS_genomic_end_1based: int = -1,
        genes = None,
    ):
        super().__init__(
            id=id,
            chrom=chrom,
            strand=strand,
            genomic_start=genomic_start,
            genomic_end_1based=genomic_end_1based,
            CDS_genomic_start=CDS_genomic_start,
            CDS_genomic_end_1based=CDS_genomic_end_1based,
            txStruct=txStruct,
            genes=genes
        )
        self.associated_data: DefaultDict[str,Set[str]] =defaultdict(set)

    def from_genePred_record(record:genePredRecord, cls_record : pd.Series = None):
        if cls_record is not None and cls_record["coding"]=="coding":
            CDS_genomic_start=int(cls_record["CDS_genomic_start"])
            CDS_genomic_end_1based=int(cls_record["CDS_genomic_end_1based"])
        else:
            CDS_genomic_start=record.cdsStart
            CDS_genomic_end_1based=record.cdsEnd

        return SharedTranscriptIsoform(
            id=record.id,
            chrom=record.chrom,
            strand=record.strand,
            genomic_start=record.txStart,
            genomic_end_1based=record.txEnd,
            CDS_genomic_start=CDS_genomic_start,
            CDS_genomic_end_1based=CDS_genomic_end_1based,
            txStruct=record,
            genes=[record.gene]
        )

    def update_associated_data(self,key,value,mode):
        _update_associated_data(self.associated_data,key,value,mode)
    
    def reset_associated_data(self,key):
        _reset_associated_data(self.associated_data,key)
class SUPPA2Event(TranscriptomicFeature):
    exported_fields=[
        "id",
        "chrom",
        "strand"
    ]

    def __init__(
        self,
        chrom: str,
        genomic_start: str,
        genomic_end_1based: str,
        strand: str,
        id: str,
        event_type: str,
        positive_transcripts: List[str],
        total_transcripts: List[str]
    ):
        super().__init__()
        self.chrom=chrom
        self.genomic_start=genomic_start
        self.genomic_end_1based=genomic_end_1based
        self.strand=strand
        self.id=id
        self.event_type=event_type
        self.positive_transcripts=set(positive_transcripts)
        self.total_transcripts=set(total_transcripts)
        assert all(t in self.total_transcripts for t in self.positive_transcripts)
        self.associated_data: DefaultDict[Str,Set[Str]] =defaultdict(set)
    @staticmethod
    def from_line(line) -> "SUPPA2Event":
        line=line.rstrip()
        fields=line.split("\t")
        seqname=fields[0]
        gene_id=fields[1]
        first_part,second_part=fields[2].split(";")
        id=second_part
        second_part_fields=second_part.split(":")
        event_type=second_part_fields[0]
        strand=second_part_fields[-1]
        intervals=[field.split("-") for field in second_part_fields[2:-1] if '-' in field]
        start_min=min(int(interval[0])-1 for interval in intervals) # make 0-based
        end_1based_max=max(int(interval[1]) for interval in intervals)
        positive_transcripts=fields[3].split(",")
        total_transcripts=fields[4].split(",")
        return SUPPA2Event(
            chrom=seqname,
            genomic_start=start_min,
            genomic_end_1based=end_1based_max,
            strand=strand,
            id=id,
            event_type=event_type,
            positive_transcripts=positive_transcripts,
            total_transcripts=total_transcripts
        )
    def __str__(self):
        return str((self.id,self.event_type,self.chrom,self.genomic_start,self.genomic_end_1based,self.strand))
    def __repr__(self):
        return self.__str__()

    def update_associated_data(self,key,value,mode):
        _update_associated_data(self.associated_data,key,value,mode)
    
    def reset_associated_data(self,key):
        _reset_associated_data(self.associated_data,key)