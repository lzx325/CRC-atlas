from .base_classes import Exon,Junction,TSS,TTS,SpliceSite,TranscriptomicFeature
from .ref_classes import RefTranscriptIsoform
from utils.isoform_classes import genePredRecord
from typing import Dict, List, DefaultDict, Set
from collections import defaultdict
import pandas as pd

def _update_associated_data(associated_data, key: str, value, mode):
    if mode=="append":
        associated_data[key].add(value)
    elif mode=="assign":
        associated_data[key]=value
    elif mode=="update":
        associated_data[key].update(value)

def _reset_associated_data(associated_data ,key: str):
    if key in associated_data:
        del associated_data[key]
class QueryTranscriptIsoform(RefTranscriptIsoform):
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
        cls_record: pd.Series = None
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
        self.cls_record=cls_record
        self.associated_data: DefaultDict[str,Set[str]] =defaultdict(set)
    @staticmethod
    def from_genePred_record(record:genePredRecord, cls_record: pd.Series =None):
        if cls_record["coding"]=="coding":
            CDS_genomic_start=int(cls_record["CDS_genomic_start"])
            CDS_genomic_end_1based=int(cls_record["CDS_genomic_end_1based"])
        else:
            CDS_genomic_start=record.cdsStart
            CDS_genomic_end_1based=record.cdsEnd
        return QueryTranscriptIsoform(
            id=record.id,
            chrom=record.chrom,
            strand=record.strand,
            genomic_start=record.txStart,
            genomic_end_1based=record.txEnd,
            CDS_genomic_start=CDS_genomic_start,
            CDS_genomic_end_1based=CDS_genomic_end_1based,
            txStruct=record,
            genes=[record.gene],
            cls_record=cls_record
        )

    def update_associated_data(self,key,value,mode):
        _update_associated_data(self.associated_data,key,value,mode)

    def reset_associated_data(self,key):
        _reset_associated_data(self.associated_data,key)

