from logging import StreamHandler
from Bio.Seq import Seq
import pyfaidx
from typing import Dict, List
from collections import defaultdict
from typing import Dict, List
from utils.isoform_classes import genePredRecord



class TranscriptomicFeature(object):
    def as_dict(self):
        d = {
            k:getattr(self,k) for k in type(self).exported_fields
        }
        
        return d

    def update_fields(self,update_dict):
        for k,v in update_dict.items():
            assert k in type(self).exported_fields, k
            setattr(self,k,v)

class Exon(TranscriptomicFeature):
    exported_fields=[
        "id",
        "chrom",
        "genomic_start",
        "genomic_end",
        "strand",
    ]

    def __init__(
        self,
        chrom: str,
        genomic_start: int,
        genomic_end_1based: int,
        strand: str,
        id: str = None
    ):
        super().__init__()
        self.chrom=chrom
        self.genomic_start=genomic_start
        self.genomic_end_1based=genomic_end_1based
        assert strand in ("+","-")
        self.strand=strand
        if id!= None:
            self.id=id
        else:
            self.id=f"EXON__{self.chrom}:{self.genomic_start}-{self.genomic_end_1based}({self.strand})"

    def __str__(self):
        return str((self.id,self.chrom,self.genomic_start,self.genomic_end_1based,self.strand))
    
    def __repr__(self):
        return self.__str__()

    @property
    def exon_head(self):
        if self.strand=="+":
            return self.genomic_start
        elif self.strand=="-":
            return self.genomic_end_1based-1
    
    @property
    def exon_tail(self):
        if self.strand=="+":
            return self.genomic_end_1based-1
        if self.strand=="-":
            return self.genomic_start
    
    @property
    def length(self):
        return self.genomic_end_1based-self.genomic_start
    
    @property
    def genomic_interval(self):
        return (self.genomic_start,self.genomic_end_1based)

class Junction(TranscriptomicFeature):
    exported_fields=[
        "id",
        "chrom",
        "genomic_start",
        "genomic_end",
        "strand",
    ]
    def __init__(
        self,
        chrom: str,
        genomic_start: int,
        genomic_end_1based: int,
        strand: str,
        id: str = None
    ):
        super().__init__()
        self.chrom=chrom
        self.genomic_start=genomic_start
        self.genomic_end_1based=genomic_end_1based
        assert strand in ("+","-")
        self.strand=strand
        if id!= None:
            self.id=id
        else:
            self.id=f"JUNC__{self.chrom}:{self.genomic_start}-{self.genomic_end_1based}({self.strand})"
        self.get_splice_sites()

    def get_splice_sites(self):
        self.splice_donor=SpliceSite(
            ss_type="donor",
            chrom=self.chrom,
            genomic_loc=self.genomic_start if self.strand=="+" else self.genomic_end_1based-1,
            strand=self.strand
        )

        self.splice_acceptor=SpliceSite(
            ss_type="acceptor",
            chrom=self.chrom,
            genomic_loc=self.genomic_end_1based-1 if self.strand=="+" else self.genomic_start,
            strand=self.strand
        )
    @property
    def junction_head(self):
        if self.strand=="+":
            return self.genomic_start
        elif self.strand=="-":
            return self.genomic_end_1based-1
    
    @property
    def junction_tail(self):
        if self.strand=="+":
            return self.genomic_end_1based-1
        if self.strand=="-":
            return self.genomic_start
    
    @property
    def length(self):
        return self.genomic_end_1based-self.genomic_start
    
    @property
    def genomic_interval(self):
        return (self.genomic_start,self.genomic_end_1based)

    def get_splice_site_seq(self, genome_dict):
        """
        Return the donor-acceptor site (ex: GTAG) for the i-th junction
        :param i: 0-based junction index
        :param genome_dict: dict of chrom --> SeqRecord
        :return: splice site pattern, ex: "GTAG", "GCAG" etc
        """
        chrom_seq=genome_dict[self.chrom]
        if isinstance(chrom_seq,Seq):
            seq_d = genome_dict[self.chrom].seq[self.genomic_start:self.genomic_start+2]
            seq_a = genome_dict[self.chrom].seq[(self.genomic_end_1based-2):self.genomic_end_1based]
        elif isinstance(chrom_seq,pyfaidx.FastaRecord):
            seq_d = Seq(genome_dict[self.chrom][self.genomic_start:self.genomic_start+2].seq)
            seq_a = Seq(genome_dict[self.chrom][(self.genomic_end_1based-2):self.genomic_end_1based].seq)
        elif isinstance(chrom_seq,str):
            seq_d = Seq(genome_dict[self.chrom][self.genomic_start:self.genomic_start+2])
            seq_a = Seq(genome_dict[self.chrom][(self.genomic_end_1based-2):self.genomic_end_1based])

        if self.strand == '+':
            return (str(seq_d)+str(seq_a)).upper()
        else:
            return (str(seq_a.reverse_complement())+str(seq_d.reverse_complement())).upper()

    def __str__(self):
        return str((self.id,self.chrom,self.genomic_start,self.genomic_end_1based,self.strand))

    def __repr__(self):
        return self.__str__()

class TranscriptStructure(TranscriptomicFeature):
    def __init__(
        self,
        chrom: str,
        strand: str,
        junction_intervals: tuple,
        id: str = None
    ):
        super().__init__()
        assert strand in ("+","-")
        self.chrom=chrom
        self.strand=strand
        self.junction_intervals=junction_intervals
        self.junction_intervals_str=','.join("{}-{}".format(s,e) for s,e in self.junction_intervals)
        if id != None:
            self.id = id
        else:
            self.id = f"TS__{self.chrom}:{self.junction_intervals_str}({self.strand})"
    def __str__(self):
        return str((self.id,self.chrom,self.junction_intervals_str,self.strand))
    def __repr__(self):
        return self.__str__()
        
class SpliceSite(TranscriptomicFeature):
    exported_fields=[
        "id",
        "ss_type",
        "chrom",
        "genomic_loc",
        "strand"
    ]
    def __init__(
        self,
        ss_type: str,
        chrom: str,
        genomic_loc: str,
        strand: str,
        id: str = None
    ):
        super().__init__()
        assert ss_type in ("donor","acceptor")
        self.ss_type=ss_type
        self.chrom=chrom
        self.genomic_loc=genomic_loc
        self.strand=strand
        if id != None:
            self.id = id
        else:
            if self.ss_type=="donor":
                self.id = f"SD__{self.chrom}:{self.genomic_loc}({self.strand})"
            elif self.ss_type=="acceptor":
                self.id = f"SA__{self.chrom}:{self.genomic_loc}({self.strand})"
    def __str__(self):
        return str((self.id,self.chrom,self.genomic_loc,self.strand,self.ss_type))
    def __repr__(self):
        return self.__str__()

class TSS(TranscriptomicFeature):
    exported_fields=[
        "id",
        "chrom"
        "genomic_center",
        "strand"
        "genomic_window_start",
        "genomic_window_end_1based",
    ]

    def __init__(
        self,
        chrom: str,
        genomic_center: int,
        strand: str,
        genomic_window_start: int = None,
        genomic_window_end_1based: int = None,
        id: str = None
    ):
        super().__init__()
        self.chrom=chrom
        self.genomic_center=genomic_center
        assert strand in ("+","-")
        self.strand=strand
        self.genomic_window_start=genomic_window_start
        self.genomic_window_end_1based=genomic_window_end_1based

        if self.genomic_window_start!=None and self.genomic_window_end_1based!=None:
            assert self.genomic_window_start<=self.genomic_center<self.genomic_window_end_1based
        if id!=None:
            self.id=id
        else:
            if self.genomic_window_start!=None and self.genomic_window_end_1based!=None:
                self.id=f"TSS__{self.chrom}:{self.genomic_center}({self.strand})[{genomic_window_start}-{self.genomic_window_end_1based}]"
            else:
                self.id=f"TSS__{self.chrom}:{self.genomic_center}({self.strand})"

    def __str__(self):
        if self.genomic_window_start!=None and self.genomic_window_end_1based!=None:
            return str((self.id,self.chrom,self.genomic_center,self.strand,self.genomic_window_start,self.genomic_window_end_1based))
        else:
            return str((self.id,self.chrom,self.genomic_center,self.strand))

    def __repr__(self):
        return self.__str__()

class TTS(TranscriptomicFeature):
    exported_fields=[
        "id",
        "chrom"
        "genomic_center",
        "strand"
        "genomic_window_start",
        "genomic_window_end_1based",
    ]

    def __init__(
        self,
        chrom: str,
        genomic_center: int,
        strand: str,
        genomic_window_start: int = None,
        genomic_window_end_1based: int = None,
        id: str = None
    ):
        super().__init__()
        self.chrom=chrom
        self.genomic_center=genomic_center
        assert strand in ("+","-")
        self.strand=strand
        self.genomic_window_start=genomic_window_start
        self.genomic_window_end_1based=genomic_window_end_1based

        if self.genomic_window_start!=None and self.genomic_window_end_1based!=None:
            assert self.genomic_window_start<=self.genomic_center<self.genomic_window_end_1based
        if id!=None:
            self.id=id
        else:
            if self.genomic_window_start!=None and self.genomic_window_end_1based!=None:
                self.id=f"TTS__{self.chrom}:{self.genomic_center}({self.strand})[{genomic_window_start}-{self.genomic_window_end_1based}]"
            else:
                self.id=f"TTS__{self.chrom}:{self.genomic_center}({self.strand})"

    def __str__(self):
        if self.genomic_window_start!=None and self.genomic_window_end_1based!=None:
            return str((self.id,self.chrom,self.genomic_center,self.strand,self.genomic_window_start,self.genomic_window_end_1based))
        else:
            return str((self.id,self.chrom,self.genomic_center,self.strand))

    def __repr__(self):
        return self.__str__()
    
