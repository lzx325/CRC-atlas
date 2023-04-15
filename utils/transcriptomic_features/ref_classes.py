from .base_classes import Exon,Junction,TSS,TTS,SpliceSite,TranscriptomicFeature,TranscriptStructure
from utils.isoform_classes import genePredRecord
from typing import Dict, List
from Bio.Seq import Seq
import pyfaidx
from utils.utils import interval_overlap
class RefTranscriptExon(Exon):
    exported_fields=[
        "id",
        "associated_transcript_id",
        "exon_index",
        "chrom",
        "genomic_start",
        "genomic_end",
        "strand",
    ]
    def __init__(
        self,
        associated_transcript_id: str,
        exon_index: int,
        chrom: str,
        genomic_start: int,
        genomic_end_1based: int,
        strand: str,
        id = None,
    ):
        super().__init__(chrom=chrom,genomic_start=genomic_start,genomic_end_1based=genomic_end_1based,strand=strand)
        self.associated_transcript_id=associated_transcript_id
        self.exon_index=exon_index
        if id==None:
            self.id="EXON__{}__{}".format(self.associated_transcript_id,"exon_"+str(self.exon_index))
    def __eq__(self,rhs:"RefTranscriptExon"):
        assert isinstance(rhs,RefTranscriptExon)
        return self.associated_transcript_id==rhs.associated_transcript_id and self.chrom==rhs.chrom and self.genomic_start==rhs.genomic_start and self.genomic_end_1based==rhs.genomic_end_1based and self.strand==rhs.strand

class RefTranscriptJunction(Junction):
    exported_fields=[
        "id",
        "associated_transcript_id",
        "junction_index",
        "chrom",
        "genomic_start",
        "genomic_end",
        "strand",
    ]
    def __init__(
        self,
        associated_transcript_id: str,
        junction_index: int,
        chrom: str,
        genomic_start: int,
        genomic_end_1based: int,
        strand: str,
        id = None,
    ):
        self.associated_transcript_id=associated_transcript_id
        self.junction_index=junction_index
        super().__init__(chrom=chrom,genomic_start=genomic_start,genomic_end_1based=genomic_end_1based,strand=strand)
        if id==None:
            self.id="JUNC__{}__{}".format(self.associated_transcript_id,"junction_"+str(self.junction_index))
    def __eq__(self,rhs:"RefTranscriptJunction"):
        assert isinstance(rhs,RefTranscriptJunction)
        return self.associated_transcript_id==rhs.associated_transcript_id and self.chrom==rhs.chrom and self.genomic_start==rhs.genomic_start and self.genomic_end_1based==rhs.genomic_end_1based and self.strand==rhs.strand

    def get_splice_sites(self):
        self.splice_donor=SpliceSite(
            ss_type="donor",
            chrom=self.chrom,
            genomic_loc=self.genomic_start if self.strand=="+" else self.genomic_end_1based-1,
            strand=self.strand,
            id="SD__{}__{}".format(self.associated_transcript_id,"junction_"+str(self.junction_index))
        )

        self.splice_acceptor=SpliceSite(
            ss_type="acceptor",
            chrom=self.chrom,
            genomic_loc=self.genomic_end_1based-1 if self.strand=="+" else self.genomic_start,
            strand=self.strand,
            id="SA__{}__{}".format(self.associated_transcript_id,"junction_"+str(self.junction_index))
        )
class RefTranscriptTSS(TSS):
    exported_fields=TSS.exported_fields+[
        "associated_transcript_id"
    ]

    def __init__(
        self,
        associated_transcript_id: str,
        chrom: str,
        genomic_center: int,
        strand: str,
        genomic_window_start: int = None,
        genomic_window_end_1based: int = None,
        id: str = None
    ):
        super().__init__(chrom,genomic_center,strand,genomic_window_start,genomic_window_end_1based,id)
        self.associated_transcript_id=associated_transcript_id
        if id==None:
            self.id="TSS__{}".format(self.associated_transcript_id)

    def center_equals(self,rhs:"RefTranscriptTSS",tol:int=0):
        assert isinstance(rhs,RefTranscriptTSS)
        return self.associated_transcript_id==rhs.associated_transcript_id \
            and self.chrom==rhs.chrom \
            and abs(self.genomic_center-rhs.genomic_center)<=tol \
            and self.strand==rhs.strand

class RefTranscriptTTS(TTS):
    exported_fields=TTS.exported_fields+[
        "associated_transcript_id"
    ]

    def __init__(
        self,
        associated_transcript_id: str,
        chrom: str,
        genomic_center: int,
        strand: str,
        genomic_window_start: int = None,
        genomic_window_end_1based: int = None,
        id: str = None
    ):
        super().__init__(chrom,genomic_center,strand,genomic_window_start,genomic_window_end_1based,id)
        self.associated_transcript_id=associated_transcript_id
        if id==None:
            self.id="TTS__{}".format(self.associated_transcript_id)

    def center_equals(self,rhs:"RefTranscriptTTS",tol:int=0):
        assert isinstance(rhs,RefTranscriptTTS)
        return self.associated_transcript_id==rhs.associated_transcript_id \
            and self.chrom==rhs.chrom \
            and abs(self.genomic_center-rhs.genomic_center)<=tol \
            and self.strand==rhs.strand

class RefTranscriptIsoform(TranscriptomicFeature):
    exported_fields=[
        "id",
        "num_exons",
        "genomic_span",
        "transcript_length",
        "chrom",
        "strand",
        "exon_intervals",
        "junction_intervals",
        "genes"
    ]
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
        genes = None
    ):
        super().__init__()
        self.id= id 
        self.genomic_start=genomic_start
        self.genomic_end_1based=genomic_end_1based
        self.chrom=chrom
        self.strand=strand
        self.genes=genes
        self.exons: Dict[int,RefTranscriptExon] = dict()
        self.junctions: Dict[int,RefTranscriptJunction] = dict()
        if isinstance(txStruct,list) or isinstance(txStruct,tuple):
            assert isinstance(txStruct[0],tuple), type(txStruct[0])
        elif isinstance(txStruct,genePredRecord):
            assert txStruct.chrom==self.chrom
            assert txStruct.strand==self.strand
            assert txStruct.txStart==txStruct.exonStarts[0]==self.genomic_start
            assert txStruct.txEnd==txStruct.exonEnds[-1]==self.genomic_end_1based
            txStruct=[(exon.start,exon.end) for exon in txStruct.exons]
        else:
            assert False

        for exon_genomic_index,(exon_genomic_start,exon_genomic_end_1based) in enumerate(txStruct):
            exon=RefTranscriptExon(
                associated_transcript_id=self.id,
                exon_index=exon_genomic_index,
                chrom=self.chrom,
                genomic_start=exon_genomic_start,
                genomic_end_1based=exon_genomic_end_1based,
                strand=self.strand
            )
            self.exons[exon_genomic_index]=exon

        for exon_genomic_index in range(len(txStruct)-1):
            junction_genomic_index=exon_genomic_index
            junction_genomic_start=txStruct[exon_genomic_index][1]
            junction_genomic_end=txStruct[exon_genomic_index+1][0]
            assert junction_genomic_start<junction_genomic_end
            junction=RefTranscriptJunction(
                associated_transcript_id=self.id,
                junction_index=junction_genomic_index,
                chrom=self.chrom,
                genomic_start=junction_genomic_start,
                genomic_end_1based=junction_genomic_end,
                strand=self.strand
            )
            self.junctions[junction_genomic_index]=junction

        self.CDS_genomic_start=CDS_genomic_start
        self.CDS_genomic_end_1based=CDS_genomic_end_1based
        self.get_ref_TSS_TTS()

    def get_ref_TSS_TTS(self):
        self.TSS=RefTranscriptTSS(
            associated_transcript_id=self.id,
            chrom=self.chrom,
            genomic_center=self.genomic_start if self.strand=="+" else self.genomic_end_1based-1,
            strand=self.strand
        )

        self.TTS=RefTranscriptTTS(
            associated_transcript_id=self.id,
            chrom=self.chrom,
            genomic_center=self.genomic_end_1based-1 if self.strand=="+" else self.genomic_start,
            strand=self.strand
        )

    @staticmethod
    def from_genePred_record(record:genePredRecord):
        return RefTranscriptIsoform(
            id=record.id,
            chrom=record.chrom,
            strand=record.strand,
            genomic_start=record.txStart,
            genomic_end_1based=record.txEnd,
            CDS_genomic_start=record.cdsStart,
            CDS_genomic_end_1based=record.cdsEnd,
            txStruct=record,
            genes=[record.gene]
        )
    @property
    def num_exons(self):
        return len(self.exons)

    @property
    def genomic_span(self):
        return self.genomic_end_1based-self.genomic_start
    
    @property
    def transcript_length(self):
        length=0
        for exon in self.exons.values():
            length+=exon.length
        return length
    
    @property
    def exon_intervals(self):
        return tuple((self.exons[k].genomic_start,self.exons[k].genomic_end_1based) for k in sorted(self.exons.keys()))
    
    @property
    def exon_intervals_with_chr_strand(self):
        return (self.chrom,self.strand,self.exon_intervals)

    @property
    def junction_intervals(self):
        return tuple((self.junctions[k].genomic_start,self.junctions[k].genomic_end_1based) for k in sorted(self.junctions.keys()))
    
    @property
    def junction_intervals_with_chr_strand(self):
        return (self.chrom,self.strand,self.junction_intervals)

    @property
    def transcript_structure(self):
        return TranscriptStructure(
            chrom=self.chrom,strand=self.strand,junction_intervals=self.exon_intervals,
            id="TS__{}".format(self.id)
        )
    def junction_equals(self,rhs:"RefTranscriptIsoform"):
        assert isinstance(rhs,RefTranscriptIsoform)
        return self.junction_intervals==rhs.junction_intervals
    
    def exon_equals(self,rhs:"RefTranscriptIsoform"):
        assert isinstance(rhs,RefTranscriptIsoform)
        return self.exon_intervals==rhs.exon_intervals

class TranscriptIsoformWithSeq(RefTranscriptIsoform):
    def __init__(
        self,
        id: str,
        chrom: str,
        genomic_start: int,
        genomic_end_1based: int,
        strand: str,
        txStruct,
        CDS_genomic_start: int = -1,
        CDS_genomic_end_1based: int = -1,
        genes = None,
        genome_dict = None
    ):
        super().__init__(
            id=id,
            chrom=chrom,
            genomic_start=genomic_start,
            genomic_end_1based=genomic_end_1based,
            strand=strand,
            txStruct=txStruct,
            CDS_genomic_start=CDS_genomic_start,
            CDS_genomic_end_1based=CDS_genomic_end_1based,
            genes=genes
        )

        chrom_seq=genome_dict[self.chrom]
        if isinstance(chrom_seq,Seq):
            self.genomic_seq = genome_dict[self.chrom].seq[self.genomic_start:self.genomic_end_1based]
        elif isinstance(chrom_seq,pyfaidx.FastaRecord):
            self.genomic_seq = Seq(genome_dict[self.chrom][self.genomic_start:self.genomic_end_1based].seq)
        elif isinstance(chrom_seq,str):
            self.genomic_seq = Seq(genome_dict[self.chrom][self.genomic_start:self.genomic_end_1based])
        if self.strand=="+":
            self.CDS_start_transcript=self.genomic_to_transcript(self.CDS_genomic_start,stick="5p")
            self.CDS_end_transcript_1based=self.genomic_to_transcript(self.CDS_genomic_end_1based-1,stick="5p")+1
        else:
            self.CDS_start_transcript=self.genomic_to_transcript(self.CDS_genomic_end_1based-1,stick="5p")
            if self.CDS_genomic_start==self.CDS_genomic_end_1based:
                self.CDS_end_transcript_1based=0
            else:
                self.CDS_end_transcript_1based=self.genomic_to_transcript(self.CDS_genomic_start,stick="5p")+1

    def get_coding_protein_sequence(self):
        if self.genomic_start<=self.CDS_genomic_start<self.CDS_genomic_end_1based<=self.genomic_end_1based:
            return self.seq_from_transcript_coords(self.CDS_start_transcript,self.CDS_end_transcript_1based).translate()
        else:
            return Seq("")

    def seq_from_genomic_coords(self, start, end = None):
        assert self.genomic_start <= start <= end <= self.genomic_end_1based
        if end == None:
            subseq = self.genomic_seq[start-self.genomic_start]
        else:
            subseq = self.genomic_seq[(start-self.genomic_start):(end-self.genomic_start)]
            
        if self.strand=="-":
            subseq = subseq.reverse_complement()
        return subseq
    
    def genomic_location_test(self,loc):
        if loc<self.genomic_start:
            return "OOB"
        elif self.genomic_start <= loc < self.genomic_end_1based:
            for exon_interval in self.exon_intervals:
                if exon_interval[0]<=loc<exon_interval[1]:
                    return "exon"
            return "intron"
        else:
            return "OOB"
            
    def seq_from_transcript_coords(self, start, end_1based = None, return_intervals = False):
        assert 0<=start<=end_1based<=self.transcript_length
        if self.strand=="+":
            genomic_start=self.transcript_to_genomic(start)
            genomic_end_1based=self.transcript_to_genomic(end_1based)
        elif self.strand=="-":
            genomic_start=self.transcript_to_genomic(end_1based)+1
            genomic_end_1based=self.transcript_to_genomic(start)+1

        intervals=list()
        for interval in self.exon_intervals:
            interval,length=interval_overlap(interval,(genomic_start,genomic_end_1based))
            if length > 0:
                intervals.append(interval)
        seq=Seq("")
        if self.strand=="+":
            for interval in intervals:
                seq+=self.seq_from_genomic_coords(interval[0],interval[1])
        elif self.strand=="-":
            for interval in reversed(intervals):
                seq+=self.seq_from_genomic_coords(interval[0],interval[1])
        if return_intervals:
            return seq, intervals
        else:
            return seq

    def transcript_to_genomic(self,loc):
        assert 0<=loc<=self.transcript_length
        if self.strand=="+":
            if loc == self.transcript_length:
                genomic_loc = self.genomic_end_1based
            else:
                exon_i = 0
                while loc >= self.exons[exon_i].length: 
                    loc -= self.exons[exon_i].length
                    exon_i += 1
                genomic_loc = self.exons[exon_i].genomic_start + loc
        else:
            if loc == self.transcript_length:
                genomic_loc = self.genomic_start - 1
            else:
                exon_i = self.num_exons - 1
                while loc >= self.exons[exon_i].length:
                    loc -= self.exons[exon_i].length
                    exon_i -= 1
                genomic_loc = self.exons[exon_i].genomic_end_1based - 1 - loc
        return genomic_loc

    def genomic_to_transcript(self,loc,stick):
        assert stick in ("5p","3p")
        if self.strand == "+":
            assert self.genomic_start<=loc<=self.genomic_end_1based
            if loc == self.genomic_end_1based:
                transcript_loc=self.transcript_length
            else:
                interval_boundaries=list()
                for interval in self.exon_intervals:
                    interval_boundaries.extend(interval)
                i = 0
                transcript_loc = 0
                while not interval_boundaries[i]<=loc<interval_boundaries[i+1]:
                    if i % 2 == 0: # the current interval is an exon
                        transcript_loc += interval_boundaries[i+1]-interval_boundaries[i]
                    i+=1
                    
                if i % 2 == 0: # loc in exon
                    transcript_loc += loc - interval_boundaries[i]
                else: # loc in junction
                    if stick == "5p":
                        transcript_loc -= 1
                    elif stick == "3p":
                        pass
        else:
            assert self.genomic_start-1<=loc<=self.genomic_end_1based-1
            if loc == self.genomic_start-1:
                transcript_loc=self.transcript_length
            else:
                interval_boundaries=list()
                for interval in self.exon_intervals:
                    interval_boundaries.extend(interval)
                i = len(interval_boundaries)-2
                transcript_loc = 0
                while not interval_boundaries[i]<=loc<interval_boundaries[i+1]:
                    if i % 2 ==0:
                        transcript_loc += interval_boundaries[i+1]-interval_boundaries[i]
                    i -= 1
                if i % 2 ==0 :
                    transcript_loc += interval_boundaries[i+1] - 1 - loc
                else:
                    if stick == "5p":
                        transcript_loc -= 1
                    elif stick == "3p":
                         pass

        return transcript_loc
        
    @staticmethod
    def from_genePred_record(record:genePredRecord,genome_dict):
        return TranscriptIsoformWithSeq(
            id=record.id,
            chrom=record.chrom,
            strand=record.strand,
            genomic_start=record.txStart,
            genomic_end_1based=record.txEnd,
            CDS_genomic_start=record.cdsStart,
            CDS_genomic_end_1based=record.cdsEnd,
            txStruct=record,
            genes=[record.gene],
            genome_dict=genome_dict
        )
    
    @staticmethod
    def from_isoform(isoform:RefTranscriptIsoform, genome_dict):
        return TranscriptIsoformWithSeq(
            id = isoform.id,
            chrom = isoform.chrom,
            strand = isoform.strand,
            genomic_start = isoform.genomic_start,
            genomic_end_1based= isoform.genomic_end_1based,
            CDS_genomic_start = isoform.CDS_genomic_start,
            CDS_genomic_end_1based = isoform.CDS_genomic_end_1based,
            txStruct=isoform.exon_intervals,
            genes=isoform.genes,
            genome_dict=genome_dict
        )
        
    def protein_to_transcript(self,loc):
        pass

    def transcript_to_protein(self,loc):
        pass

    def protein_to_genomic(self,loc):
        pass