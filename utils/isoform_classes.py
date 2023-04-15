import pandas as pd
from math import isnan
from bx.intervals import Interval, IntervalTree
from typing import Dict, List
from collections import defaultdict
from Bio.Seq import Seq
import pyfaidx
class genePredReader(object):
    def __init__(self, filename):
        self.f = open(filename)

    def __iter__(self):
        return self

    def __next__(self):
        line = self.f.readline().strip()
        if len(line) == 0:
            raise StopIteration
        return genePredRecord.from_line(line)

    def update_transcript_isoform(self,isoform_dict: Dict[str,"QueryTranscriptIsoform"], use_existing_cds=True):
        update_isoform_dict={rec.id:rec for rec in self}
        key_intersection=set(isoform_dict.keys()).intersection(update_isoform_dict.keys())
        print("In existing isoform collection:",len(isoform_dict))
        print("In genePred:",len(update_isoform_dict))
        print("To update:",len(key_intersection))
        for k in key_intersection:
            isoform_dict[k].update_fields(dict(txStruct=update_isoform_dict[k]))
            if use_existing_cds:
                isoform_dict[k].txStruct.cdsStart = isoform_dict[k].CDS_start
                isoform_dict[k].txStruct.cdsEnd = isoform_dict[k].CDS_end
        return isoform_dict


class genePredRecord(object):
    def __init__(self, id, chrom, strand, txStart, txEnd, cdsStart, cdsEnd, exonCount, exonStarts, exonEnds, gene=None):
        self.id = id
        self.chrom = chrom
        self.strand = strand
        self.txStart = txStart         # 0-based start [original comment is incorrect]
        self.txEnd = txEnd             # 1-based end
        self.cdsStart = cdsStart       # 0-based start [original comment is incorrect]
        self.cdsEnd = cdsEnd           # 1-based end
        self.exonCount = exonCount
        self.exonStarts = exonStarts   # 0-based starts
        self.exonEnds = exonEnds       # 1-based ends
        self.gene = gene

        self.length = 0
        self.exons = []

        for s,e in zip(exonStarts, exonEnds):
            self.length += e-s
            self.exons.append(Interval(s, e))

        # junctions are stored (1-based last base of prev exon, 1-based first base of next exon)
        self.junctions = [(self.exonEnds[i],self.exonStarts[i+1]) for i in range(self.exonCount-1)]

    @property
    def segments(self):
        return self.exons


    @classmethod
    def from_line(cls, line):
        raw = line.strip().split('\t')
        return cls(id=raw[0],
                  chrom=raw[1],
                  strand=raw[2],
                  txStart=int(raw[3]),
                  txEnd=int(raw[4]),
                  cdsStart=int(raw[5]),
                  cdsEnd=int(raw[6]),
                  exonCount=int(raw[7]),
                  exonStarts=[int(x) for x in raw[8][:-1].split(',')],  #exonStarts string has extra , at end
                  exonEnds=[int(x) for x in raw[9][:-1].split(',')],     #exonEnds string has extra , at end
                  gene=raw[11] if len(raw)>=12 else None,
                  )

    def get_splice_site(self, genome_dict, i):
        """
        Return the donor-acceptor site (ex: GTAG) for the i-th junction
        :param i: 0-based junction index
        :param genome_dict: dict of chrom --> SeqRecord
        :return: splice site pattern, ex: "GTAG", "GCAG" etc
        """
        assert 0 <= i < self.exonCount-1

        d = self.exonEnds[i]
        a = self.exonStarts[i+1]

        seq_d = genome_dict[self.chrom].seq[d:d+2]
        seq_a = genome_dict[self.chrom].seq[a-2:a]

        if self.strand == '+':
            return (str(seq_d)+str(seq_a)).upper()
        else:
            return (str(seq_a.reverse_complement())+str(seq_d.reverse_complement())).upper()

class SQANTICLSReader(object):
    FIELDS_CLASS_MAPPING = {
        'isoform': ("id",str),
         'chrom': ("chrom",str),
         'strand': ("strand",str),
         'length': ("length",int),
         'exons': ("num_exons",int),
         'structural_category': ("str_class",str),
         'associated_gene': ("genes",list), # needs special treatment
         'associated_transcript': ("transcripts",list), # needs special treatment
         'ref_length': ("refLen",int),
         'ref_exons': ("refExons",int),
         'diff_to_TSS': ("tss_diff",int),
         'diff_to_TTS': ("tts_diff",int),
         'diff_to_gene_TSS': ("tss_gene_diff",int),
         'diff_to_gene_TTS': ("tts_gene_diff",int),
         'subcategory': ("subtype",str),
         'RTS_stage': ("RT_switching",bool),
         'all_canonical': ("canonical",bool),
         'min_sample_cov': ("min_samp_cov",float),
         'min_cov': ("min_cov",float),
         'min_cov_pos': ("min_cov_pos",int),
         'sd_cov': ("sd",float),
         'FL': ("FL",int),
         'n_indels': ("nIndels",int),
         'n_indels_junc': ("nIndelsJunc",int),
         'bite': ("bite",bool),
         'iso_exp': ("isoExp",float),
         'gene_exp': ("geneExp",float),
         'FSM_class': ("FSM_class",str),
         'coding': ("coding",str),
         'ORF_length': ("ORFlen",int),
         'ORF_seq': ("ORFseq",str),
         'CDS_start': ("CDS_start",int),
         'CDS_end': ("CDS_end",int),
         'CDS_genomic_start': ("CDS_genomic_start",int),
         'CDS_genomic_end': ("CDS_genomic_end",int),
         'predicted_NMD': ("is_NMD",bool),
         'perc_A_downstream_TTS': ("percAdownTTS",float),
         'seq_A_downstream_TTS': ("seqAdownTTS",str),
         'dist_to_cage_peak': ("dist_cage",int),
         'within_cage_peak': ("within_cage",bool),
         'dist_to_polya_site': ("dist_polya_site",int),
         'within_polya_site': ("within_polya_site",bool),
         'polyA_motif': ("polyA_motif",str),
         'polyA_dist': ("polyA_dist",int),
    }
    def __init__(self, filename):
        self.df = pd.read_csv(filename,sep="\t")
        self.df=self.df[list(SQANTICLSReader.FIELDS_CLASS_MAPPING.keys())]
        self.df["associated_gene"]=self.df["associated_gene"].str.split("__")
        self.df["associated_transcript"]=self.df["associated_transcript"].str.split("__")
        self.df["CDS_genomic_start"]-=1 # make 0-based
        self.df['CDS_start']-=1 # make 0-based

    def update_transcript_isoform(self,isoform_dict=None) -> Dict[str,"QueryTranscriptIsoform"]:
        def parse_dict(d):
            ret_d=dict()
            for k,v in d.items():
                internal_name,internal_dtype=SQANTICLSReader.FIELDS_CLASS_MAPPING[k]
                if type(v)==float and isnan(v):
                    if internal_dtype==str:
                        v="NA"
                    if internal_dtype in (int,float,bool):
                        v=float("nan")
                else:
                    v=internal_dtype(v)
                ret_d[internal_name]=v
            return ret_d
        if isoform_dict is None: # create new mode
            isoform_dict=dict()
            for _,row in self.df.iterrows():
                d=parse_dict(row.to_dict())
                isoform=QueryTranscriptIsoform()
                isoform.update_fields(d)
                isoform_dict[isoform.id]=isoform
        else: # update mode
            for _,row in self.df.iterrows():
                d=parse_dict(row.to_dict())
                assert d["id"] in isoform_dict
                isoform_dict[d["id"]].update_fields(d)
        return isoform_dict

class JunctionReader:
    FIELDS_CLASS_MAPPING = {
        'isoform': ("isoform",str),
        'junction_number': ("junction_number",int), # need special treatment
        "chrom": ("chrom",str),
        "strand": ("strand",str),
        "genomic_start_coord": ("genomic_start_coord",int),  # need special treatment to convert to 0-based
        "genomic_end_coord": ("genomic_end_coord",int),      # already is 1-based end
        "junction_category": ("junction_category",str),
        "start_site_category": ("start_site_category",str),
        "end_site_category": ("end_site_category",str),
        "diff_to_Ref_start_site": ("diff_to_Ref_start_site",str),
        "diff_to_Ref_end_site": ("diff_to_Ref_end_site",str),
        "bite_junction": ("bite_junction",bool),
        "splice_site": ("splice_site",str),
        "canonical": ("canonical",str),
        "RTS_junction": ("RTS_junction",bool), # First write ???? in _tmp, later is TRUE/FALSE
        "indel_near_junct": ("indel_near_junct",str),
        "phyloP_start": ("phyloP_start",str),
        "phyloP_end": ("phyloP_end",str),
        "sample_with_cov": ("sample_with_cov",int),
        "total_coverage_unique": ("total_coverage_unique",float),
        "total_coverage_multi": ("total_coverage_multi",float)
    }
    def __init__(self,filename):
        self.filename=filename
        self.df = pd.read_csv(filename,sep="\t")
        self.df=self.df[list(JunctionReader.FIELDS_CLASS_MAPPING.keys())]
        self.df["junction_number"]=self.df["junction_number"].map(lambda x: int(x.split('_')[1])-1).astype("int64")
        self.df["genomic_start_coord"]-=1 # make 0-based
    def update_transcript_isoform(self,isoform_dict: Dict[str,"QueryTranscriptIsoform"]):
        def parse_dict(d):
            ret_d=dict()
            for k,v in d.items():
                internal_name,internal_dtype=JunctionReader.FIELDS_CLASS_MAPPING[k]
                if type(v)==float and isnan(v):
                    if internal_dtype==str:
                        v="NA"
                    if internal_dtype in (int,float,bool):
                        v=float("nan")
                else:
                    v=internal_dtype(v)
                ret_d[internal_name]=v
            return ret_d

        key_intersection=set(isoform_dict.keys()).intersection(self.df["isoform"])
        print("In existing isoform collection:",len(isoform_dict))
        print("In Junction File:",self.df["isoform"].nunique())
        print("To update:",len(key_intersection))

        for _,row in self.df.iterrows():
            d=parse_dict(row.to_dict())
            junction_obj=QueryJunction()
            junction_obj.update_fields(d)
            if d["isoform"] in isoform_dict:
                isoform_dict[d["isoform"]].update_junction_info(junction_obj)
        
        return isoform_dict
        
class CAGEPeak:
    def __init__(self, cage_bed_filename):
        self.cage_bed_filename = cage_bed_filename
        self.cage_peaks = defaultdict(lambda: IntervalTree()) # (chrom,strand) --> intervals of peaks

        self.read_bed()

    def read_bed(self):
        for line in open(self.cage_bed_filename):
            raw = line.strip().split()
            chrom = raw[0]
            start0 = int(raw[1])
            end1 = int(raw[2])
            strand = raw[5]
            tss0 = int(raw[6])
            self.cage_peaks[(chrom,strand)].insert(start0, end1, (tss0, start0, end1))

    def find(self, chrom, strand, query, search_window=10000):
        """
        :param start0: 0-based start of the 5' end to query
        :return: <True/False falls within a cage peak>, <nearest dist to TSS>
        dist to TSS is 0 if right on spot
        dist to TSS is + if downstream, - if upstream (watch for strand!!!)
        """
        within_peak, dist_peak = False, 'NA'
        for (tss0,start0,end1) in self.cage_peaks[(chrom,strand)].find(query-search_window, query+search_window):
 # Skip those cage peaks that are downstream the detected TSS because degradation just make the transcript shorter
            if strand=='+' and start0>int(query) and end1>int(query):
                continue
            if strand=='-' and start0<int(query) and end1<int(query):
                continue
##
            if not within_peak:
                within_peak, dist_peak = (start0<=query<end1), (query - tss0) * (-1 if strand=='-' else +1)
            else:
                d = (query - tss0) * (-1 if strand=='-' else +1)
                if abs(d) < abs(dist_peak):
                    within_peak, dist_peak = (start0<=query<end1), d
        return within_peak, dist_peak

class PolyAPeak:
    def __init__(self, polya_bed_filename):
        self.polya_bed_filename = polya_bed_filename
        self.polya_peaks = defaultdict(lambda: IntervalTree()) # (chrom,strand) --> intervals of peaks

        self.read_bed()

    def read_bed(self):
        for line in open(self.polya_bed_filename):
            raw = line.strip().split()
            chrom = raw[0]
            # start0 = int(raw[1])
            # end1 = int(raw[2])
            polyA_peak_loc=int(raw[3].split(":")[1])-1
            start0=polyA_peak_loc
            end1=polyA_peak_loc+1

            strand = raw[5]
            self.polya_peaks[(chrom,strand)].insert(start0, end1, (start0, end1))

    def find(self, chrom, strand, query, search_window=100):
        """
        :param start0: 0-based start of the 5' end to query
        :return: <True/False falls within some distance to polyA>, distance to closest
        + if downstream, - if upstream (watch for strand!!!)
        """
        assert strand in ('+', '-')
        hits = self.polya_peaks[(chrom,strand)].find(query-search_window, query+search_window)
        if len(hits) == 0:
            return False, None
        else:
            s0, e1 = hits[0]
            min_dist = query - s0
            for s0, e1 in hits[1:]:
                d = query - s0
                if abs(d) < abs(min_dist):
                    min_dist = d
            if strand == '-':
                min_dist = -min_dist
            return True, min_dist

def reference_genePred_TSS_TTS_parser(referenceFiles):
    known_TSS_TTS_lists_by_gene = defaultdict(lambda: {'TSS':list(), 'TTS': list()})
    num_dropped_transcript_isoforms=0
    for i,r in enumerate(genePredReader(referenceFiles)):
        if i%100==0:
            print(f"[{i}/?]",end="\r")
        if r.length < 200:
            num_dropped_transcript_isoforms+=1
            continue
        ref_ti=RefTranscriptIsoform(
            id=r.id,
            chrom=r.chrom,
            genomic_start=r.txStart,
            genomic_end_1based=r.txEnd,
            strand=r.strand,
            txStruct=r,
            CDS_genomic_start=r.cdsStart,
            CDS_genomic_end_1based=r.cdsEnd,
            genes=[r.gene]
        )

        known_TSS_TTS_lists_by_gene[r.gene]["TSS"].append(ref_ti.TSS)
        known_TSS_TTS_lists_by_gene[r.gene]["TTS"].append(ref_ti.TTS)
    
    return known_TSS_TTS_lists_by_gene

if __name__=="__main__":
    pass
    

    