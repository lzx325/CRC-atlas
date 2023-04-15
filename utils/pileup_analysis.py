from collections import namedtuple

from utils.read_utils import get_read_info


GenomicRegion=namedtuple("GenomicRegion",["chr","start","end"])
GenomicLocation=namedtuple("GenomicLocation",["chr","loc"])
SNP=namedtuple("SNP",["loc","type","op"])
BAM_FUNMAP=4
BAM_FSECONDARY=256
BAM_FQCFAIL=512
BAM_FDUP=1024

def assert_valid_region(region):
    assert region.end>region.start
    assert region.start>=0

def get_expected_sequence(ref,snp,g_region):
    assert_valid_region(g_region)
    assert 0<=snp.loc.loc<len(ref[snp.loc.chr])
    assert g_region.end<=len(ref[g_region.chr])
    seq=ref[g_region.chr][g_region.start:g_region.end]
    if snp.loc.chr==g_region.chr:
        op_start=snp.loc.loc-g_region.start
        if snp.type=="SNP":
            if op_start<len(seq):
                oldnt,newnt=snp.op.split(">")
                seq=str(seq)
                assert seq[op_start]==oldnt
                seq=seq[:op_start]+newnt+seq[(op_start+1):]
            return seq
        elif snp.type=="INS":
            if op_start<len(seq):
                oldnt,newnt=snp.op.split(">")
                seq=str(seq)
                assert oldnt=="-"
                seq=seq[:(op_start+1)]+newnt+seq[(op_start+1):]
            return seq
        elif snp.type=="DEL":
            if op_start<len(seq):
                
                oldnt,newnt=snp.op.split(">")
                op_end=op_start+len(oldnt)
                seq=str(seq)
                assert str(ref[snp.loc.chr][snp.loc.loc:(snp.loc.loc+len(oldnt))])==oldnt
                assert newnt=="-"
                op_end=min(op_end,len(seq))
                
                seq=seq[:(op_start)]+seq[op_end:]
            return seq
        elif snp.type=="DNP":
            if op_start<len(seq):
                oldnt,newnt=snp.op.split(">")
                op_end=op_start+len(oldnt)
                seq=str(seq)
                assert str(ref[snp.loc.chr][snp.loc.loc:(snp.loc.loc+len(oldnt))])==oldnt
                op_end=min(op_end,len(seq))
                seq=seq[:(op_start)]+newnt+seq[op_end:]
            return seq


    
def get_supporting_reads(samfile,ref,snp,strand=None):
    result_dict={
        "read_name":list(),
        "cell_barcode":list(),
        "UMI":list(),
        "cell_barcode_from":list(),
        "UMI_from":list(),
        "UMI_qual":list(),
        "mapping_quality":list(),
        "flag":list()
    }
    assert strand in (None,"+","-")
    if snp.type=="SNP":
        oldnt,newnt=snp.op.split(">")
        assert ref[snp.loc.chr][snp.loc.loc]==oldnt
        for pileupcolumn in samfile.pileup(
                        snp.loc.chr, 
                        start=snp.loc.loc, 
                        stop=snp.loc.loc+1,
                        truncate=True,
                        min_base_quality=0,
                        ignore_orphans=False,
                        ignore_overlaps=False,
                        flag_filter=0
#                         fastafile=fastafile
                    ):
            read_count=0
            for pileupread in pileupcolumn.pileups:
                if not pileupread.is_del \
                and not pileupread.is_refskip \
                and (strand is None or (strand == "+" and not pileupread.alignment.is_reverse) or (strand == "-" and pileupread.alignment.is_reverse)) \
                and pileupread.alignment.query_sequence[pileupread.query_position]==newnt:
                    info=get_read_info(pileupread.alignment)
                    [result_dict[k].append(v) for k,v in info.items()]

    elif snp.type=="INS" or snp.type=="DEL":
        oldnt,newnt=snp.op.split(">")
        if snp.type=="INS":
            assert oldnt=="-"
        else:
            assert oldnt==ref[snp.loc.chr][snp.loc.loc:(snp.loc.loc+len(oldnt))] and newnt=="-"
        look_around=10
        candidates=list()
        if snp.type=="INS":
            candidates.append(snp.loc.loc)
            cur=snp.loc.loc
            while cur>=snp.loc.loc-look_around and newnt==ref[snp.loc.chr][(cur-len(newnt)):cur]:
                cur-=len(newnt)
                candidates.append(cur)
            cur=snp.loc.loc
            while cur<snp.loc.loc+1+look_around and newnt==ref[snp.loc.chr][(cur+len(newnt)):(cur+2*len(newnt))]:
                cur+=len(newnt)
                candidates.append(cur)
        elif snp.type=="DEL":
            candidates.append(snp.loc.loc-1)
            cur=snp.loc.loc
            while cur>=snp.loc.loc-look_around and oldnt==ref[snp.loc.chr][(cur-len(newnt)):cur]:
                cur-=len(newnt)
                candidates.append(cur-1)
            cur=snp.loc.loc
            while cur<snp.loc.loc+1+look_around and oldnt==ref[snp.loc.chr][(cur+len(newnt)):(cur+2*len(newnt))]:
                cur+=len(newnt)
                candidates.append(cur-1)
        candidates=list(sorted(candidates))
        
        for pileupcolumn in samfile.pileup(
                        snp.loc.chr, 
                        start=min(candidates), 
                        stop=max(candidates)+1,
                        truncate=True,
                        min_base_quality=0,
                        ignore_orphans=False,
                        ignore_overlaps=False,
                        flag_filter=0
#                         fastafile=fastafile
                    ):
            
            if pileupcolumn.reference_pos in candidates:
                for i,pileupread in enumerate(pileupcolumn.pileups):
                    
                    if (
                        snp.type=="INS" and pileupread.indel==len(newnt) \
                        or snp.type=="DEL" and pileupread.indel==-len(oldnt) \
                    ) \
                    and not pileupread.is_del \
                    and not pileupread.is_refskip \
                    and (strand is None or (strand == "+" and pileupread.alignment.is_forward) or (strand == "-" and pileupread.alignment.is_reverse)) \
                    and (
                        snp.type=="INS" and pileupread.alignment.query_sequence[(pileupread.query_position+1):(pileupread.query_position+1+len(newnt))] == newnt \
                        or snp.type=="DEL" \
                    ):
                        # if pileupread.alignment.qname=="A01018:86:HG7GNDSXY:2:2245:21215:22639":
                        #     import ipdb; ipdb.set_trace()
                        info=get_read_info(pileupread.alignment)
                        
                        [result_dict[k].append(v) for k,v in info.items()]
    return result_dict