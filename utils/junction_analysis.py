from Bio.Seq import Seq
import pysam
from argparse import Namespace
import os
import pandas as pd
import numpy as np
from utils.read_utils import get_read_info

def find_intron(samfile,start,length,ROI):
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
    for i,seq in enumerate(samfile.fetch(ROI.ref,ROI.start,ROI.end)):
        pos=seq.pos
        for tup in seq.cigartuples:
            if pos==start and tup[0]==3 and tup[1]==length:
                info=get_read_info(seq)
                [result_dict[k].append(v) for k,v in info.items()]
            if tup[0] in [0,2,3,7,8]:
                pos+=tup[1]
    return result_dict