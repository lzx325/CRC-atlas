import pandas as pd
def line_count(fp):
    line_count=0
    with open(fp) as f:
        for line in f:
            line_count+=1
    return line_count

def iprint(s="",*args,ident=0,**kwargs):
	print('\t'*ident+s,*args,**kwargs)

def write_bed(bed_df,fp,gffTags=False):
    with open(fp,'w') as f:
        if gffTags:
            f.write("#gffTags\n")
        bed_df.to_csv(f,header=None,sep="\t",index=None)

def read_bed(fp):
    bed_df=pd.read_csv(fp,sep="\t",header=None)
    bed_df.columns=["chrom","start","end","name","score","strand","thick_start", "thick_end","color"]
    return bed_df

