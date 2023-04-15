import pandas as pd
import re
import csv
from collections import defaultdict
from tqdm import tqdm

def my_gtf_read(gtf_fp,skiprows=None):
    table=pd.read_csv(gtf_fp,comment="#",sep="\t",header=None,skiprows=skiprows)
    table.columns=['seqname', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'tags']
    def parse_gtf_tags(x):
        kv_dict=dict()
        for kv in x.split("; "):
            res=re.match(r'(\w+) "(.*)"',kv)
            if res:
                k=res.group(1)
                v=res.group(2)
                kv_dict[k]=v
        return kv_dict
    table["tags"]=table["tags"].map(parse_gtf_tags)
    return table

def my_gtf_write(table,gtf_fp,mode="gtf",header=None,escape_characters={"+":"%2B",";":"%3B","=":"%3D","&":"%26",",":"%2C"}):
    table=table.copy()
    if header==None:
        header="# Produced by function my_gtf_write\n"

    def escape(s,escape_characters):
        for c in escape_characters:
            s=s.replace(c,escape_characters[c])
        return s

    def deparse_gtf_tags(kv_dict,mode="gtf"):
        assert mode in ("gtf","gff3")
        if mode=="gtf":
            kv_dict_to_string=""
            for k in sorted(kv_dict.keys()):
                v=kv_dict[k]
                kv_dict_to_string+=f'{k} "{v}"; '
            kv_dict_to_string=kv_dict_to_string.rstrip(' ')
        elif mode=="gff3":
            kv_dict_to_string=""
            for k in sorted(kv_dict.keys()):
                v=kv_dict[k]
                k=escape(str(k),escape_characters=escape_characters)
                v=escape(str(v),escape_characters=escape_characters)
                kv_dict_to_string+=f'{k}={v}; '
            kv_dict_to_string=kv_dict_to_string.rstrip(' ')
        return kv_dict_to_string
        
    table["tags"]=table["tags"].map(lambda x: deparse_gtf_tags(x,mode=mode))
    with open(gtf_fp,'w') as f:
        f.write(header)
        table.to_csv(f,sep="\t",header=None,index=None,quoting=csv.QUOTE_NONE)

class CategorizedGTF(object):
    def __init__(self):
        self.categorized_gtf=defaultdict(
            lambda: dict(
                gene_metadata=None,
                transcripts=defaultdict(
                    lambda: dict(
                        transcript_metadata=None,
                        transcript_feature_collections=defaultdict(list)
                    )
                )
            )
        )
    def update(self,parsed_gtf_df: pd.DataFrame, update_mode="everything"):
        assert update_mode in ("everything","assert_existing","ignore_unknown")
        skipped_row_index_list=list()
        for index, row in tqdm(parsed_gtf_df.iterrows(),total=parsed_gtf_df.shape[0]):
            gene_key=(row["seqname"],row["strand"],row["tags"]["gene_id"])
            if row["feature"]=="gene":
                gene_id = row["tags"]["gene_id"]
                if update_mode=="everything" or gene_key in self.categorized_gtf:
                    assert self.categorized_gtf[gene_key]["gene_metadata"] is None
                    self.categorized_gtf[gene_key]["gene_metadata"]=row
                else:
                    if update_mode=="assert_existing":
                        assert False, "{} is not in current collection".format(gene_key)
                    skipped_row_index_list.append(index)

            elif row["feature"]=="transcript":
                transcript_id=row["tags"]["transcript_id"]
                if update_mode=="everything" or gene_key in self.categorized_gtf and transcript_id in self.categorized_gtf[gene_key]["transcripts"]:
                    assert self.categorized_gtf[gene_key]["transcripts"][transcript_id]["transcript_metadata"] is None
                    self.categorized_gtf[gene_key]["transcripts"][transcript_id]["transcript_metadata"]=row
                else:
                    if update_mode=="assert_existing":
                        if gene_key not in self.categorized_gtf:
                            assert False, "{} is not in current collection".format(gene_key)
                        elif transcript_id in self.categorized_gtf[gene_key]["transcripts"]:
                            assert False, "{} is not in current collection".format(transcript_id)
                    skipped_row_index_list.append(index)

            elif row["feature"] in ["exon","CDS","UTR","start_codon","stop_codon","Selenocysteine"]:
                feature_name=row["feature"]
                transcript_id=row["tags"]["transcript_id"]

                if update_mode=="everything" or gene_key in self.categorized_gtf and transcript_id in self.categorized_gtf[gene_key]["transcripts"]:
                    self.categorized_gtf[gene_key]["transcripts"][transcript_id]["transcript_feature_collections"][feature_name].append(row)
                else:
                    if update_mode=="assert_existing":
                        if gene_key not in self.categorized_gtf:
                            assert False, "{} is not in current collection".format(gene_key)
                        elif transcript_id in self.categorized_gtf[gene_key]["transcripts"]:
                            assert False, "{} is not in current collection".format(transcript_id)
                    skipped_row_index_list.append(index)
        gene_start_updated=list()
        gene_end_updated=list()
        for gene_key in self.categorized_gtf.keys():
            if len(self.categorized_gtf[gene_key]["transcripts"])>0:
                transcript_start_min=min(
                    transcript["transcript_metadata"]["start"] for transcript in self.categorized_gtf[gene_key]["transcripts"].values()
                )
                transcript_end_max=max(
                    transcript["transcript_metadata"]["end"] for transcript in self.categorized_gtf[gene_key]["transcripts"].values()
                )

                if self.categorized_gtf[gene_key]["gene_metadata"] is not None:
                    if transcript_start_min<self.categorized_gtf[gene_key]["gene_metadata"]["start"]:
                        self.categorized_gtf[gene_key]["gene_metadata"]["start"]=transcript_start_min
                        gene_start_updated.append(gene_key)
                    if transcript_end_max>self.categorized_gtf[gene_key]["gene_metadata"]["end"]:
                        self.categorized_gtf[gene_key]["gene_metadata"]["end"]=transcript_end_max
                        gene_end_updated.append(gene_key)

        print("update finished: {} rows skipped, {} gene starts updated, {} gene ends updated".format(
            len(skipped_row_index_list),len(gene_start_updated),len(gene_end_updated)))

        return skipped_row_index_list


    def concatenate(self):
        concatenated_df_row_list=list()
        
        for gene_key, gene_data in tqdm(self.categorized_gtf.items()):
            gene_metadata=gene_data["gene_metadata"]
            if gene_metadata is not None:
                concatenated_df_row_list.append(gene_metadata)
            for transcript_id, transcript_data in gene_data["transcripts"].items():
                transcript_metadata=transcript_data["transcript_metadata"]
                if transcript_metadata is not None:
                    concatenated_df_row_list.append(transcript_metadata)
                for feature_name, feature_data in transcript_data["transcript_feature_collections"].items():
                    concatenated_df_row_list.extend(feature_data)
        concatenated_df=pd.concat(concatenated_df_row_list,axis=1,ignore_index=True).T
        return concatenated_df

if __name__=="__main__":
    tmp_gtf1=my_gtf_read("tmp_gtf/tmp.gtf")
    tmp_gtf2=my_gtf_read("tmp_gtf/tmp2.gtf")

    cat_gtf=CategorizedGTF()
    cat_gtf.update(tmp_gtf1)
    cat_gtf.update(tmp_gtf2,update_mode="ignore_unknown")
    concat_table=cat_gtf.concatenate()
    my_gtf_write(concat_table,"tmp_gtf/tmp_combined.gtf")
