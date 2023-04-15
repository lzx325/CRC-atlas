import re
def get_read_info(alignment):
    if alignment.has_tag("CB"):
        cell_barcode=alignment.get_tag("CB")
        cell_barcode_from="CB"
    elif alignment.has_tag("XC"):
        cell_barcode=alignment.get_tag("XC")
        cell_barcode_from="XC"
    elif alignment.has_tag("CR"):
        cell_barcode=alignment.get_tag("CR")
        cell_barcode_from="CR"
    else:
        cell_barcode=""
        cell_barcode_from=""
    if alignment.has_tag("UB"):
        UMI=alignment.get_tag("UB")
        UMI_from="UB"
    elif alignment.has_tag("XM"):
        UMI=alignment.get_tag("XM")
        UMI_from="XM"
    elif alignment.has_tag("UR"):
        UMI=alignment.get_tag("UR")
        UMI_from="UR"
    else:
        UMI=""
        UMI_from=""
    if alignment.has_tag("UY"):
        UMI_qual=alignment.get_tag("UY")
    else:
        UMI_qual=""
    cell_barcode=re.sub("-\d+","",cell_barcode)
    UMI=re.sub("-\d+","",UMI)
    return {
        "read_name":alignment.query_name,
        "cell_barcode":cell_barcode,
        "UMI":UMI,
        "cell_barcode_from":cell_barcode_from,
        "UMI_from":UMI_from,
        "UMI_qual":UMI_qual,
        "mapping_quality":alignment.mapping_quality,
        "flag":alignment.flag
    }