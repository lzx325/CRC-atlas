def canonical_chr_regex():
    return r'^chr(\d+|X|Y|M)$'

def interval_overlap(interval1,interval2):
    start=max(interval1[0],interval2[0])
    end=min(interval1[1],interval2[1])
    return (start,end),end-start 