from collections import defaultdict

def compute_subset_weight(subset,identity_weights):
    return sum(identity_weights[i] for i in subset)

def select_best_kmer(candidate_kmers,kmer_binding_sets,covered_sets,identity_weights):
    if candidate_kmers:
        best_kmer=None
        best_weight_gain=0
        for kmer in candidate_kmers:
            new_kmer_binding_set=kmer_binding_sets[kmer]
            difference_set=new_kmer_binding_set-covered_sets
            weight_gain=compute_subset_weight(difference_set,identity_weights)
            if weight_gain>best_weight_gain:
                best_weight_gain=weight_gain
                best_kmer=kmer
    else:
        best_kmer,best_weight_gain=None,0
    return best_kmer,best_weight_gain

def select_kmer_greedy_iterative(kmer_binding_sets,identity_weights,n_kmers=None, kmers_by_transcripts=None, n_kmers_per_transcript=3):
    identity_set=set()
    for s in kmer_binding_sets.values():
        identity_set.update(s)
    assert all(s in identity_weights for s in identity_set)
    
    if kmers_by_transcripts is not None:
        kmer_set=set()
        for k in kmers_by_transcripts.values():
            kmer_set.update(k)
        assert all(k in kmer_binding_sets for k in kmer_set)
    
    print("# kmers:",len(kmer_binding_sets))
    print("# identities:",len(identity_weights))
    
    selected_kmers=list()
    unselected_kmers=list(kmer_binding_sets.keys())
    covered_sets=set()
    if kmers_by_transcripts is not None:
        selected_kmers_by_isoform=defaultdict(list)
    
    if n_kmers is None:
        n_kmers=n_kmers_per_transcript*len(kmers_by_transcripts)

    if kmers_by_transcripts is not None:
        transcripts_list=list(kmers_by_transcripts)
    for i in range(n_kmers):
        if kmers_by_transcripts is None:
            kmer,weight_gain=select_best_kmer(
                unselected_kmers,
                kmer_binding_sets,
                covered_sets,
                identity_weights
            )
        else:
            transcript_id=transcripts_list[i//n_kmers_per_transcript]
            candidate_kmers=[kmer for kmer in sorted(kmers_by_transcripts[transcript_id]) if kmer in unselected_kmers]
            
            kmer,weight_gain=select_best_kmer(
                candidate_kmers,
                kmer_binding_sets,
                covered_sets,
                identity_weights
            )

        if kmer is not None:
            if kmers_by_transcripts is not None:
                print("From {} selecting kmer: {}, weight_gain: {}".format(transcript_id,kmer,weight_gain))
                selected_kmers_by_isoform[transcript_id].append(kmer)
            else:
                print("selecting kmer: {}, weight_gain: {}".format(kmer,weight_gain))
            unselected_kmers.remove(kmer)
            selected_kmers.append(kmer)
            covered_sets=covered_sets|kmer_binding_sets[kmer]
        if len(covered_sets)==len(identity_weights):
            print("covered sets full, resetting")
            covered_sets=set()
            
    final_covered_sets=set()
    for kmer in selected_kmers:
        final_covered_sets.update(kmer_binding_sets[kmer])
        
    covered_sets_weight=compute_subset_weight(final_covered_sets,identity_weights)
    print("weight of covered set:",covered_sets_weight)
    
    if kmers_by_transcripts is not None:
        selected_kmers_by_isoform=[(k,v) for k,v in selected_kmers_by_isoform.items()]
        return selected_kmers_by_isoform,covered_sets
    else:
        return selected_kmers,covered_sets