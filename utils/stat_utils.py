import os
import sys
import time
import contextlib
import joblib
import numpy as np
from scipy.stats import fisher_exact
from tqdm import tqdm
from threadpoolctl import threadpool_limits
@contextlib.contextmanager
def tqdm_joblib(tqdm_object):
    """Context manager to patch joblib to report into tqdm progress bar given as argument"""
    class TqdmBatchCompletionCallback(joblib.parallel.BatchCompletionCallBack):
        def __call__(self, *args, **kwargs):
            tqdm_object.update(n=self.batch_size)
            return super().__call__(*args, **kwargs)

    old_batch_callback = joblib.parallel.BatchCompletionCallBack
    joblib.parallel.BatchCompletionCallBack = TqdmBatchCompletionCallback
    try:
        yield tqdm_object
    finally:
        joblib.parallel.BatchCompletionCallBack = old_batch_callback
        tqdm_object.close()

def fisher_exact_test(self,count_ad,feature_id=None,group1="Epi_tumor",group2="Epi_normal",log1p=True, n_threads=0):
    if type(group1)==str:
        group1=self.illumina_sample_info.cluster_midway_sets["Epi_tumor"]
    if type(group2)==str:
        group2=self.illumina_sample_info.cluster_midway_sets["Epi_normal"]
            
    group1=list(group1)
    group2=list(group2)
    
    if feature_id is None:
        feature_id_list=list(count_ad.var.index)
    if type(feature_id)==str:
        feature_id_list=[feature_id]
    
    odds_ratio_list=list()
    pvalue_list=list()
    ctable_list=list()

    
    def loop_fn(fid):
        group1_pos=count_ad[group1,fid].X.nnz
        group1_neg=len(group1)-group1_pos
        group2_pos=count_ad[group2,fid].X.nnz
        group2_neg=len(group2)-group2_pos
        ctable=np.array([group1_pos,group1_neg,group2_pos,group2_neg])
        oddsr, p=fisher_exact(ctable.reshape((2,2)),alternative="two-sided")
        return (ctable,oddsr,p)

    print("Performing Fisher exact test...")
    if n_threads>0:
        with joblib.parallel_backend(backend="loky"):
            with joblib.Parallel(n_jobs=n_threads) as parallel, tqdm_joblib(tqdm(total=len(feature_id_list))) as progress_bar: 
                with threadpool_limits(limits=1, user_api="openmp"), threadpool_limits(limits=1, user_api="blas"):
                    results=parallel(
                        joblib.delayed(loop_fn)(fid) for fid in feature_id_list
                    )

                for ctable,oddsr,p in results:
                    odds_ratio_list.append(oddsr)
                    pvalue_list.append(p)
                    ctable_list.append(ctable)

    else:
        for fid in tqdm(feature_id_list):
            ctable,oddsr,p=loop_fn(fid)
            odds_ratio_list.append(oddsr)
            pvalue_list.append(p)
            ctable_list.append(ctable)
    
    if len(pvalue_list)==1:
        odds_ratio_list=odds_ratio_list[0]
        pvalue_list=pvalue_list[0]
        ctable_list=ctable_list[0]
        
    return {
        "oddsr":odds_ratio_list,
        "pvalue":pvalue_list,
        "ctable":ctable_list
    }