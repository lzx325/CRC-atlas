import os
import sys
from os.path import join,basename,dirname,splitext
from pathlib import Path
from collections import OrderedDict
import types
import numpy as np
import scipy
from scipy import io
import scipy.sparse as sp
import matplotlib.pyplot as plt
import matplotlib as mpl
import pandas as pd
import seaborn as sns
from pprint import pprint
class CNVPlotter(object):
    def __init__(self,rm_grouped,expr_df,meta_data_df):
        self.rm_grouped=rm_grouped
        self.expr_df=expr_df
        self.meta_data_df=meta_data_df

    @staticmethod
    def sort_by_custom_order(aux_data_df,sort_by_dict):
        assert isinstance(sort_by_dict,OrderedDict)
        for col_name,mapping in reversed(sort_by_dict.items()):
            if mapping is None:
                aux_data_df=aux_data_df.sort_values(col_name,kind="stable")
            elif type(mapping)==list:
                mapping_dict=dict(zip(mapping,range(len(mapping))))
                mapping_fn=lambda s: s.apply(lambda x: mapping_dict[x])
                aux_data_df=aux_data_df.sort_values(col_name,key=mapping_fn,kind="stable")
            elif type(mapping)==dict:
                mapping_dict=mapping
                mapping_fn=lambda s: s.apply(lambda x: mapping_dict[x])
                aux_data_df=aux_data_df.sort_values(col_name,key=mapping_fn,kind="stable")
            elif type(mapping)==types.FunctionType:
                aux_data_df=aux_data_df.sort_values(col_name,key=mapping,kind="stable")

        return aux_data_df
    @staticmethod
    def filter_df(aux_data_df,filter_dict):
        for col_name,subset in filter_dict.items():
            aux_data_df=aux_data_df[np.isin(aux_data_df[col_name],subset)]
        return aux_data_df
    def gather_data(self,aux_features,sort_by_dict,filter_dict,supplementary_df):
        all_unique_features=set(aux_features).union(sort_by_dict.keys()).union(filter_dict.keys())
        all_columns=list()
        for feature in all_unique_features:
            if feature in self.expr_df.columns:
                all_columns.append(self.expr_df[feature])
            elif feature in self.meta_data_df.columns:
                all_columns.append(self.meta_data_df[feature])
            elif type(supplementary_df)==pd.DataFrame and feature in supplementary_df.columns:
                all_columns.append(supplementary_df[feature])
            else:
                assert False, feature+" not found"
        aux_data_df=pd.concat(all_columns,axis=1,join="inner")
        aux_data_df=self.filter_df(aux_data_df,filter_dict)
        aux_data_df=self.sort_by_custom_order(aux_data_df,sort_by_dict)
        return self.rm_grouped,aux_data_df

    @staticmethod
    def draw_categorical_legends(ax,labels,cmap):
        from matplotlib.lines import Line2D
        if cmap==str:
            cmap=plt.get_cmap(cmap)
        assert type(cmap)==mpl.colors.ListedColormap
        legend_elements = [
            Line2D([0], [0], color=cmap(i%(len(cmap.colors))), marker="s", lw=0, label=labels[i], markersize=15)
            for i in range(len(labels))
        ]
        ax.legend(handles=legend_elements, loc='center',bbox_to_anchor=(0.0,0.0, 1,1))
        ax.axis('off')
    def plot_cnv_heatmap(
        self,
        chr_to_plot=["chr1"],
        aux_features=["library_id"],
        sort_by_dict=OrderedDict([("library_id",None)]),
        filter_dict=OrderedDict([]),
        aux_plot_percentage=0.02,
        aux_categorical_legend_percentage=0.05,
        aux_continuous_legend_percentage=0.03,
        main_legend_percentage=0.03,
        cell_sample=None,
        categorical_cmap="Set3",
        continuous_cmap="coolwarm",
        custom_categorical_order=dict(),
        supplementary_df=None,
        vmin=None,
        center=0,
        vmax=None,
    ):
        assert len(aux_features)<=2
        aux_feature_names=list()
        for i in range(len(aux_features)):
            if type(aux_features[i])==str:
                aux_feature_names.append(aux_features[i])
            elif type(aux_features[i])==tuple:
                aux_feature_names.append(aux_features[i][0])

        rm_grouped,aux_data_df=self.gather_data(aux_feature_names,sort_by_dict,filter_dict,supplementary_df)

        if type(categorical_cmap)==str:
            categorical_cmap=plt.get_cmap(categorical_cmap)
        if type(continuous_cmap)==str:
            continuous_cmap=plt.get_cmap(continuous_cmap)

        assert type(categorical_cmap)==mpl.colors.ListedColormap
 
        aux_feature_dtypes=list()
        aux_width_ratios=list()
        for i in range(len(aux_features)):
            if type(aux_features[i])==str:
                if type(aux_data_df[aux_features[i]][0])==str:
                    aux_feature_dtype="categorical"
                    aux_width_ratios+=[aux_categorical_legend_percentage,aux_plot_percentage]
                elif issubclass(type(aux_data_df[aux_features[i]][0]),np.number):
                    aux_feature_dtype="continuous"
                    aux_width_ratios+=[aux_continuous_legend_percentage,aux_plot_percentage]
                else:
                    assert False, type(aux_data_df[aux_features[i]][0])

            elif type(aux_features[i])==tuple:
                aux_feature_dtype=aux_features[i][1]
                assert aux_feature_dtype in ["categorical","continuous"], aux_feature_dtype
                if aux_feature_dtype == "categorical":
                    aux_width_ratios+=[aux_categorical_legend_percentage,aux_plot_percentage]
                elif aux_feature_dtype == "continuous" :
                    aux_width_ratios+=[aux_continuous_legend_percentage,aux_plot_percentage]
            aux_feature_dtypes.append(aux_feature_dtype)
            
        main_width_ratios=np.array([rm_grouped[ch].shape[1] for ch in chr_to_plot])
        main_width_ratios_total=1-sum(aux_width_ratios)-main_legend_percentage
        main_width_ratios_total=main_width_ratios_total*main_width_ratios/np.sum(main_width_ratios)
        width_ratios=np.concatenate([aux_width_ratios,main_width_ratios_total,[main_legend_percentage]])
        fig,axes=plt.subplots(1,len(width_ratios),gridspec_kw={'width_ratios': width_ratios})
        fig.set_size_inches((70,20))

        if cell_sample and cell_sample<aux_data_df.shape[0]:
            np.random.seed(0)
            selected_indices=np.random.choice(aux_data_df.shape[0],cell_sample,replace=False)
            selected_indices.sort()
            selected_barcodes=aux_data_df.iloc[selected_indices].index
            aux_data_df=aux_data_df.loc[selected_barcodes]
        else:
            selected_barcodes=aux_data_df.index
            aux_data_df=aux_data_df

        # draw qualitative legends
        for i in range(len(aux_feature_names)):
            feature_dtype=aux_feature_dtypes[i]

            if feature_dtype=="categorical":
                if aux_feature_names[i] in custom_categorical_order:
                    codes=pd.Categorical(aux_data_df[aux_feature_names[i]],categories=custom_categorical_order[aux_feature_names[i]]).codes
                    labels=custom_categorical_order[aux_feature_names[i]]
                else:
                    codes=pd.Categorical(aux_data_df[aux_feature_names[i]]).codes
                    labels=aux_data_df[aux_feature_names[i]].unique()
                    labels.sort()
                self.draw_categorical_legends(axes[2*i],labels,cmap=categorical_cmap)
                # draw qualitative map
                plt.sca(axes[2*i+1])
                sns.heatmap(np.mod(codes[:,None],len(categorical_cmap.colors)),
                                cmap=categorical_cmap,xticklabels=False,yticklabels=False,cbar=False,vmin=0,vmax=len(categorical_cmap.colors))
                # return np.mod(codes[:,None],len(categorical_cmap.colors))
            elif feature_dtype=="continuous":
                plt.sca(axes[2*i+1])
                sns.heatmap(aux_data_df[aux_feature_names[i]].values[:,None],
                                cmap=continuous_cmap,xticklabels=False,yticklabels=False,cbar=True,cbar_ax=axes[2*i])
            
        
        main_plot_axis_start_index=len(aux_width_ratios)
        for i,ch in enumerate(chr_to_plot,main_plot_axis_start_index):
            print(i,ch)
            plt.sca(axes[i])
            expr_rm=rm_grouped[ch].loc[selected_barcodes].values
            expr_rm=np.where((expr_rm>-0.3)&(expr_rm<0.3),0.0,expr_rm)
            print("main heatmap",expr_rm.shape)
            if i==main_plot_axis_start_index+len(chr_to_plot)-1:
                
                sns.heatmap(expr_rm,
                cmap=continuous_cmap,xticklabels=False,yticklabels=False,cbar=True,cbar_ax=axes[i+1],vmin=vmin,vmax=vmax,center=center)
            else:
                sns.heatmap(expr_rm,
                            cmap=continuous_cmap,xticklabels=False,yticklabels=False,cbar=False,vmin=vmin,vmax=vmax,center=center)

            axes[i].set_xlabel(ch.replace("chr",""))