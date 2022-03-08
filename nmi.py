import pandas as pd
import numpy as np
import scanpy as sc
import sys
import scirpy as ir
import os
os.chdir('C:/Users/Gheorghe Pascu/OneDrive - tum.de/WiSe 21-22/Computational_Methods_in_Single_Cell_Biology/T_cell_project/') # where the folder 'data' can be found
cwd = os.getcwd()
sys.path.insert(0, '../mvTCR')
import tcr_embedding.utils_training as utils
from tcr_embedding.evaluation.Clustering import predict_clustering, get_clustering_scores

adata = utils.load_data('09_tcr_annotation_A_B_with_gender_data.h5ad')
latent = sc.read('./data/11_latent_moe.h5ad')
cluster_results = []

name_label = ['cell types', 'Age']
for resolution in [0.01, 0.1, 1.0]:
    print("Results for resolution ", resolution)
    print(" ")
    cluster_params={'resolution': resolution, 'num_neighbors': 5}
    for i in range(2):
        if i == 0:
            print("Results for cell types:")
        else:
            print("Results for Age:")
        print(" ")
        labels_true_adata = adata.obs[name_label[i]].to_numpy()
        labels_true_latent = latent.obs[name_label[i]].to_numpy()
        labels_predicted_adata = predict_clustering(adata, cluster_params, visualize = False, name_label = name_label[i])
        labels_predicted_latent = predict_clustering(latent, cluster_params, visualize = False, name_label = name_label[i])
        scores_adata = get_clustering_scores(adata.X, labels_true_adata, labels_predicted_adata)
        scores_latent = get_clustering_scores(latent.X, labels_true_latent, labels_predicted_latent)
        print("Scores for adata:")
        print(" ")
        print(scores_adata)
        print("Scores for latent:")
        print(" ")
        print(scores_latent)