import scanpy as sc
import numpy as np
import sys
sys.path.insert(0, '../mvTCR')

from tcr_embedding.models.model_selection import run_model_selection
import tcr_embedding.utils_training as utils
from tcr_embedding.utils_preprocessing import group_shuffle_split, stratified_group_shuffle_split

import os
import argparse

'''
If you don't have them already, please install the following packages before running this code:

- optuna
- PyTorch
- comet_ml

'''

random_seed = 42
utils.fix_seeds(random_seed)

parser = argparse.ArgumentParser()
parser.add_argument('--rna_weight', type=int, default=1)
parser.add_argument('--model', type=str, default='rna') #Changes here: moe/poe/rna
parser.add_argument('--split', type=int, default=0)
parser.add_argument('--gpus', type=int, default=1)
args = parser.parse_args()

# Put the .h5ad file read below in a folder named 'data' and put the folder 'data' in the 'mvTCR' folder.
# This way you won't get an error from the utils.load_data function
adata = utils.load_data('09_tcr_annotation_A_B_with_gender_data.h5ad')

random_seed = args.split

# sub, non_sub = group_shuffle_split(adata, group_col='clonotype', val_split=0.2, random_seed=random_seed)
# train, val = group_shuffle_split(adata, group_col='clonotype', val_split=0.20, random_seed=random_seed)

# Create Train-Val and Test set
train, val = stratified_group_shuffle_split(adata.obs, stratify_col='cell types', group_col='clonotype', val_split=0.20, random_seed=random_seed)
# Split Train-Val into Train and Val set
#train, val = stratified_group_shuffle_split(train_val, stratify_col='cell types', group_col='clonotype', val_split=0.25, random_seed=random_seed)
adata.obs['set'] = 'train'
#adata.obs.loc[non_sub.obs.index, 'set'] = '-'
adata.obs.loc[val.index, 'set'] = 'val'
#adata = adata[adata.obs['set'].isin(['train', 'val'])]
adata.obs['set'] = adata.obs['set'].astype('category')

params_experiment = {
    'study_name': f'TCR_rna_split_{args.split}_{args.rna_weight}',
    'comet_workspace': None, 
    'model_name': 'rna', #Changed moe to rna
    'balanced_sampling': 'clonotype',
    'metadata': ['clonotype', 'cell types', 'Gender'],
    'save_path': os.path.join(os.path.dirname(__file__), 'optuna', 
                              f'TCR_rna_split_{args.split}_{args.rna_weight}')
    # In the current directory, create a folder called 'optuna' so that this works
}

if args.model == 'rna':
    params_experiment['balanced_sampling'] = None

params_optimization = {
    'name': 'pseudo_metric',
    'prediction_labels': 
    {'clonotype': 1, 
    'cell types': args.rna_weight}
}


timeout = (2 * 24 * 60 * 60) - 300
run_model_selection(adata, params_experiment, params_optimization, None, timeout, args.gpus)


