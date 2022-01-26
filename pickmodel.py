import scanpy as sc
import numpy as np
import sys
sys.path.insert(0, './mvTCR')

from tcr_embedding.models.model_selection import run_model_selection
import tcr_embedding.utils_training as utils
from tcr_embedding.utils_preprocessing import group_shuffle_split

import os
import argparse

'''
If you don't have them already, please install the following packages before running this code:

- optuna
- PyTorch
- comet_ml

'''

path_tcr_anno_with_gender = './data/07_tcr_annotation_with_gender_data.h5ad'
adata = sc.read(path_tcr_anno_with_gender)
randindices = np.random.choice(adata.shape[0], size = 100, replace = False)
adata = adata[randindices] # Randomly pick only 100 observations from our data

utils.fix_seeds(42)

parser = argparse.ArgumentParser()
parser.add_argument('--model', type=str, default='poe')
parser.add_argument('--split', type=int, default=0)
parser.add_argument('--gpus', type=int, default=1)
args = parser.parse_args()

# Put the .h5ad file read below in a folder named 'data' and put the folder 'data' in the 'mvTCR' folder.
# This way you won't get an error from the utils.load_data function
adata = utils.load_data('07_tcr_annotation_with_gender_data.h5ad')


random_seed = args.split
train, val = group_shuffle_split(adata, group_col='clonotype', val_split=0.20, random_seed=random_seed)
adata.obs['set'] = 'train'
adata.obs.loc[val.obs.index, 'set'] = 'val'
adata = adata[adata.obs['set'].isin(['train', 'val'])]

params_experiment = {
    'study_name': f'TCR_moe_split_{args.split}',
    'comet_workspace': None, 
    'model_name': 'moe',
    'balanced_sampling': 'clonotype',
    'metadata': ['clonotype', 'cell types', 'Gender'],
    'save_path': os.path.join(os.path.dirname(__file__), 'optuna', 
                              f'TCR_moe_split_{args.split}')
    # In the current directory, create a folder called 'optuna' so that this works
}

params_optimization = {
    'name': 'pseudo_metric',
    'prediction_labels': ['clonotype', 'cell types'],
}


timeout = (2 * 24 * 60 * 60) - 300
run_model_selection(adata, params_experiment, params_optimization, None, timeout, args.gpus)

