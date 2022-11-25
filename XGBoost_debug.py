#!/usr/bin/env python

from globals import parse_input_reads, load_xgb_model
import globals
from overlaps import load_overlaps, Overlap
from typing import *
from XGBoost_features import generate_input_features, aux_data
import numpy as np
from XGBoost_infer import construct_sequences
import xgboost as xgb
from Bio.SeqRecord import SeqRecord
'''
Functions to:
1. initialize global reads. Initialize overlaps
2. load model
3. generate input features for a specific read (by id)
4. use model to get prediction and construct sequence from it 

---------------------------------------
5. Get the list of overlapping reads to a target and print an alignment
'''

def initialize(reads_path:str, paf_path:str, model_path:str)->Dict[str, List[Overlap]]:
    parse_input_reads(reads_path)
    load_xgb_model(model_path)
    return load_overlaps(paf_path)

def generate_features_for_read(read_id:str, overlaps:Dict[str, List[Overlap]])-> Tuple[np.ndarray, aux_data]:
    one_item_dict = {read_id:overlaps[read_id]}
    input_features, aux_data_list = generate_input_features(globals.reads, one_item_dict)
    return input_features, aux_data_list[0]

def predict_and_construct_sequence(input_features:np.ndarray, aux:aux_data)->SeqRecord:
    input_features = xgb.DMatrix(data=input_features)
    
    # 1-d contiguous array of predictions
    pred = globals.xgb_model.predict(input_features)
    pred = np.asarray([np.argmax(logit) for logit in pred])

    # construct sequences from predictions
    sequences = construct_sequences(pred, [aux])
    return sequences[0], pred
