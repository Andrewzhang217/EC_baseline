#!/usr/bin/env python
import os
os.environ["OMP_THREAD_LIMIT"] = "1"  # So that each process only uses 1 cpu
import argparse
from XGBoost_features import generate_input_features, aux_data
from sequences import get_reads
from concurrent.futures import ProcessPoolExecutor
from tqdm import tqdm
from concurrent.futures import as_completed
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import xgboost as xgb
import numpy as np
from overlaps import *
import time
import globals
decoder = {
    0 : 'A',
    1 : 'C',
    2 : 'G',
    3 : 'T',
    4 : ''
    #5 : ''
}

def decode(i):
    return decoder[i]

def stitch_sequence(predictions, start_idx, end_idx):
    return ''.join([decoder[pred] for pred in predictions[start_idx:end_idx]])

'''
Construct sequences from model predictions.
'''
def construct_sequences(prediction_array, aux_data_list:aux_data):
    starts = list(np.cumsum([aux.num_pos for aux in aux_data_list]))
    starts = [0] + starts[:-1]
    sequences = [SeqRecord(Seq(stitch_sequence(prediction_array, start, start + aux.num_pos)),
                    id=globals.reads[aux.tname].id,
                    name=globals.reads[aux.tname].name,
                    description=globals.reads[aux.tname].description)
                    for start, aux in zip(starts, aux_data_list)]
    
    return sequences
    
    

'''
Function given to pool. 
Generate features for the given dict of target reads and overlaps and use the given model to do inference.
Returns the corrected reads for the batch.
'''
def generate_features_and_infer(overlaps: dict[str, list[Overlap]]):  

    '''
    batch_input_features is (n, 5) np array
    '''
    batch_input_features, aux_data_list = generate_input_features(globals.reads, overlaps)
    batch_input_features = xgb.DMatrix(data=batch_input_features)
    
    # 1-d contiguous array of predictions
    pred = globals.xgb_model.predict(batch_input_features)
    pred = np.asarray([np.argmax(logit) for logit in pred])

    # construct sequences from predictions
    sequences = construct_sequences(pred, aux_data_list)
     
    return sequences
    

'''
'main' routine - given various input paths, do inference and write corrected 
reads to files.
'''
def do_inference(model_path:str, paf_path:str, reads_path:str, output_path:str, num_proc:int):
    
    overlaps = load_overlaps(paf_path)
    globals.parse_input_reads(reads_path)
    globals.load_xgb_model(model_path)


    overlap_keys = list(overlaps)
    overlap_list = overlaps.items()
    
    futures_corrected_reads = []

    time5 = time.time()
    
    with ProcessPoolExecutor(max_workers=num_proc) as executor:
        # step = int(len(overlap_map) / workers)
        step = 10
        # managed_reads = manager.dict(reads)
        # for i in range(0, 1000, step):
        
        for i in range(0, len(overlap_list), step):
            end = min(i + step, len(overlap_list))
            curr_dict = {k: overlaps[k] for k in overlap_keys[i:end]}

            f = executor.submit(generate_features_and_infer, curr_dict)
            futures_corrected_reads.append(f)
            
        with open(output_path, 'w') as f_out:
            for result in tqdm(as_completed(futures_corrected_reads)):
                # seq_lst.extend(result.result())
                SeqIO.write(result.result(), f_out, 'fasta')

            # Write uncorrected reads
            for r_id, record in reads.items():
                if r_id not in overlaps:
                    f_out.write(f'>{r_id} {record.description}\n')
                    f_out.write(f'{record.seq}\n')
        
    time6 = time.time()
    print(f'Time taken for correction : {time6 - time5}')

    

def parse_args():
    parser = argparse.ArgumentParser(
        description='Do reads correction with XGBoost.')
    parser.add_argument('-i',
                        '--input',
                        type=str,
                        help='path of the input FASTA/Q file.',
                        required=True)
    parser.add_argument('-p',
                        '--paf',
                        type=str,
                        help='path of the input PAF file.',
                        required=True)
    parser.add_argument('-m',
                        '--model',
                        type=str,
                        help='path of the model file.',
                        required=True)
    parser.add_argument('-o',
                        '--out',
                        type=str,
                        help='path of the output file.',
                        required=True)
    parser.add_argument('-t',
                        '--threads',
                        type=int,
                        help='number of threads(processes actually',
                        default=1)
    args = parser.parse_args()
    return args

def main():
    args = parse_args()
    do_inference(args.model, args.paf, args.input, args.out, args.threads)
                
if __name__ == '__main__':
    main()