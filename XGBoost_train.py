#!/usr/bin/env python

import argparse
from XGBoost_features import generate_training_data
from sklearn.model_selection import train_test_split
import xgboost as xgb
import time 
train_param = {
    'eta': 0.3, 
    'max_depth': 5,  
    'objective': 'multi:softprob',  
    'num_class': 5} 

train_steps = 50  # The number of training iterations





def parse_args():
    parser = argparse.ArgumentParser(
    description='Train XGBoost model for reads correction.')
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
    parser.add_argument('--truth',
                        type=str,
                        help='path of the input ground-truth sequences.',
                        required=True)
    parser.add_argument('-o',
                        '--out',
                        type=str,
                        help='path of the output model file.',
                        required=True)
    parser.add_argument('-v',
                        '--val',
                        type=float,
                        help='proportion of samples used for validation.',
                        default=0.2)
    parser.add_argument('-t',
                        '--threads',
                        type=int,
                        help='number of threads(processes actually',
                        default=1)
    parser.add_argument('-l',
                        '--limit',
                        type=int,
                        help='number of target reads to use for training. -1 for all',
                        default=-1)
    args = parser.parse_args()
    return args

def main():
    args = parse_args()
    global train_param
    train_param['n_jobs'] = args.threads
    input_features, ground_truth = generate_training_data(args.paf, args.input, args.truth, args.threads, args.limit)
    X_train, X_val, Y_train, Y_val = train_test_split(input_features, ground_truth, test_size=args.val)
    print('Loading data...')
    DM_train = xgb.DMatrix(data=X_train,label=Y_train, nthread=args.threads)
    DM_val = xgb.DMatrix(data=X_val,label=Y_val, nthread=args.threads)
    time1 = time.time()
    print('Training...')
    model = xgb.train(train_param, DM_train, train_steps, evals=[(DM_train, 'train'), (DM_val, 'eval')], early_stopping_rounds=10)
    time2 = time.time()
    print(f'Time taken for training: {time2 - time1}')
    model.save_model(args.out)
if __name__ == '__main__':
    main()