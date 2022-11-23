import time
from sequences import get_reads
import xgboost as xgb

def parse_input_reads(reads_path:str):
    global reads
    time_start = time.time()
    reads = get_reads(reads_path)
    time_end = time.time()
    print(f"Time taken for parsing input reads: {time_end-time_start}")
    
def parse_truth_reads(reads_path:str):
    global truth_reads
    time_start = time.time()
    truth_reads = get_reads(reads_path)
    time_end = time.time()
    print(f"Time taken for parsing truth reads: {time_end-time_start}")

def load_xgb_model(model_path:str):
    global xgb_model 
    xgb_model = xgb.Booster(model_file=model_path)