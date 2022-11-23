import numpy as np
from typing import *
from XGBoost_features import aux_data
from os.path import join, isdir
from os import mkdir
from XGBoost_features import INPUTS_TYPE
from XGBoost_features import LABELS_TYPE
ARRAYS_FILE = 'ARRAYS.npz'
INDEX_FILE = 'INDEX.tsv'

INPUT_SUFFIX = 'INPUTS'
LABELS_SUFFIX = 'LABELS'

AUX_SUFFIX = 'AUX'

def aux_data_to_array(auxiliary_data:List[aux_data])->np.ndarray:
    num_pos_with_ins_counts = [[len(aux.ins_counts)] + aux.ins_counts for aux in auxiliary_data]
    flattened = [item for sublist in num_pos_with_ins_counts for item in sublist]    
    return np.array(flattened)

def array_to_aux_data(arr:np.ndarray, key:str, index_dict:Dict)->List[aux_data]:
    i = 0
    idx = 0
    while i < arr.shape[0]:
        segment_len = arr[i]
        segment_start = i + 1
        ins_counts = arr[segment_start:segment_start+segment_len]
        num_pos = sum(ins_counts) + segment_len
        tname = index_dict[(key, idx)]
        idx += 1
        i += segment_start + segment_len
        yield aux_data(tname, ins_counts, num_pos)


class DataStorage:
    
    def __write_index__(self, auxiliary_data:List[aux_data], key:str)->None:
        for i, aux in enumerate(auxiliary_data):
            self.index_file.write(f'{aux.tname}\t{key}\t{i}\n')
    def __read_index__(self)->Dict:
        index_dict = {}
        for line in self.index_file:
            line = line.strip().split('\t')
            tname = line[0]
            key = line[1]
            idx = int(line[2])
            index_dict[tname] = (key, idx)
            index_dict[(key, idx)] = tname
        return index_dict
    def __init__(self, storage_dir:str, mode:str):
        self.dir = storage_dir

        if mode == 'r':
            
            self.arrays_file = open(join(storage_dir, ARRAYS_FILE), mode +'b')
            self.index_file = open(join(storage_dir, INDEX_FILE), mode)
            
            self.arrays_dict = np.load(self.arrays_file)
            self.index_dict = self.__read_index__()
        elif mode == 'w':
            if isdir(storage_dir):
                raise ValueError(f'{storage_dir} already exists!')
            mkdir(storage_dir)
            self.arrays_file = open(join(storage_dir, ARRAYS_FILE), mode+'b')
            self.index_file = open(join(storage_dir, INDEX_FILE), mode)
        else:
            raise ValueError(f'Invalid mode: {mode}')
    
    def write_batch(self, key:str, input_features:np.ndarray, auxliary_data:List[aux_data],
                    labels:np.ndarray=None)->None:
        inputs_identifier = key + INPUT_SUFFIX
        labels_identifier = key + LABELS_SUFFIX
        aux_identifier = key + AUX_SUFFIX
        
        save_dict = {}
        
        # inputs
        save_dict[inputs_identifier] = input_features
        
        # labels
        if labels:
            save_dict[labels_identifier] = labels
            
        # ins counts plus num pos
        save_dict[aux_identifier] = aux_data_to_array(auxliary_data)
        
        # name and index data
        self.__write_index__(auxliary_data, key)

        np.savez(self.arrays_file, **save_dict)
    def read_batch(self, key:str) -> Tuple[np.ndarray, np.ndarray, List[aux_data]]:
        inputs_identifier = key + INPUT_SUFFIX
        labels_identifier = key + LABELS_SUFFIX
        aux_identifier = key + AUX_SUFFIX
        labels = None
        if labels_identifier in self.arrays_dict:
             labels = self.arrays_dict[labels_identifier]
        
        aux_array = self.arrays_dict[aux_identifier]
        
        return self.arrays_dict[inputs_identifier], labels, list(array_to_aux_data(aux_array, key, self.index_dict)) 
    
    def read_target_data(self, tname:str) -> Tuple[np.ndarray, np.ndarray, aux_data]:
        key, idx = self.index_dict[tname]
        
        inputs_identifier = key + INPUT_SUFFIX
        labels_identifier = key + LABELS_SUFFIX
        aux_identifier = key + AUX_SUFFIX
        labels = None
        if labels_identifier in self.arrays_dict:
             labels = self.arrays_dict[labels_identifier]
        
        inputs_array = self.arrays_dict[inputs_identifier]
        
        aux_array = self.arrays_dict[aux_identifier]
        aux_list = list(array_to_aux_data(aux_array, key, self.index_dict))
        aux = aux_list[idx]
        start = sum([aux.num_pos for aux in aux_list[:idx]])
        return inputs_array[start:start+aux.num_pos], labels[start:start+aux.num_pos], aux
    
    def close(self):
        self.arrays_file.close()
        self.index_file.close()
    
   