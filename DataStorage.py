import h5py
import numpy as np
from typing import *
from XGBoost_features import aux_data, INPUTS_TYPE, LABELS_TYPE
'''
This file aims to support facilities to store generated feature data for reuse, reference and debugging.

Requirements:
1. Supports the storage and retrieval of all input samples and their corresponding labels if they exist.
2. Supports the storage of relevant auxiliary data to trace each sample to
   a target read and the specific position on it.
3. Supports the storage of relevant auxiliary data to construct the corrected reads from the predictions.

----
4. Supports storage/selective retrieval in batches 



Design ideas:
1. Need some kind of index object in the file to find data for a specific target read
2. Each batch should contain input features, optionally contains labels, optionally contains auxiliary 
data 
3. maybe the different batches will just be named by some index number? since they each will
consist of data for many target reads, not really anything to distinguish them (the index range of target
reads in some container is possible, but it may not be a fixed order? e.g. dictionary items)

structure:
    -root group
      -index_group
      -data_group
        -batch1_data_group
          -
        -batch2_data_group
        .
        .
        .
      

'''

'''
Since auxiliary data format can be changing a lot. Put this function here
'''
def write_auxiliary_data(auxiliary_data:List[aux_data], batch_group:h5py.Group):
    # the target read names
    string_dtype = h5py.string_dtype()
    target_names_np = np.array([aux.tname for aux in auxiliary_data], dtype=string_dtype)
    batch_group.create_dataset('tnames', data=target_names_np)
    
    # the list of lists of ins counts
    int_list_dtype = h5py.vlen_dtype(np.dtype('int32'))
    
    ins_lists_np = [np.array(aux.ins_counts) for aux in auxiliary_data]
    batch_group.create_dataset('ins_counts', data=ins_lists_np, dtype=int_list_dtype)
    
    # the total number of positions for each target
    num_pos_list_np = np.array([aux.num_pos for aux in auxiliary_data])
    batch_group.create_dataset('num_pos', data=num_pos_list_np)
    
# The lists will become np.ndarray but I guess it's fine
def read_auxiliary_data(batch_group:h5py.Group)->List[aux_data]:
    target_names_np = batch_group['tnames']
    ins_lists_np = batch_group['ins_counts']
    num_pos_list_np = batch_group['num_pos']
    return [aux_data(tname.decode(), ins_counts, num_pos) for
            tname, ins_counts, num_pos in zip(target_names_np, ins_lists_np, num_pos_list_np)]

# The lists will become np.ndarray but I guess it's fine
def read_auxiliary_data_single(batch_group:h5py.Group, idx:int)->List[aux_data]:
    target_names_np = batch_group['tnames']
    ins_lists_np = batch_group['ins_counts']
    num_pos_list_np = batch_group['num_pos']
    return aux_data(target_names_np[idx].decode(), ins_lists_np[idx], num_pos_list_np[idx])
    
class DataStorage:
    # mode is either 'w' or 'r'
    def __load_index__(self) -> Dict[str, List[int]]:
        index_dict = {}
        for batch_idx in range(self.num_batches):
            batch_group = self.data_group[f'batch{batch_idx}']
            auxiliary_data = read_auxiliary_data(batch_group)
            starts = list(np.cumsum([aux.num_pos for aux in auxiliary_data]))
            starts = [0] + starts[:-1]
            for i, aux in enumerate(auxiliary_data):
                index_dict[aux.tname] = [batch_idx, i, starts[i]]
        return index_dict
    '''
    Index sequences only for reading mode
    '''
    def __init__(self, path_to_file:str, mode:str, index:bool=False):
        self.mode = mode
        self.top_group= h5py.File(path_to_file, mode)
        
        # Currently supports one-off write and no more addition of data to the same file
        if mode == 'w':
            self.data_group = self.top_group.create_group('data')
            self.num_batches = 0 # should write this to an attribute somewhere
        elif mode == 'r':
            self.data_group = self.top_group['data']
            self.num_batches = self.top_group.attrs['num_batches']
            if index:
                self.index = self.__load_index__()
        else:
            raise ValueError(f'Invalid mode: {mode}')
    
    def close(self):
        if self.mode == 'w':
            self.top_group.attrs['num_batches'] = self.num_batches
        self.top_group.close()

    ''' 
    This function should create a group somewhere at a predefined location in a predefined
    directory structure. 
    
    'w' mode required
    '''
    def write_batch(self, input_array:np.ndarray, labels_array:np.ndarray, auxiliary_data:List[aux_data]) -> None:
        batch_group = self.data_group.create_group(f'batch{self.num_batches}')
        
        if labels_array == None:
            labels_array = h5py.Empty('f')
        
        # TODO maybe do chunking or sth?
        batch_group.create_dataset('inputs', data=input_array)
        
        batch_group.create_dataset('labels', data=labels_array)
        
        write_auxiliary_data(auxiliary_data, batch_group)
        
        self.num_batches += 1
    
    def read_batch(self, batch_id:int)->Tuple[np.ndarray, np.ndarray, List[aux_data]]:
        batch_group = self.data_group[f'batch{batch_id}']
        labels_ds = batch_group['labels']
        if labels_ds.shape == None:
            labels = None
        else:
            labels = np.asarray(labels_ds, LABELS_TYPE)
        return np.asarray(batch_group['inputs'], dtype=INPUTS_TYPE), \
            labels, read_auxiliary_data(batch_group)
    # requires index
    def get_data_for_target(self, tname:str):
        batch_idx, seq_idx, start_idx = self.index[tname]
        batch_group = self.data_group[f'batch{batch_idx}']
        aux = read_auxiliary_data_single(batch_group, seq_idx)
        labels_ds = batch_group['labels']
        if labels_ds.shape == None:
            labels = None
        else:
            labels = np.asarray(batch_group['labels'][start_idx:start_idx+aux.num_pos], dtype=LABELS_TYPE)
        return np.asarray(batch_group['inputs'][start_idx:start_idx+aux.num_pos], dtype=INPUTS_TYPE),\
            labels, aux
            
    
    
    