import numpy as np
import pandas as pd
import pickle as pkl
import struct
import idx_format as idx
import time
import utils as mu

def readSwitchFileAsDict(fname, inv=False):

    content = np.genfromtxt(fname, dtype=np.str, delimiter=',')
    if inv:
        switch_dict = {int(b):a.split('_')[-1] for a, b in content}
    else:
        switch_dict = {a.split('_')[-1]:int(b) for a, b in content}
    return switch_dict

def readAAFileAsDict(fname):

    content = np.genfromtxt(fname, dtype=np.str, delimiter=',')
    amino_map = {a:b for a,b in content}
    return amino_map

def topickle(fname, data):
    with open(fname, 'wb') as fid:
        pkl.dump(data, fid)

def frompickle(fname):
    with open(fname, 'rb') as fid:
        data = pkl.load(fid)
    return data

def totxtWithIndex(fname, counts, header=None, delimiter=','):
    with open(fname, 'w') as fid:

        if header is not None:
            header_str = delimiter.join(header)
            fid.write(f'{header_str}\n')

        for i, j in enumerate(counts):
            fid.write(f'{i}{delimiter}{j}\n')

def saveBin(fname, data):

    with open(fname, 'wb') as fid:

        total_count = len(data)
        fid.write(struct.pack('Q', total_count))

        for i in range(total_count):
            print(f'Saving metainfo {fname}: {100*(i+1)/total_count:0.2f}%', end='\r')
            fid.write(struct.pack('H', len(data[i])))
        print('')

        for i in range(total_count):
            print(f'Saving sequences {fname}: {100*(i+1)/total_count:0.2f}%', end='\r')
            data_points = data[i]
            for j in range(len(data_points)):
                fid.write(struct.pack('H', data_points[j]))
        print('')


def loadBin(fname):

    print(f'Loading file {fname}')
    with open(fname, 'rb') as fid:
        buf = fid.read()

    n_samples = np.frombuffer(buf, dtype=np.uint64, count=1, offset=0)[0]
    n_elements = np.frombuffer(buf, dtype=np.uint16, count=n_samples, offset=8)

    max_length = np.max(n_elements)

    values = []
    count = 0
    for e, i in enumerate(n_elements):
        print(f'Arranging data {100*(e+1)/n_samples:0.2f}%', end='\r')
        v = np.frombuffer(buf, dtype=np.uint16, count=i, offset=int(8+2*(n_samples+count)))
        values.append(v)
        count += i

    print('')
    return values, max_length

def saveIDX(fname, data):

    save = idx.save_idx()
    save.save(data, fname)

def loadIDX(fname):

    load = idx.load_idx()
    return load.load(fname)


def saveSkipgramCounts(fname, count_dict):

    print(f'[Info::saveSkipgramCounts] Saving counts to file {fname}')
    with open(fname, 'w') as fid:
        fid.write('word,context,label,count\n')
        for key in count_dict.keys():
            fid.write(f'{key[0]},{key[1]},{key[2]},{count_dict[key]}\n')

def saveSkipgram(skipgram_bin_fname, skipgram_count_fname, seed, vocabulary_size):

    tuple_metadata = pd.read_csv(skipgram_count_fname)
    total_tuples = tuple_metadata['count'].sum()
    metadata_size = tuple_metadata.shape[0]

    print(f'[Info::saveSkipgram] Total unique tuples/Total tuples: {tuple_metadata.shape[0]}/{total_tuples}')

    tuple_array = np.zeros((total_tuples, 3), dtype=np.int16)

    curr_index = 0
    for i, d in enumerate(tuple_metadata.values):
        print(f'[Progress::saveSkipgram] Generating full list --- {100*(i+1)/metadata_size:0.2f}%', end='\r')

        c = d[-1]

        tuple_array[curr_index:curr_index+c, 0] = d[0]
        tuple_array[curr_index:curr_index+c, 1] = d[1]
        tuple_array[curr_index:curr_index+c, 2] = d[2]

        curr_index += c

    print(f'\n[Info::saveSkipgram] Shuffling array')
    start = time.time()
    np.random.shuffle(tuple_array)
    time_taken = time.time() - start
    print(f'[Info::saveSkipgram] Time taken {mu.processTime(time_taken)}')

    total_samples_per_file = (1024*1024*1024*4)//4 # (Equivalent to 4GB - number of bytes in 4 GB/ size of one sample in bytes (tuple of 2 each with 16 bits == 4 bytes))
    min_samples_in_file = (1024*1024*1024)//4 # Equivalent to 1GB

    if total_tuples < total_samples_per_file:
        total_samples_per_file = total_tuples
        min_samples_in_file = total_tuples

    ranges = list(range(0, total_tuples, total_samples_per_file))

    if total_tuples - ranges[-1] < min_samples_in_file:
        ranges[-1] = total_tuples
    else:
        ranges.append(total_tuples)

    groups = [[ranges[i-1], ranges[i]] for i in range(1, len(ranges))]

    fnames = [[vocabulary_size, -1]]
    count = 0
    print(f'[Info::saveSkipgram] Saving skipgrams to')
    for v in groups:

        count += 1
        print(f'[Progress::saveSkipgram] Saving group {count}/{len(groups)}-({v[0]},{v[1]})')

        d = tuple_array[v[0]:v[1], :2]
        l = tuple_array[v[0]:v[1], 2].astype(np.int8)

        data_file_path = skipgram_bin_fname + f'_data_{v[0]}_{v[1]}.idx'
        label_file_path = skipgram_bin_fname + f'_label_{v[0]}_{v[1]}.idx'

        fnames.append([data_file_path, label_file_path])
        print(f'\t{data_file_path}\t{label_file_path}')

        saveIDX(data_file_path, d)
        saveIDX(label_file_path, l)

    filename_file = skipgram_bin_fname + f'_fnames.csv'
    np.savetxt(filename_file, fnames, fmt='%s', delimiter=',')
    print(f'[Info::saveSkipgram] Saved file info at path {filename_file}')

def loadSkipgram(skipgram_bin_filename_fname):

    file_info = np.genfromtxt(skipgram_bin_filename_fname, dtype=np.str, delimiter=',')
    vocabulary_size = int(file_info[0, 0])
    file_info = file_info[1:]

    datalist = []
    labellist = []

    count = 0
    for df, lf in file_info:
        count += 1
        print(f'[Progress::loadSkipgram] Loading group {count}/{file_info.shape[0]}')
        datalist.append(loadIDX(df))
        labellist.append(loadIDX(lf))

    data = np.concatenate(datalist, axis=0)
    label = np.concatenate(labellist, axis=0)

    return data, label, vocabulary_size