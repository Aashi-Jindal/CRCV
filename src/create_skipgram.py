import mio as io
import argparse
import os
import skipgram
import random
from collections import defaultdict
import utils as mu
import pandas as pd
import numpy as np
import time

def createSkipgram(ss_file, sc_file, ws, ns, seed, out_path, save_tuples):

    random.seed(seed)

    sequences, _ = io.loadBin(ss_file)
    counts =  pd.read_csv(sc_file, index_col=0)
    frequencies = (counts/np.sum(counts)).values
    probs = np.zeros_like(frequencies)

    for i in range(frequencies.shape[0]):
        if frequencies[i] > 0:
            probs[i] = (np.sqrt(frequencies[i]/ 0.001) + 1) * (0.001/frequencies[i])
            if probs[i] > 1:
                probs[i] = 1

    context_count = defaultdict(int)

    start = time.time()
    couples = []
    labels = []
    for i in range(len(sequences)):
        print(f'[Progress::createSkipgram] --- This may take some time {100*(i+1)/len(sequences):0.2f}%', end='\r')

        _couples, _labels = skipgram.skipgrams(sequences[i], len(frequencies), window_size=ws, sampling_table=probs, shuffle=False, negative_samples=ns)

        for cps, _lbs in zip(_couples, _labels):
            context_count[tuple(cps + [_lbs])] += 1

        del _couples
        del _labels
    time_taken = time.time() - start
    print(f'\n[Info::createSkipgram] Time taken {mu.processTime(time_taken)}')

    skipgram_count_file = os.path.join(out_path, f'tupleCount_ws{ws}_ns{ns}.csv')
    io.saveSkipgramCounts(skipgram_count_file, context_count)

    if save_tuples:
        skipgram_bin_path = os.path.join(out_path, 'idx')
        mu.mkdir(skipgram_bin_path)
        skipgram_file_stub = os.path.join(skipgram_bin_path, f'skipgram_ws{ws}_ns{ns}')
        io.saveSkipgram(skipgram_file_stub, skipgram_count_file, seed, len(frequencies))


if __name__ == '__main__':

    parser = argparse.ArgumentParser()

    parser.add_argument('-w', '--windowSize', dest='window_size', type=int, default=3)
    parser.add_argument('-n', '--negativeSampling', dest='negative_sampling', type=float, default=0.2)
    parser.add_argument('-c', '--switchCountFile', dest='switch_count_file', type=str, required=True)
    parser.add_argument('-r', '--randomState', dest='random_state', type=int, default=0)
    parser.add_argument('-S', '--switchSequenceFile', dest='switch_sequence_file', type=str, required=True)
    parser.add_argument('-o', '--outputPath', dest='output_path', type=str, default='./')
    parser.add_argument('-s', '--saveTuples', dest='save_tuples', action='store_true')

    args = parser.parse_args()

    window_size = args.window_size
    negative_sampling = args.negative_sampling
    switch_count_file = args.switch_count_file
    switch_sequence_file = args.switch_sequence_file
    random_state = args.random_state
    output_path = args.output_path
    save_tuples = args.save_tuples

    start = time.time()
    createSkipgram(switch_sequence_file, switch_count_file, window_size, negative_sampling, random_state, output_path, save_tuples)
    time_taken = time.time() - start
    print(f'[Info::__main__] Total Time taken {mu.processTime(time_taken)}')
