import mio as io
import os
import utils as mu
import argparse
import pandas as pd
import numpy as np
import time

def createPredictionData(seqs, splice_variants, tr_seqs, aa_dict, switch_dict, remove_synonymous, max_sequence_length, out_path):

    if remove_synonymous:
        print(f'[Info::createPredictionData] Removing synonymous mutations')
        start = time.time()
        seqs, splice_variants = mu.removeSynonymousMutations(seqs, splice_variants, aa_dict, switch_dict)
        time_taken = time.time() - start
        print(f'[Info::createPredictionData] Time taken {mu.processTime(time_taken)}')


    if max_sequence_length > 0:
        print(f'[Info::createPredictionData] Removing sequences larget than {max_sequence_length}')
        start = time.time()
        seqs, splice_variants = mu.removeLognerSequence(seqs, splice_variants, max_sequence_length)
        time_taken = time.time() - start
        print(f'[Info::createPredictionData] Time taken {mu.processTime(time_taken)}')

    print(f'[Info::createPredictionData] Removing overlap with training data')
    start = time.time()
    seqs, splice_variants = mu.removeTrainOverlap(seqs, splice_variants, tr_seqs)
    time_taken = time.time() - start
    print(f'[Info::createPredictionData] Time taken {mu.processTime(time_taken)}')

    print(f'[Info::createPredictionData] Creating groups')
    start = time.time()
    groups = mu.createGroups(splice_variants)
    time_taken = time.time() - start
    print(f'[Info::createPredictionData] Time taken {mu.processTime(time_taken)}')

    sequence_files = os.path.join(out_path, 'data.bin')
    splice_variants_file = os.path.join(out_path, 'info.csv')
    groups_file = os.path.join(out_path, 'groups.pkl')

    io.saveBin(sequence_files, seqs)
    splice_variants.to_csv(splice_variants_file)
    io.topickle(groups_file, groups)

    print(f'[Info::createPredictionData] Total remaining samples {len(seqs)}')

if __name__ == '__main__':

    parser = argparse.ArgumentParser()

    parser.add_argument('-f', '--spliceVariantFile', dest='splice_variant_file', type=str, required=True)
    parser.add_argument('-S', '--sequenceFile', dest='sequence_file', type=str, required=True)
    parser.add_argument('-i', '--indexFile', dest='index_file', type=str, default=None)
    parser.add_argument('-b', '--baseLabel', dest='base_label', type=int, default=None)
    parser.add_argument('-t', '--trainSequenceFile', dest='train_sequence_file', type=str, required=True)
    parser.add_argument('-a', '--aminoAcidMapFile', dest='amino_acid_map_file', type=str, required=True)
    parser.add_argument('-s', '--switchFile', dest='switch_file', type=str, required=True)
    parser.add_argument('-R', '--removeSynonymous', dest='remove_synonymous', action='store_true')
    parser.add_argument('-l', '--maxSequenceLength', dest='max_sequence_length', type=int, default=1500)
    parser.add_argument('-o', '--outputPath', dest='output_path', type=str, default='./')

    args = parser.parse_args()

    splice_variant_file = args.splice_variant_file
    sequence_file = args.sequence_file
    index_file = args.index_file
    base_label = args.base_label
    train_sequence_file = args.train_sequence_file
    amino_acid_map_file = args.amino_acid_map_file
    switch_file = args.switch_file
    remove_synonymous = args.remove_synonymous
    max_sequence_length = args.max_sequence_length
    output_path = args.output_path

    splice_variants = pd.read_csv(splice_variant_file)
    sequences, _ = io.loadBin(sequence_file)

    if index_file is None:
        if splice_variants.shape[0] != len(sequences):
            raise ValueError(f'row count in splice_variant_file and sequences is not same {splice_variants.shape[0]}, {len(sequences)}')
    else:
        indexes = np.genfromtxt(index_file)
        splice_variants = splice_variants.iloc[indexes].reset_index(drop=True)

    if base_label is None:
        if 'label' not in splice_variants.columns:
            raise ValueError(f'splice_variant_file missing label column, Either add a label column in splice_variant file. or provide target label (0/1).')
    else:
        label = np.array([base_label]*len(splice_variants), dtype=np.int)
        splice_variants = pd.concat([splice_variants, pd.Series(label, name='label')], axis=1)

    #train_sequences, _ = io.loadBin(train_sequence_file)

    train_sequences, _ = io.loadBin(os.path.join(train_sequence_file, f'train_data_0.bin'))
    test_sequences, _ = io.loadBin(os.path.join(train_sequence_file, f'test_data_0.bin'))

    train_sequences = train_sequences + test_sequences

    seq_lengths = [len(s) for s in sequences]
    splice_variants = pd.concat([splice_variants, pd.Series(seq_lengths, name='lengths')], axis=1)

    aa_dict = io.readAAFileAsDict(amino_acid_map_file)
    switch_dict = io.readSwitchFileAsDict(switch_file, inv=True)

    createPredictionData(sequences, splice_variants, train_sequences, aa_dict, switch_dict, remove_synonymous, max_sequence_length, output_path)