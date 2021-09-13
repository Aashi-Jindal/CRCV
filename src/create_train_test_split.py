import numpy as np
import pandas as pd
import sys
import os
import utils as mu
import mio as io
import argparse
import time
import warnings
import sklearn.model_selection as skms


def createTrainTestSplit(E_svf, E_sf, E_if, C_svf, C_sf, C_if, train_ratio, max_seq_len, min_gene_freq, random_state, aa_dict, switch_dict, remove_synonymous, remove_label_noise, out_path):

    #Inner function
    def sampleGeneIndexes(variants_info, gene_list):

        indexes = []
        genes = variants_info['Gene'].astype(np.str)

        for g in gene_list:
            inds = np.where(genes == g)[0]
            if len(inds) > 0:
                indexes += inds.tolist()

        return indexes

    print(f'[Info::createTrainTestSplit] Loading files')
    ExAC_spliced_variants = pd.read_csv(E_svf, index_col=0)
    ExAC_indexes = np.genfromtxt(E_if, dtype=np.int)
    ExAC_spliced_variants = ExAC_spliced_variants.iloc[ExAC_indexes].reset_index(drop=True)
    ExAC_sequences, ExAC_max_sequence_length = io.loadBin(E_sf)

    COSMIC_spliced_variants = pd.read_csv(C_svf, index_col=0)
    COSMIC_indexes = np.genfromtxt(C_if, dtype=np.int)
    COSMIC_spliced_variants = COSMIC_spliced_variants.iloc[COSMIC_indexes].reset_index(drop=True)
    COSMIC_sequences, COSMIC_max_sequence_length = io.loadBin(C_sf)

    ExAC_uq_gene, ExAC_uq_gene_count = np.unique(ExAC_spliced_variants['Gene'].astype(np.str), return_counts=True)
    COSMIC_uq_gene, COSMIC_uq_gene_count = np.unique(COSMIC_spliced_variants['Gene'].astype(np.str), return_counts=True)

    ExAC_gene_keep = ExAC_uq_gene[np.where(ExAC_uq_gene_count > min_gene_freq)[0]]
    COSMIC_gene_keep = COSMIC_uq_gene[np.where(COSMIC_uq_gene_count > min_gene_freq)[0]]

    unique_genes = np.union1d(ExAC_gene_keep, COSMIC_gene_keep)

    if train_ratio > 0 and train_ratio < 1:
        train_count = int(len(unique_genes)*train_ratio)
        np.random.seed(random_state)
        unique_genes_perm = np.random.permutation(unique_genes)

        train_genes_list = [unique_genes_perm[:train_count]]
        test_genes_list = [unique_genes_perm[train_count:]]
    elif train_ratio > 1:
        if isinstance(train_ratio, float):
            warnings.warn(f'train_split is more than 1 and float. Only the integer part of the number will be taken.')
            train_ratio = int(train_ratio)
            if train_ratio == 1:
                raise ValueError('train_split is 1. Which is not possible. Valid ranges are (0, 1) and {2,..}')

            else:
                kfold = skms.KFold(n_splits=train_ratio)
                train_genes_list = []
                test_genes_list = []
                for train_index, test_index in kfold.split(unique_genes):
                    train_genes_list.append(unique_genes[train_index])
                    test_genes_list.append(unique_genes[test_index])


    fold_id = -1
    for train_genes, test_genes in zip(train_genes_list, test_genes_list):
        fold_id += 1
        ExAC_train_indexes = sampleGeneIndexes(ExAC_spliced_variants, train_genes)
        COSMIC_train_indexes = sampleGeneIndexes(COSMIC_spliced_variants, train_genes)

        ExAC_test_indexes = sampleGeneIndexes(ExAC_spliced_variants, test_genes)
        COSMIC_test_indexes = sampleGeneIndexes(COSMIC_spliced_variants, test_genes)

        ExAC_train_variants = ExAC_spliced_variants.iloc[ExAC_train_indexes].reset_index(drop=True)
        COSMIC_train_variants = COSMIC_spliced_variants.iloc[COSMIC_train_indexes].reset_index(drop=True)

        ExAC_test_variants = ExAC_spliced_variants.iloc[ExAC_test_indexes].reset_index(drop=True)
        COSMIC_test_variants = COSMIC_spliced_variants.iloc[COSMIC_test_indexes].reset_index(drop=True)

        train_data = []
        train_label = []
        train_lengths = []
        for i in ExAC_train_indexes:
            s = ExAC_sequences[i]
            train_data.append(s)
            train_label.append(0)
            train_lengths.append(len(s))

        for i in COSMIC_train_indexes:
            s = COSMIC_sequences[i]
            train_data.append(s)
            train_label.append(1)
            train_lengths.append(len(s))

        test_data = []
        test_label = []
        test_lengths = []
        for i in ExAC_test_indexes:
            s = ExAC_sequences[i]
            test_data.append(s)
            test_label.append(0)
            test_lengths.append(len(s))

        for i in COSMIC_test_indexes:
            s = COSMIC_sequences[i]
            test_data.append(s)
            test_label.append(1)
            test_lengths.append(len(s))

        train_info = pd.concat([ExAC_train_variants, COSMIC_train_variants], axis=0, sort=True).reset_index(drop=True)
        test_info = pd.concat([ExAC_test_variants, COSMIC_test_variants], axis=0, sort=True).reset_index(drop=True)

        train_info = pd.concat([train_info, pd.Series(train_label, name='label'), pd.Series(train_lengths, name='lengths')], axis=1)
        test_info = pd.concat([test_info, pd.Series(test_label, name='label'), pd.Series(test_lengths, name='lengths')], axis=1)

        if remove_synonymous:
            print(f'[Info::createTrainTestSplit] Fold ID-{fold_id} Remove synonymous from train data')
            start = time.time()
            train_data, train_info = mu.removeSynonymousMutations(train_data, train_info, aa_dict, switch_dict)
            time_taken = time.time() - start
            print(f'[Info::createTrainTestSplit] Fold ID-{fold_id} Time taken {mu.processTime(time_taken)}')

            print(f'[Info::createTrainTestSplit] Fold ID-{fold_id} Remove synonymous from test data')
            start = time.time()
            test_data, test_info = mu.removeSynonymousMutations(test_data, test_info, aa_dict, switch_dict)
            time_taken = time.time() - start
            print(f'[Info::createTrainTestSplit] Fold ID-{fold_id} Time taken {mu.processTime(time_taken)}')


        if remove_label_noise:
            print(f'[Info::createTrainTestSplit] Fold ID-{fold_id} Remove label noise from train data')
            start = time.time()
            train_data, train_info = mu.removeLabelNoise(train_data, train_info)
            time_taken = time.time() - start
            print(f'[Info::createTrainTestSplit] Fold ID-{fold_id} Time taken {mu.processTime(time_taken)}')


            print(f'[Info::createTrainTestSplit] Fold ID-{fold_id} Remove label noise from test data')
            start = time.time()
            test_data, test_info = mu.removeLabelNoise(test_data, test_info)
            time_taken = time.time() - start
            print(f'[Info::createTrainTestSplit] Fold ID-{fold_id} Time taken {mu.processTime(time_taken)}')


        if max_seq_len > 0:
            print(f'[Info::createTrainTestSplit] Fold ID-{fold_id} Remove sequences longer than {max_seq_len} from train data')
            start = time.time()
            train_data, train_info = mu.removeLognerSequence(train_data, train_info, max_seq_len)
            time_taken = time.time() - start
            print(f'[Info::createTrainTestSplit] Fold ID-{fold_id} Time taken {mu.processTime(time_taken)}')


            print(f'[Info::createTrainTestSplit] Fold ID-{fold_id} Remove sequences longer than {max_seq_len} from test data')
            start = time.time()
            test_data, test_info = mu.removeLognerSequence(test_data, test_info, max_seq_len)
            time_taken = time.time() - start
            print(f'[Info::createTrainTestSplit] Fold ID-{fold_id} Time taken {mu.processTime(time_taken)}')


        print(f'[Info::createTrainTestSplit] Fold ID-{fold_id} Creating Train groups')
        start = time.time()
        train_group = mu.createGroups(train_info)
        time_taken = time.time() - start
        print(f'[Info::createTrainTestSplit] Fold ID-{fold_id} Time taken {mu.processTime(time_taken)}')

        print(f'[Info::createTrainTestSplit] Fold ID-{fold_id} Creating Test groups')
        start = time.time()
        test_group = mu.createGroups(test_info)
        time_taken = time.time() - start
        print(f'[Info::createTrainTestSplit] Fold ID-{fold_id} Time taken {mu.processTime(time_taken)}')


        train_data_file = os.path.join(out_path, f'train_data_{fold_id}.bin')
        train_info_file = os.path.join(out_path, f'train_info_{fold_id}.csv')
        train_group_file = os.path.join(out_path, f'train_group_{fold_id}.pkl')

        test_data_file = os.path.join(out_path, f'test_data_{fold_id}.bin')
        test_info_file = os.path.join(out_path, f'test_info_{fold_id}.csv')
        test_group_file = os.path.join(out_path, f'test_group_{fold_id}.pkl')

        io.saveBin(train_data_file, train_data)
        train_info.to_csv(train_info_file)
        io.topickle(train_group_file, train_group)

        io.saveBin(test_data_file, test_data)
        test_info.to_csv(test_info_file)
        io.topickle(test_group_file, test_group)


        train_genes = train_info['Gene'].unique().shape[0]
        train_mrna = train_info['#name'].unique().shape[0]
        train_samples = train_info.shape[0]
        train_imbalance = train_info['label'].mean()

        test_genes = test_info['Gene'].unique().shape[0]
        test_mrna = test_info['#name'].unique().shape[0]
        test_samples = test_info.shape[0]
        test_imbalance = test_info['label'].mean()

        print(f'[Info::createTrainTestSplit] Fold ID-{fold_id} Train Gene {train_genes}\tTrain splice variants {train_mrna}\tTrain samples {train_samples}\tTrain imbalance {train_imbalance:0.4f}')
        print(f'[Info::createTrainTestSplit] Fold ID-{fold_id} Test Gene {test_genes}\tTest splice variants {test_mrna}\tTest samples {test_samples}\t test imbalance {test_imbalance:0.4f}')


if __name__ == '__main__':

    parser = argparse.ArgumentParser()

    parser.add_argument('-e', '--ExACSplicedVariantsFile', dest='ExAC_spliced_variant_file', type=str, required=True)
    parser.add_argument('-E', '--ExACSequencesFile', dest='ExAC_sequences_file', type=str, required=True)
    parser.add_argument('-i', '--ExACIndexes', dest='ExAC_index_file', type=str, required=True)
    parser.add_argument('-c', '--COSMICSplicedVariantFile', dest='COSMIC_spliced_variant_file', type=str, required=True)
    parser.add_argument('-C', '--COSMICSequencesFile', dest='COSMIC_sequences_file', type=str, required=True)
    parser.add_argument('-I', '--COSMICIndexFile', dest='COSMIC_index_file', type=str, required=True)
    parser.add_argument('-f', '--trainSplit', dest='train_split_fraction', type=float, default=0.7)
    parser.add_argument('-l', '--maxSequenceLength', dest='max_sequence_length', type=int, default=1500)
    parser.add_argument('-m', '--minGeneFrequencey', dest='min_gene_frequencey', type=int, default=200)
    parser.add_argument('-r', '--randomState', dest='random_state', type=int, default=0)
    parser.add_argument('-a', '--aminoAcidMapFile', dest='amino_acid_map_file', type=str, required=True)
    parser.add_argument('-s', '--switchFile', dest='switch_file', type=str, required=True)
    parser.add_argument('-R', '--removeSynonymous', dest='remove_synonymous', action='store_true')
    parser.add_argument('-u', '--removeLabelNoise', dest='remove_label_noise', action='store_true')
    parser.add_argument('-o', '--output_path', dest='output_path', type=str, default='./')

    args = parser.parse_args()

    ExAC_spliced_variant_file = args.ExAC_spliced_variant_file
    ExAC_sequences_file = args.ExAC_sequences_file
    ExAC_index_file = args.ExAC_index_file
    COSMIC_spliced_variant_file = args.COSMIC_spliced_variant_file
    COSMIC_sequences_file = args.COSMIC_sequences_file
    COSMIC_index_file = args.COSMIC_index_file
    train_split_fraction = args.train_split_fraction
    max_sequence_length = args.max_sequence_length
    min_gene_frequencey = args.min_gene_frequencey
    random_state = args.random_state
    amino_acid_map_file = args.amino_acid_map_file
    switch_file = args.switch_file
    remove_synonymous = args.remove_synonymous
    remove_label_noise = args.remove_label_noise
    output_path = args.output_path

    aa_dict = io.readAAFileAsDict(amino_acid_map_file)
    switch_dict = io.readSwitchFileAsDict(switch_file, inv=True)

    createTrainTestSplit(ExAC_spliced_variant_file, ExAC_sequences_file, ExAC_index_file, COSMIC_spliced_variant_file, COSMIC_sequences_file, COSMIC_index_file, train_split_fraction, max_sequence_length, min_gene_frequencey, random_state, aa_dict, switch_dict, remove_synonymous, remove_label_noise, output_path)
