import preprocess_steps as ps
import numpy as np
import mio as io
import utils as mu
import argparse
import os
import time

def processExacFiles(chromosome, in_filename, mrna_locs, reference_genome, switch_dict, embedding_fraction, out_path, random_state):

    start = time.time()
    filtered_vcf = ps.filterExACVCF(in_filename, chromosome) #112374
    time_taken = time.time() - start
    print(f'[Info::processExacFiles] Time taken {mu.processTime(time_taken)}')

    start = time.time()
    splicing_data = ps.createAlternateSplicingData(filtered_vcf, mrna_locs)
    time_taken = time.time() - start
    print(f'[Info::processExacFiles] Time taken {mu.processTime(time_taken)}')
    del filtered_vcf

    start = time.time()
    splicing_data_with_length = ps.getLengthInformation(splicing_data, mrna_locs)
    time_taken = time.time() - start
    print(f'[Info::processExacFiles] Time taken {mu.processTime(time_taken)}')
    del splicing_data

    start = time.time()
    splicing_out_fname = os.path.join(out_path, f'chr{chromosome}_ExAC_spliced_variants_data.csv')
    flatten_variants = ps.flattenSplicingData(splicing_data_with_length, splicing_out_fname)
    time_taken = time.time() - start
    print(f'[Info::processExacFiles] Time taken {mu.processTime(time_taken)}')
    del splicing_data_with_length

    start = time.time()
    sequence_out_fstub = os.path.join(out_path, f'chr{chromosome}_ExAC_sequence')
    indexes = ps.createSequences(flatten_variants, mrna_locs, reference_genome, sequence_out_fstub)
    time_taken = time.time() - start
    print(f'[Info::processExacFiles] Time taken {mu.processTime(time_taken)}')
    index_count = len(indexes)
    del flatten_variants

    embedding_indexes_count = int(index_count * embedding_fraction)
    classification_indexes_count = index_count - embedding_indexes_count

    print(f'[Info::processExacFiles] Embedding Fraction {embedding_fraction} | Embedding Sequence Count {embedding_indexes_count} | Classification sequence count {classification_indexes_count}')

    np.random.seed(random_state)
    permuted_indexes = np.random.permutation(index_count)
    embedding_indexes = np.sort(permuted_indexes[:embedding_indexes_count])
    classification_indexes = np.sort(permuted_indexes[embedding_indexes_count:])

    embedding_index_fname = os.path.join(out_path, f'chr{chromosome}_ExAC_embedding_index.txt')
    embedding_bin_fstub = os.path.join(out_path, f'chr{chromosome}_ExAC_embedding_{embedding_fraction}')
    np.savetxt(embedding_index_fname, indexes[embedding_indexes], fmt='%d')
    print(f'[Info::processExacFiles] Creating switch sequences for Embedding')
    start = time.time()
    ps.createSwitchSequences(sequence_out_fstub, embedding_indexes, embedding_bin_fstub, switch_dict, dump_count=True)
    time_taken = time.time() - start
    print(f'[Info::processExacFiles] Time taken {mu.processTime(time_taken)}')


    classification_index_fname = os.path.join(out_path, f'chr{chromosome}_ExAC_classificaiton_index.txt')
    classification_bin_fstub = os.path.join(out_path, f'chr{chromosome}_ExAC_classificaiton')
    print(f'[Info::processExacFiles] Creating switch sequences for classification')
    start = time.time()
    kept_indexes = ps.createSwitchSequences(sequence_out_fstub, classification_indexes, classification_bin_fstub, switch_dict, remove_duplicates=True)
    time_taken = time.time() - start
    print(f'[Info::processExacFiles] Time taken {mu.processTime(time_taken)}')
    np.savetxt(classification_index_fname, indexes[classification_indexes][kept_indexes], fmt='%d')


def processCOSMICFiles(chromosome, in_filename, mrna_locs, reference_genome, switch_dict, out_path):

    start = time.time()
    filtered_vcf = ps.loadCOSMICVCF(in_filename, chromosome)
    time_taken = time.time() - start
    print(f'[Info::processCOSMICFiles] Time taken {mu.processTime(time_taken)}')


    start = time.time()
    splicing_data = ps.createAlternateSplicingData(filtered_vcf, mrna_locs)
    time_taken = time.time() - start
    print(f'[Info::processCOSMICFiles] Time taken {mu.processTime(time_taken)}')
    del filtered_vcf

    start = time.time()
    splicing_data_with_length = ps.getLengthInformation(splicing_data, mrna_locs)
    time_taken = time.time() - start
    print(f'[Info::processCOSMICFiles] Time taken {mu.processTime(time_taken)}')
    del splicing_data

    start = time.time()
    splicing_out_fname = os.path.join(out_path, f'chr{chromosome}_COSMIC_spliced_mutation_data.csv')
    flatten_variants = ps.flattenSplicingData(splicing_data_with_length, splicing_out_fname)
    time_taken = time.time() - start
    print(f'[Info::processCOSMICFiles] Time taken {mu.processTime(time_taken)}')
    del splicing_data_with_length

    sequence_out_fstub = os.path.join(out_path, f'chr{chromosome}_COSMIC_sequence')
    start = time.time()
    indexes = ps.createSequences(flatten_variants, mrna_locs, reference_genome, sequence_out_fstub)
    time_taken = time.time() - start
    print(f'[Info::processCOSMICFiles] Time taken {mu.processTime(time_taken)}')
    index_count = len(indexes)
    del flatten_variants

    all_indexes = np.arange(index_count, dtype=np.int)

    index_fname = os.path.join(out_path, f'chr{chromosome}_COSMIC_index.txt')
    bin_fstub = os.path.join(out_path, f'chr{chromosome}_COSMIC_classificaiton')
    start = time.time()
    kept_indexes = ps.createSwitchSequences(sequence_out_fstub, all_indexes, bin_fstub, switch_dict, remove_duplicates=True)
    time_taken = time.time() - start
    print(f'[Info::processCOSMICFiles] Time taken {mu.processTime(time_taken)}')

    np.savetxt(index_fname, indexes[kept_indexes], fmt='%d')

def processDBSNPFiles(chromosome, in_filename, mrna_locs, reference_genome, switch_dict, out_path):

    start = time.time()
    filtered_vcf = ps.filterDBSNPVCF(in_filename, chromosome, out_path)
    time_taken = time.time() - start
    print(f'[Info::processDBSNPFiles] Time taken {mu.processTime(time_taken)}')

    start = time.time()
    splicing_data = ps.createAlternateSplicingData(filtered_vcf, mrna_locs)
    time_taken = time.time() - start
    print(f'[Info::processDBSNPFiles] Time taken {mu.processTime(time_taken)}')
    del filtered_vcf

    start = time.time()
    splicing_data_with_length = ps.getLengthInformation(splicing_data, mrna_locs)
    time_taken = time.time() - start
    print(f'[Info::processDBSNPFiles] Time taken {mu.processTime(time_taken)}')
    del splicing_data

    splicing_out_fname = os.path.join(out_path, f'chr{chromosome}_dbSNP_spliced_mutation_data.csv')
    start = time.time()
    flatten_variants = ps.flattenSplicingData(splicing_data_with_length, splicing_out_fname)
    time_taken = time.time() - start
    print(f'[Info::processDBSNPFiles] Time taken {mu.processTime(time_taken)}')
    del splicing_data_with_length

    sequence_out_fstub = os.path.join(out_path, f'chr{chromosome}_dbSNP_sequence')
    start = time.time()
    indexes = ps.createSequences(flatten_variants, mrna_locs, reference_genome, sequence_out_fstub)
    time_taken = time.time() - start
    print(f'[Info::processDBSNPFiles] Time taken {mu.processTime(time_taken)}')
    index_count = len(indexes)
    del flatten_variants

    all_indexes = np.arange(index_count, dtype=np.int)

    index_fname = os.path.join(out_path, f'chr{chromosome}_dbSNP_index.txt')
    bin_fstub = os.path.join(out_path, f'chr{chromosome}_dbSNP_classificaiton')
    start = time.time()
    kept_indexes = ps.createSwitchSequences(sequence_out_fstub, all_indexes, bin_fstub, switch_dict, remove_duplicates=True)
    time_taken = time.time() - start
    print(f'[Info::processDBSNPFiles] Time taken {mu.processTime(time_taken)}')

    np.savetxt(index_fname, indexes[kept_indexes], fmt='%d')

def processMETFiles(chromosome, in_filename, mrna_locs, reference_genome, switch_dict, out_path):

    start = time.time()
    filtered_vcf = ps.filterMETData(in_filename, chromosome)
    time_taken = time.time() - start
    print(f'[Info::processMETFiles] Time taken {mu.processTime(time_taken)}')

    start = time.time()
    splicing_data = ps.createAlternateSplicingData(filtered_vcf, mrna_locs)
    time_taken = time.time() - start
    print(f'[Info::processMETFiles] Time taken {mu.processTime(time_taken)}')
    del filtered_vcf

    start = time.time()
    splicing_data_with_length = ps.getLengthInformation(splicing_data, mrna_locs)
    time_taken = time.time() - start
    print(f'[Info::processMETFiles] Time taken {mu.processTime(time_taken)}')
    del splicing_data

    start = time.time()
    splicing_out_fname = os.path.join(out_path, f'chr{chromosome}_mET_spliced_mutation_data.csv')
    flatten_variants = ps.flattenSplicingData(splicing_data_with_length, splicing_out_fname)
    time_taken = time.time() - start
    print(f'[Info::processMETFiles] Time taken {mu.processTime(time_taken)}')
    del splicing_data_with_length

    sequence_out_fstub = os.path.join(out_path, f'chr{chromosome}_mET_sequence')
    start = time.time()
    indexes = ps.createSequences(flatten_variants, mrna_locs, reference_genome, sequence_out_fstub)
    time_taken = time.time() - start
    print(f'[Info::processMETFiles] Time taken {mu.processTime(time_taken)}')
    index_count = len(indexes)
    del flatten_variants

    all_indexes = np.arange(index_count, dtype=np.int)

    index_fname = os.path.join(out_path, f'chr{chromosome}_mET_index.txt')
    bin_fstub = os.path.join(out_path, f'chr{chromosome}_mET_classificaiton')
    start = time.time()
    kept_indexes = ps.createSwitchSequences(sequence_out_fstub, all_indexes, bin_fstub, switch_dict, remove_duplicates=True)
    time_taken = time.time() - start
    print(f'[Info::processMETFiles] Time taken {mu.processTime(time_taken)}')

    np.savetxt(index_fname, indexes[kept_indexes], fmt='%d')


if __name__ == '__main__':

    parser = argparse.ArgumentParser()

    parser.add_argument('-C', '--chromosome', dest='_chr', type=str, default='X')
    parser.add_argument('-e', '--ExACVCF', dest='ExAC_vcf_path', type=str, default=None)
    parser.add_argument('-c', '--COSMICVCF', dest='COSMIC_vcf_path', type=str, default=None)
    parser.add_argument('-d', '--dbSnpVCF', dest='dbSNP_vcf_path', type=str, default=None)
    parser.add_argument('-m', '--mETVCF', dest='mET_vcf_path', type=str, default=None)
    parser.add_argument('-u', '--ucscKg', dest='ucsc_kg_path', type=str, required=True)
    parser.add_argument('-U', '--ucscKgref', dest='ucsc_kg_ref_path', type=str, required=True)
    parser.add_argument('-R', '--refGenomeFasta', dest='reference_genome_path', type=str, required=True)
    parser.add_argument('-s', '--switchFile', dest='switch_file', type=str, required=True)
    parser.add_argument('-E', '--embeddingFraction', dest='embedding_fraction', type=float, default=0.4)
    parser.add_argument('-r', '--randomState', dest='random_state', type=int, default=0)
    parser.add_argument('-o', '--outputPath', dest='output_path', type=str, default='./')

    args = parser.parse_args()

    files_to_process = {}

    if args.ExAC_vcf_path is not None:
        files_to_process['ExAC'] = args.ExAC_vcf_path
    if args.COSMIC_vcf_path is not None:
        files_to_process['COSMIC'] = args.COSMIC_vcf_path
    if args.dbSNP_vcf_path is not None:
        files_to_process['dbSNP'] = args.dbSNP_vcf_path
    if args.mET_vcf_path is not None:
        files_to_process['mET'] = args.mET_vcf_path

    if len(files_to_process) == 0:
        raise ValueError(f'Define atleat one of ExAC_vcf_path(-e|--ExACVCF), COSMIC_vcf_path(-c|--COSMICVCF), dbSNP_vcf_path(-d|--dbSnpVCF), mET_vcf_path(-m|--mETVCF)')


    chromosome = args._chr
    output_path = args.output_path
    ucsc_kg_path = args.ucsc_kg_path
    ucsc_kg_ref_path = args.ucsc_kg_ref_path
    reference_genome_path = args.reference_genome_path
    switch_file = args.switch_file
    embedding_fraction = args.embedding_fraction
    random_state = args.random_state

    print(f'[Info::Preprocess] Extracting Genome')
    start = time.time()
    reference_genome = ps.splitRefGenome(reference_genome_path, chromosome)
    time_taken = time.time() - start
    print(f'[Info::Preprocess] Time taken {mu.processTime(time_taken)}')

    start = time.time()
    mrna_locations = ps.processUCSCFiles(ucsc_kg_path, ucsc_kg_ref_path, chromosome)
    time_taken = time.time() - start
    print(f'[Info::Preprocess] Time taken {mu.processTime(time_taken)}')

    switch_dict = io.readSwitchFileAsDict(switch_file)
    #Now process ucsc files

    if 'ExAC' in files_to_process.keys():
        processExacFiles(chromosome, files_to_process['ExAC'], mrna_locations, reference_genome, switch_dict, embedding_fraction, output_path, random_state)
    if 'COSMIC' in files_to_process.keys():
        processCOSMICFiles(chromosome, files_to_process['COSMIC'], mrna_locations, reference_genome, switch_dict, output_path)
    if 'dbSNP' in files_to_process.keys():
        processDBSNPFiles(chromosome, files_to_process['dbSNP'], mrna_locations, reference_genome, switch_dict, output_path)
    if 'mET' in files_to_process.keys():
        processMETFiles(chromosome, files_to_process['mET'], mrna_locations, reference_genome, switch_dict, output_path)

    #parser.add_argument('-r', '--reshape_3d', dest='should_reshape', action='store_true')



