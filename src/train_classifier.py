import mio as io
import os
import model as md
import utils as mu
import argparse
import pandas as pd
import numpy as np
import warnings

def train_combined_f(embeddings, seq, label, group, start_e, end_e, peek, batch_size, wreg, areg, max_seq_length, out_path_base):

    out_path = os.path.join(out_path_base, f'combined')
    mu.mkdir(out_path)

    batches = mu.createBatch(seq, label, group, max_seq_length=max_seq_length, min_batch_size=batch_size)
    total_batches = len(batches)
    labels_from_batch = []
    for v in batches:
        labels_from_batch += v[1].tolist()

    labels_from_batch = np.array(labels_from_batch)
    model = md.classifier(embeddings, wreg, areg)

    for i in range(start_e, end_e+1):
        if i > 1 and i == start_e:
            weight_path_file = os.path.join(out_path, f'model_epoch{i-1}', 'weights.pkl')
            weights = io.loadpickle(weight_path_file)
            model.set_weights(weights)

        generator_for_training = mu.get_generator_shuffle(batches)

        print(f'Train epochs {i}/{end_e}')
        model.fit_generator(generator=generator_for_training, steps_per_epoch=total_batches, epochs=1, verbose=1, shuffle=True)

        if (i%peek == 0) or (i%valid_ind == 0):
            weights_list = model.get_weights()

            train_generator = mu.get_generator(batches, test=True)
            prediction_train = model.predict_generator(generator=train_generator, steps=total_batches, verbose=1)
            train_preds = pd.DataFrame(np.c_[labels_from_batch, prediction_train], columns=['label','preds'])

            weight_folder = os.path.join(out_path, f'model_epoch{i}')
            mu.mkdir(weight_folder)

            weights_file = os.path.join(weight_folder, 'weights.pkl')
            io.topickle(weights_file, weights_list)

            train_prediction_file = os.path.join(weight_folder, 'train_preds.csv')
            train_preds.to_csv(train_prediction_file)

            train_metric_file = os.path.join(weight_folder, 'train_metric.txt')
            mu.computeMetrics(labels_from_batch, train_preds['preds'], train_metric_file)


def train(embeddings, tr_seq_list, tr_l_list, va_seq_list, va_l_list, tg_list, vg_list, strt_e, end_e, peek, valid_ind, batch_size, wreg, areg, max_seq_length_list, out_path_base):

    fold_id = -1
    n_folds = len(tr_seq_list)
    for tr_seq, tr_l, va_seq, va_l, tg, vg, max_seq_length in zip(tr_seq_list, tr_l_list, va_seq_list, va_l_list, tg_list, vg_list, max_seq_length_list):

        fold_id += 1

        out_path = os.path.join(out_path_base, f'fold{fold_id}')
        mu.mkdir(out_path)

        train_batches = mu.createBatch(tr_seq, tr_l, tg, max_seq_length=max_seq_length, min_batch_size=batch_size)
        validation_batches = mu.createBatch(va_seq, va_l, vg, max_seq_length=max_seq_length, min_batch_size=batch_size)

        total_train_batches = len(train_batches)
        total_validation_batches = len(validation_batches)

        train_labels_from_batch = []
        for v in train_batches:
            train_labels_from_batch += v[1].tolist()

        validation_labels_from_batch = []
        for v in validation_batches:
            validation_labels_from_batch += v[1].tolist()

        train_labels_from_batch = np.array(train_labels_from_batch)
        validation_labels_from_batch = np.array(validation_labels_from_batch)

        model = md.classifier(embeddings, wreg, areg)

        for i in range(strt_e, end_e+1):
            if i > 1 and i == strt_e:
                weight_path_file = os.path.join(out_path, f'model_epoch{i-1}', 'weights.pkl')
                weights = io.loadpickle(weight_path_file)
                model.set_weights(weights)

            generator_for_training = mu.get_generator_shuffle(train_batches)

            print(f'[{fold_id+1}/{n_folds}]Train epochs {i}/{end_e}')
            model.fit_generator(generator=generator_for_training, steps_per_epoch=total_train_batches, epochs=1, verbose=1, shuffle=True)

            if (i%peek == 0) or (i%valid_ind == 0):
                weights_list = model.get_weights()

                train_generator = mu.get_generator(train_batches, test=True)
                prediction_train = model.predict_generator(generator=train_generator, steps=total_train_batches, verbose=1)
                train_preds = pd.DataFrame(np.c_[train_labels_from_batch, prediction_train], columns=['label','preds'])

                weight_folder = os.path.join(out_path, f'model_epoch{i}')
                mu.mkdir(weight_folder)

                weights_file = os.path.join(weight_folder, 'weights.pkl')
                io.topickle(weights_file, weights_list)

                train_prediction_file = os.path.join(weight_folder, 'train_preds.csv')
                train_preds.to_csv(train_prediction_file)

                train_metric_file = os.path.join(weight_folder, 'train_metric.txt')
                mu.computeMetrics(train_labels_from_batch, train_preds['preds'], train_metric_file)

                if i%valid_ind == 0:

                    validation_generator = mu.get_generator(validation_batches, test=True)
                    prediction_valid = model.predict_generator(generator=validation_generator, steps=total_validation_batches, verbose=1)
                    valid_preds = pd.DataFrame(np.c_[validation_labels_from_batch, prediction_valid], columns=['label', 'preds'])

                    valid_prediction_file = os.path.join(weight_folder, 'valid_preds.csv')
                    valid_preds.to_csv(valid_prediction_file)

                    valid_metric_file = os.path.join(weight_folder, 'valid_metric.txt')
                    mu.computeMetrics(validation_labels_from_batch, valid_preds['preds'], valid_metric_file)


if __name__ == '__main__':

    parser = argparse.ArgumentParser()

    parser.add_argument('-e', '--embeddings', dest='embedding_file', type=str, required=True)
    parser.add_argument('-t', '--splitFileLocations', dest='train_sequences_file', type=str, required=True)
    parser.add_argument('-f', '--trainSplit', dest='train_split_fraction', type=float, default=0.7)
    parser.add_argument('-F', '--fileType', dest='file_type', type=str, default='pkl')
    parser.add_argument('-S', '--startEpoch', dest='start_epoch', type=int, default=1)
    parser.add_argument('-E', '--endEpoch', dest='end_epoch', type=int, default=200)
    parser.add_argument('-p', '--peekInterval', dest='peek_interval', type=int, default=1)
    parser.add_argument('-P', '--valdationInterval', dest='validation_interval', type=int, default=1)
    parser.add_argument('-b', '--batchSize', dest='batch_size', type=int, default=25)
    parser.add_argument('-w', '--weightReg', dest='weight_regularizer', type=float, default=1e-5)
    parser.add_argument('-a', '--attentionReg', dest='attention_regularizer', type=float, default=1e-6)
    parser.add_argument('-c', '--trainCombined', dest='train_combined', action='store_true')
    parser.add_argument('-o', '--outputPath', dest='output_path', type=str, default='./')

    args = parser.parse_args()

    embedding_file = args.embedding_file
    split_file_location = args.train_sequences_file
    train_split = args.train_split_fraction
    file_type = args.file_type
    start_epoch = args.start_epoch
    end_epoch = args.end_epoch
    peek_interval = args.peek_interval
    validation_interval = args.validation_interval
    batch_size = args.batch_size
    weight_regularizer = args.weight_regularizer
    attention_regularizer = args.attention_regularizer
    train_combined = args.train_combined
    output_path = args.output_path


    if not train_combined:
        if train_split > 0 and train_split < 1:
            train_split = [0]
        elif isinstance(train_split, float):
            warnings.warn(f'train_split is more than 1 and float. Only the integer part of the number will be taken.')
            train_split = int(train_split)
            if train_split == 1:
                raise ValueError('train_split is 1. Which is not possible. Valid ranges are (0, 1) and {2,..}')

            train_split = list(range(train_split))

    if file_type == 'pkl':
        embeddings = io.frompickle(embedding_file)[0]
    elif file_type == 'txt':
        embeddings = np.genfromtxt(embedding_file)

    #train_sequences, max_train_sequence = io.loadBin(train_sequences_file)

    if train_combined:

        train_sequences, max_train_sequence = io.loadBin(os.path.join(split_file_location, f'train_data_0.bin'))
        test_sequences, max_test_sequence = io.loadBin(os.path.join(split_file_location, f'test_data_0.bin'))

        sequences = train_sequences + test_sequences
        max_seq_length = max(max_train_sequence, max_test_sequence)

        train_info = pd.read_csv(os.path.join(split_file_location, f'train_info_0.csv'))
        test_info = pd.read_csv(os.path.join(split_file_location, f'test_info_0.csv'))
        info = pd.concat((train_info, test_info), axis=0).reset_index(drop=True)

        labels = info['label'].values
        groups = mu.createGroups(info)

        train_combined_f(embeddings, sequences, labels, groups, start_epoch, end_epoch, peek_interval, batch_size, weight_regularizer, attention_regularizer, max_seq_length, output_path)

    else:
        sequenes = [io.loadBin(os.path.join(split_file_location, f'train_data_{i}.bin')) for i in train_split]
        train_sequences, max_train_sequence = list(map(list, list(zip(*sequenes))))
        train_label = [pd.read_csv(os.path.join(split_file_location, f'train_info_{i}.csv'))['label'].values for i in train_split]

        sequenes = [io.loadBin(os.path.join(split_file_location, f'test_data_{i}.bin')) for i in train_split]
        validation_sequences, max_validation_sequence = list(map(list, list(zip(*sequenes))))
        validation_label = [pd.read_csv(os.path.join(split_file_location, f'test_info_{i}.csv'))['label'].values for i in train_split]

        train_groups = [io.frompickle(os.path.join(split_file_location, f'train_group_{i}.pkl')) for i in train_split]
        validation_groups = [io.frompickle(os.path.join(split_file_location, f'test_group_{i}.pkl')) for i in train_split]

        #max_seq_length = max(max_train_sequence, max_validation_sequence)
        max_seq_length = [max(i, j) for i, j in zip(max_train_sequence, max_validation_sequence)]
        train(embeddings, train_sequences, train_label, validation_sequences, validation_label, train_groups, validation_groups, start_epoch, end_epoch, peek_interval, validation_interval, batch_size, weight_regularizer, attention_regularizer, max_seq_length, output_path)
