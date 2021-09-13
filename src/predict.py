import mio as io
import os
import model as md
import utils as mu
import argparse
import pandas as pd
import numpy as np

def predict(embeddings, weights, sequences, labels, groups, batch_size, wreg, atnreg, max_seq_len, out_path):

    batches = mu.createBatchFromDF(sequences, labels, groups, max_seq_length=max_seq_len, min_batch_size=batch_size)
    total_batches = len(batches)

    labels_from_batch = []
    ordering = []
    for v in batches:
        df = pd.concat(v[1], axis=1).T
        labels_from_batch += df['label'].values.tolist()
        ordering.append(df)

    labels_from_batch = np.array(labels_from_batch)
    ordering = pd.concat(ordering).reset_index(drop=True)

    model = md.classifier(embeddings, wreg, atnreg)
    model.set_weights(weights)

    generator = mu.get_generator(batches, test=True)
    predicted_values = model.predict_generator(generator=generator, steps=total_batches, verbose=1)

    prediction_metric_file = os.path.join(out_path, 'prediction_metric.txt')
    mu.computeMetrics(labels_from_batch, predicted_values, prediction_metric_file)

    prediction_file = os.path.join(out_path, 'predictions.csv')
    ordering = pd.concat([ordering, pd.Series(predicted_values.flatten(), name='preds')], axis=1)
    ordering.to_csv(prediction_file)


if __name__ == '__main__':

    parser = argparse.ArgumentParser()

    parser.add_argument('-e', '--embeddings', dest='embedding_file', type=str, required=True)
    parser.add_argument('-m', '--modelFile', dest='model_file', type=str, required=True)
    parser.add_argument('-t', '--sequenceFile', dest='sequences_file', type=str, required=True)
    parser.add_argument('-f', '--fileType', dest='file_type', type=str, default='pkl')
    parser.add_argument('-l', '--labelFile', dest='label_file', type=str, required=True)
    parser.add_argument('-T', '--groupFile', dest='group_file', type=str, required=True)
    parser.add_argument('-b', '--batchSize', dest='batch_size', type=int, default=25)
    parser.add_argument('-w', '--weightReg', dest='weight_regularizer', type=float, default=1e-5)
    parser.add_argument('-a', '--attentionReg', dest='attention_regularizer', type=float, default=1e-6)
    parser.add_argument('-o', '--outputPath', dest='output_path', type=str, default='./')

    args = parser.parse_args()

    embedding_file = args.embedding_file
    model_file = args.model_file
    sequences_file = args.sequences_file
    file_type = args.file_type
    label_file = args.label_file
    group_file = args.group_file
    batch_size = args.batch_size
    wreg = args.weight_regularizer
    atnreg = args.attention_regularizer
    output_path = args.output_path

    if file_type == 'pkl':
        embeddings = io.frompickle(embedding_file)[0]
    elif file_type == 'txt':
        embeddings = np.genfromtxt(embedding_file)

    weights = io.frompickle(model_file)
    sequences, max_seq_len = io.loadBin(sequences_file)
    labels = pd.read_csv(label_file, index_col=0)
    groups = io.frompickle(group_file)

    predict(embeddings, weights, sequences, labels, groups, batch_size, wreg, atnreg, max_seq_len, output_path)
