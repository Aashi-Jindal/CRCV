import mio as io
import os
import model as md
import utils as mu
import argparse

def train(info_file, emb_size, disp_int, start_epoch, end_epoch, batch_size, mid_path, end_path):

    data, label, vocabulary_size = io.loadSkipgram(info_file)
    model = md.skipgramModel(vocabulary_size, emb_size)

    for i in range(start_epoch, end_epoch+1):

        print(f'Epoch {i}/{end_epoch}')
        if i > 1 and i == start_epoch:
            weight_file = os.path.join(end_path, f'weights_emb{emb_size}_itr_{start_epoch - 1}.pkl')
            weights = io.frompickle(weight_file)
            model.set_weights(weights)

        model.fit([data[:, 0], data[:, 1]], label, epochs=1, verbose=1, batch_size=batch_size, shuffle=True)

        if i%disp_int == 0:
            output_weights_file = os.path.join(mid_path, f'weights_emb{emb_size}_itr_{i}.pkl')
            weights = model.get_weights()
            io.topickle(output_weights_file, weights)

    output_weights_file = os.path.join(end_path, f'weights_emb{emb_size}_itr_{end_epoch}.pkl')
    weights_list = model.get_weights()
    io.topickle(output_weights_file, weights_list)


if __name__ == '__main__':

    parser = argparse.ArgumentParser()

    parser.add_argument('-t', '--tupleCountInfoFile', dest='tuple_count_info_file', type=str, required=True)
    parser.add_argument('-e', '--embeddingSize', dest='embedding_size', type=int, default=300)
    parser.add_argument('-p', '--peekInterval', dest='peek_interval', type=int, default=10)
    parser.add_argument('-S', '--startEpoch', dest='start_epoch', type=int, default=0)
    parser.add_argument('-E', '--endEpoch', dest='end_epoch', type=int, default=200)
    parser.add_argument('-b', '--batchSize', dest='batch_size', type=int, default=8192)
    parser.add_argument('-o', '--outputPath', dest='output_path', type=str, default='./')

    args = parser.parse_args()

    tuple_count_info_file = args.tuple_count_info_file
    embedding_size = args.embedding_size
    peek_interval = args.peek_interval
    start_epoch = args.start_epoch
    end_epoch = args.end_epoch
    batch_size = args.batch_size
    output_path = args.output_path

    midway_folder = os.path.join(output_path, 'midway')
    end_folder = os.path.join(output_path, 'finish')

    mu.mkdir(midway_folder)
    mu.mkdir(end_folder)

    train(tuple_count_info_file, embedding_size, peek_interval, start_epoch, end_epoch, batch_size, midway_folder, end_folder)
