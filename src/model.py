import numpy as np
import os

import tensorflow.compat.v1 as tf
tf.disable_v2_behavior()

import tensorflow
physical_devices = tensorflow.config.experimental.list_physical_devices('GPU')
assert len(physical_devices) > 0, "Not enough GPU hardware devices available"
tensorflow.config.experimental.set_memory_growth(physical_devices[0], True)

os.environ["TF_ENABLE_GPU_GARBAGE_COLLECTION"] = "false"

from tensorflow.keras.layers import dot, Input
from tensorflow.compat.v1.keras.models import Sequential
from tensorflow.compat.v1.keras.layers import Dense, Reshape
from tensorflow.compat.v1.keras.layers import CuDNNLSTM
from tensorflow.compat.v1.keras import initializers as initializers, regularizers
from tensorflow.compat.v1.keras.layers import Bidirectional, TimeDistributed, BatchNormalization
from tensorflow.compat.v1.keras.layers import Embedding
from tensorflow.compat.v1.keras import optimizers
from tensorflow.compat.v1.keras.models import Model
from atnlayer import AttentionWithContext

def skipgramModel(vocab_size, emb_size):

    input_target = Input((1,))
    input_context = Input((1,))

    embedding = Embedding(vocab_size, emb_size, input_length=1, name='embedding')
    target = embedding(input_target)
    target = Reshape((emb_size, 1))(target)
    context = embedding(input_context)
    context = Reshape((emb_size, 1))(context)

    dot_product = dot([target, context], normalize=False, axes=1)
    dot_product = Reshape((1,))(dot_product)

    output = Dense(1, activation='sigmoid')(dot_product)

    model = Model(inputs=[input_target, input_context], outputs=output)
    model.compile(loss='binary_crossentropy', optimizer='adam')
    print(model.summary())

    return model

def classifier(embeddings, wreg, atnreg):

    model = Sequential()
    model.add(Embedding(embeddings.shape[0], embeddings.shape[1], input_length=None, weights=[embeddings], trainable=False, dtype=np.float32))
    model.add(Bidirectional(CuDNNLSTM(300, return_sequences=True, kernel_regularizer=regularizers.l2(wreg), recurrent_regularizer=regularizers.l2(wreg), bias_regularizer=regularizers.l2(wreg), dtype=np.float32)))
    model.add(BatchNormalization())
    model.add(Bidirectional(CuDNNLSTM(300, return_sequences=True, kernel_regularizer=regularizers.l2(wreg), recurrent_regularizer=regularizers.l2(wreg), bias_regularizer=regularizers.l2(wreg), dtype=np.float32)))
    model.add(BatchNormalization())
    model.add(TimeDistributed(Dense(100, activation='relu', kernel_regularizer=regularizers.l2(wreg), bias_regularizer=regularizers.l2(wreg), dtype=np.float32)))
    model.add(BatchNormalization())
    model.add(AttentionWithContext(W_regularizer=regularizers.l2(atnreg), u_regularizer=regularizers.l2(atnreg), b_regularizer=regularizers.l2(atnreg), name='atnlayer'))
    model.add(Dense(1, activation='sigmoid', kernel_regularizer=regularizers.l2(wreg), bias_regularizer=regularizers.l2(wreg), dtype=np.float32))

    model.compile(loss='binary_crossentropy', optimizer='adam', metrics=['accuracy'])
    print(model.summary())

    return model