#!/bin/bash

EMBEDDING_FILE='models/models_skipgram_weights_emb300_itr_200.txt'
MODEL_FILE='models/weights_model.pkl'
mET_SEQUENCE_FILE='output/prediction_data/CHRX/mET/data.bin'
EMBEDDING_FILE_TYPE='txt'
LABEL_FILE='output/prediction_data/CHRX/mET/info.csv'
GROUP_FILE='output/prediction_data/CHRX/mET/groups.pkl'
BATCH_SIZE=25
WEIGHT_REG=1e-5
ATN_REG=1e-6
OUTPUT_PATH='output/predictions/CHRX/mET'

mkdir -p ${OUTPUT_PATH}

CUDA_VISIBLE_DEVICES=0 python src/predict.py -e ${EMBEDDING_FILE} -m ${MODEL_FILE} -t ${mET_SEQUENCE_FILE} -f ${EMBEDDING_FILE_TYPE} -l ${LABEL_FILE} -T ${GROUP_FILE} -b ${BATCH_SIZE} -w ${WEIGHT_REG} -a ${ATN_REG} -o ${OUTPUT_PATH}
