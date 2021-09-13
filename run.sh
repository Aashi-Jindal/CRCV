#!/bin/bash

echo "Preprocessing......"

CHROMOSOME=X
ExAC_VCF_FILE=data/ExAC.r1.sites.vep.vcf.gz
COSMIC_VCF_FILE=data/CosmicCodingMuts_37_cosmic_build_89.vcf.gz
dbSNP_VCF_FILE=data/GCF_000001405.25.gz
mET_VCF_FILE=data/41586_2019_1689_MOESM10_ESM.xlsx
UCSC_KG_FILE=data/exone_locations_hg19_ucsc
UCSC_KGxREF_FILE=data/exone_locations_hg19_ref_ucsc
REFERENCE_GENOME_FILE=data/Homo_sapiens_assembly19.fasta
SWITCH_FILE=data/snp_switched_dict.csv
EMBEDDING_FRACTION=0.4
RANDOM_STATE=0
PREPROCESS_OUTPUTPATH=output/preprocess/CHR${CHROMOSOME}

mkdir -p ${PREPROCESS_OUTPUTPATH}
python src/preprocess.py -C ${CHROMOSOME} -e ${ExAC_VCF_FILE} -c ${COSMIC_VCF_FILE} -d ${dbSNP_VCF_FILE} -m ${mET_VCF_FILE} -u ${UCSC_KG_FILE} -U ${UCSC_KGxREF_FILE} -R ${REFERENCE_GENOME_FILE} -s ${SWITCH_FILE} -E ${EMBEDDING_FRACTION} -r ${RANDOM_STATE} -o ${PREPROCESS_OUTPUTPATH}

echo "Creating skipgram....."

WINDOW_SIZE=3
NEGATIVE_SAMPLING=0.2
SWITCH_COUNT_FILE=${PREPROCESS_OUTPUTPATH}/chrX_ExAC_embedding_${EMBEDDING_FRACTION}_count.txt
SWITCH_SEQUENCE_FILE=${PREPROCESS_OUTPUTPATH}/chrX_ExAC_embedding_${EMBEDDING_FRACTION}_sequences.bin
CREATED_SKIPGRAM_SAVE_PATH=output/skipgram/CHR${CHROMOSOME}

mkdir -p ${CREATED_SKIPGRAM_SAVE_PATH}

python src/create_skipgram.py -w ${WINDOW_SIZE} -n ${NEGATIVE_SAMPLING} -c ${SWITCH_COUNT_FILE} -r ${RANDOM_STATE} -S ${SWITCH_SEQUENCE_FILE} -o ${CREATED_SKIPGRAM_SAVE_PATH} -s

echo "Training Embedding....."

TUPLE_COUNT_INFO_FILE=${CREATED_SKIPGRAM_SAVE_PATH}/idx/skipgram_ws3_ns0.2_fnames.csv
EMBEDDING_SIZE=300
EMB_PEEK_INTERVAL=10
EMB_START_EPOCH=1
EMB_END_EPOCH=200
EMB_BATCH_SIZE=8192
EMB_OUTPUTPATH=output/skipgram/CHR${CHROMOSOME}/embedding_${EMBEDDING_SIZE}

mkdir -p ${EMB_OUTPUTPATH}

CUDA_VISIBLE_DEVICES=0 python src/train_skipgram.py -t ${TUPLE_COUNT_INFO_FILE} -e ${EMBEDDING_SIZE} -p ${EMB_PEEK_INTERVAL} -S ${EMB_START_EPOCH} -E ${EMB_END_EPOCH} -b ${EMB_BATCH_SIZE} -o ${EMB_OUTPUTPATH}

echo "Creating train test split....."

ExAC_SPLICED_VARIANT_FILE=${PREPROCESS_OUTPUTPATH}/chrX_ExAC_spliced_variants_data.csv
ExAC_SEQUENCE_FILE=${PREPROCESS_OUTPUTPATH}/chrX_ExAC_classificaiton_sequences.bin
ExAC_INDEX_FILE=${PREPROCESS_OUTPUTPATH}/chrX_ExAC_classificaiton_index.txt
COSMIC_SPLICED_VARIANT_FILE=${PREPROCESS_OUTPUTPATH}/chrX_COSMIC_spliced_mutation_data.csv
COSMIC_SEQUENCE_FILE=${PREPROCESS_OUTPUTPATH}/chrX_COSMIC_classificaiton_sequences.bin
COSMIC_INDEX_FILE=${PREPROCESS_OUTPUTPATH}/chrX_COSMIC_index.txt
TRAIN_SPLIT_FRACTION=4
MAX_SEQUENCE_LENGTH=1500
MIN_GENE_FREQUENCEY=200
AMINO_ACID_MAP_FILE=data/aminoMap.csv
TRAIN_TEST_SPLIT_OUTPUT_FILE=output/train_test_split/CHR${CHROMOSOME}

mkdir -p ${TRAIN_TEST_SPLIT_OUTPUT_FILE}

python src/create_train_test_split.py -e ${ExAC_SPLICED_VARIANT_FILE} -E ${ExAC_SEQUENCE_FILE} -i ${ExAC_INDEX_FILE} -c ${COSMIC_SPLICED_VARIANT_FILE} -C ${COSMIC_SEQUENCE_FILE} -I ${COSMIC_INDEX_FILE} -f ${TRAIN_SPLIT_FRACTION} -l ${MAX_SEQUENCE_LENGTH} -m ${MIN_GENE_FREQUENCEY} -r ${RANDOM_STATE} -a ${AMINO_ACID_MAP_FILE} -s ${SWITCH_FILE} -o ${TRAIN_TEST_SPLIT_OUTPUT_FILE} -R -u


echo "Training CV model....."

EMBEDDING_FILE_PATH=${EMB_OUTPUTPATH}/finish/weights_emb${EMBEDDING_SIZE}_itr_${EMB_END_EPOCH}.pkl
SPLIT_FILE_LOCATION=${TRAIN_TEST_SPLIT_OUTPUT_FILE}
TRAINING_EMBEDDING_FILE_TYPE='pkl'
TRAIN_START_EPOCH=1
TRAIN_END_EPOCH=200
TRAIN_PEEK_INTERVAL=1
TRAIN_VALIDATION_INTERVAL=1
TRAIN_BATCH_SIZE=25
WEIGHT_REG=1e-5
ATN_REG=1e-6
TRAIN_OUTPUT_PATH=output/train/CHR${CHROMOSOME}

mkdir -p ${TRAIN_OUTPUT_PATH}

CUDA_VISIBLE_DEVICES=0 python src/train_classifier.py -e ${EMBEDDING_FILE_PATH} -t ${SPLIT_FILE_LOCATION} -f {TRAIN_SPLIT_FRACTION} -F ${TRAINING_EMBEDDING_FILE_TYPE} -S ${TRAIN_START_EPOCH} -E ${TRAIN_END_EPOCH} -p ${TRAIN_PEEK_INTERVAL} -p ${TRAIN_VALIDATION_INTERVAL} -b ${TRAIN_BATCH_SIZE} -w ${WEIGHT_REG} -a ${ATN_REG} -o ${TRAIN_OUTPUT_PATH}

echo "Training Combined model...."

CUDA_VISIBLE_DEVICES=0 python src/train_classifier.py -e ${EMBEDDING_FILE_PATH} -t ${SPLIT_FILE_LOCATION} -f {TRAIN_SPLIT_FRACTION} -F ${TRAINING_EMBEDDING_FILE_TYPE} -S ${TRAIN_START_EPOCH} -E ${TRAIN_END_EPOCH} -p ${TRAIN_PEEK_INTERVAL} -p ${TRAIN_VALIDATION_INTERVAL} -b ${TRAIN_BATCH_SIZE} -w ${WEIGHT_REG} -a ${ATN_REG} -o ${TRAIN_OUTPUT_PATH} -c

echo "Creating prediction data......"

mET_SPLICED_VARIANT_FILE=${PREPROCESS_OUTPUTPATH}/chrX_mET_spliced_mutation_data.csv
mET_SEQUENCE_FILE=${PREPROCESS_OUTPUTPATH}/chrX_mET_classificaiton_sequences.bin
mET_INDEX_FILE=${PREPROCESS_OUTPUTPATH}/chrX_mET_index.txt
mET_BASE_LABEL=1
mET_OUTPUT_PATH=output/prediction_data/CHR${CHROMOSOME}/mET/

mkdir -p ${mET_OUTPUT_PATH}

python src/create_prediction_data.py -f ${mET_SPLICED_VARIANT_FILE} -S ${mET_SEQUENCE_FILE} -i ${mET_INDEX_FILE} -b ${mET_BASE_LABEL} -t ${SPLIT_FILE_LOCATION} -a ${AMINO_ACID_MAP_FILE} -s ${SWITCH_FILE} -l ${MAX_SEQUENCE_LENGTH} -o ${mET_OUTPUT_PATH} -R

dbSNP_SPLICED_VARIANT_FILE=${PREPROCESS_OUTPUTPATH}/chrX_dbSNP_spliced_mutation_data.csv
dbSNP_SEQUENCE_FILE=${PREPROCESS_OUTPUTPATH}/chrX_dbSNP_classificaiton_sequences.bin
dbSNP_INDEX_FILE=${PREPROCESS_OUTPUTPATH}/chrX_dbSNP_index.txt
dbSNP_BASE_LABEL=0
dbSNP_OUTPUT_PATH=output/prediction_data/CHR${CHROMOSOME}/dbSNP/

mkdir -p ${dbSNP_OUTPUT_PATH}

python src/create_prediction_data.py -f ${dbSNP_SPLICED_VARIANT_FILE} -S ${dbSNP_SEQUENCE_FILE} -i ${dbSNP_INDEX_FILE} -b ${dbSNP_BASE_LABEL} -t ${SPLIT_FILE_LOCATION} -a ${AMINO_ACID_MAP_FILE} -s ${SWITCH_FILE} -l ${MAX_SEQUENCE_LENGTH} -o ${dbSNP_OUTPUT_PATH} -R

echo "Predict...."

MODEL_FILE=${TRAIN_OUTPUT_PATH}/combined/model_epoch200/weights.pkl
mET_SEQUENCE_FILE=${mET_OUTPUT_PATH}/data.bin
mET_LABEL_FILE=${mET_OUTPUT_PATH}/info.csv
mET_GROUP_FILE=${mET_OUTPUT_PATH}/groups.pkl
EMBEDDING_FILE_TYPE=pkl
PREDICT_BATCH_SIZE=25
PREDICT_OUTPUT_PATH=output/predictions/CHR${CHROMOSOME}/mET

mkdir -p ${PREDICT_OUTPUT_PATH}

CUDA_VISIBLE_DEVICES=0 python src/predict.py -e ${EMBEDDING_FILE_PATH} -m ${MODEL_FILE} -t ${mET_SEQUENCE_FILE} -f ${EMBEDDING_FILE_TYPE} -l ${mET_LABEL_FILE} -T ${mET_GROUP_FILE} -b ${PREDICT_BATCH_SIZE} -w ${WEIGHT_REG} -a ${ATN_REG} -o ${PREDICT_OUTPUT_PATH}

