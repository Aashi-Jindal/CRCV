# CRCV

## Contents
[Introduction](#introduction) <br />
[Prerequisites](#prerequisites) <br />
[Resources](#resources) <br />
[Usage](#usage) <br />
-[Preprocess](#preprocess) <br />
-[Create skipgram](#cskipgram) <br />
-[Train skipgram](#tskipgram) <br />
-[Train-Test split](#trtesplit) <br />
-[Train deep network](#train) <br />
-[Create prediction data](#create) <br />
-[Prediction](#prediction) <br />
-[Predict with pre-trained model](#pretrain) <br />

<a name="introduction"></a>
## Introduction

Tested on UBUNTU 18.04 LTS.

All results in manuscript have been generated using `python 3.0`.

<a name="prerequisites"></a>
## Prerequisites

<h4> Required python modules </h4>

```python
    numpy >= 1.17.4
    pandas >= 0.25.3
    pickle >= 4.0
    tensorflow >= 2.0.0
    keras >= 2.2.4
    Bio >= 1.73
```

<a name="resources"></a>
## Resources

Reference genome vhg19/Grch37 is used.

<h4>a. ExAC VCF (release 1)</h4>
<a href="https://console.cloud.google.com/storage/browser/gnomad-public/legacy/exacv1_downloads/release1">ExAC browser</a> 

Alternatively, it can also be downloaded from <a href="https://gnomad.broadinstitute.org/downloads">gnomAD</a>

<h4>b. COSMIC VCF (v89)</h4> 
<a href="https://cancer.sanger.ac.uk/cosmic/download">COMSIC VCF files (coding mutations)</a>

<h4>c. UCSC knownGene table </h4>
UCSC table browser is available at: <a href="http://genome.ucsc.edu/cgi-bin/hgTables">UCSC Table</a>

![](images/exon.png?raw=true)

<h4>d. UCSC kgXref table </h4>

![](images/ref.png?raw=true)

<h4>e. dbSNP VCF </h4>
<a href="ftp://ftp.ncbi.nih.gov/snp/latest_release/VCF%7D">dbSNP</a>

<h4>f. mET VCF (Supplementary Table 3)</h4>
<a href ="https://www.nature.com/articles/s41586-019-1689-y#MOESM1">mET </a>

<h4>g. Driver genes information </h4>

<a href="https://www.oncokb.org/cancerGenes">Precision Oncology Knowledge Base (OncoKB)</a>

<a href="https://www.intogen.org/download">Interactive Onco Genomics (IntOGen)</a>

<a href ="https://www.cancergenomeinterpreter.org/mutations">Cancer Genome Interpreter (CGI)</a>

Curated driver genes from above databases is provided in `output/GENEINFO.xls` file. 

<a name="usage"></a>
## Usage

<a name="preprocess"></a>
<h4>1. Preprocess </h4>

```bash
python src/preprocess.py 
usage: preprocess.py [-h] [-C _CHR] [-e EXAC_VCF_PATH] [-c COSMIC_VCF_PATH]
                     [-d DBSNP_VCF_PATH] [-m MET_VCF_PATH] -u UCSC_KG_PATH -U
                     UCSC_KG_REF_PATH -R REFERENCE_GENOME_PATH -s SWITCH_FILE
                     [-E EMBEDDING_FRACTION] [-r RANDOM_STATE]
                     [-o OUTPUT_PATH]
```

|Parameter | Required or Optional| Datatype | Default Value | Description |
| -----:| -----:| -----:|-----:|-----:|
|-C|Optional|`str`|`X`|Chromosome|
|-e|Optional|`str`|None|Path of ExAC VCF|
|-c|Optional|`str`|None|Path of COSMIC VCF|
|-d|Optional|`str`|None|Path of dbSNP VCF|
|-m|Optional|`str`|None|Path of mET VCF|
|-u|Required|`str`|None|Path of knownGene table from UCSC|
|-U|Required|`str`|None|Path of kgXref table from UCSC|
|-R|Required|`str`|None|Path of hg19 reference genome|
|-s|Required|`str`|None|Path of switch file|
|-E|Optional|`float`|0.4|Percentage of ExAC variants for embedding|
|-r|Optional|`int`|0|Random state|
|-o|Optional|`str`|Current Directory|Output path|

Note: It is mandatory to provide value for at least one of the parameters `-e`, `-c`, `-d`, `-m`.

The overall summary of preprocessing pipeline is illustrated as follows:

![](images/preprocessing.png?raw=true)

Example:

```bash
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

```

To process __ExAC__ file:

```bash
python src/preprocess.py -C ${CHROMOSOME} -e ${ExAC_VCF_FILE} -u ${UCSC_KG_FILE} -U ${UCSC_KGxREF_FILE} -R ${REFERENCE_GENOME_FILE} -s ${SWITCH_FILE} -E ${EMBEDDING_FRACTION} -r ${RANDOM_STATE} -o ${PREPROCESS_OUTPUTPATH}
```

Expected output:

```
[Info::Preprocess] Extracting Genome
[Info::Preprocess] Time taken 16.71s
[Progress::processUCSCFiles] --- 100.00%
[Info::Preprocess] Time taken 3.09s
[Message::filterExACVCF] Reading file data/ExAC.r1.sites.vep.vcf.gz [This step may take few minutes]
[Message::filterExACVCF] Filtering ...
[Info::processExacFiles] Time taken 5.85m
[Progress::createAlternateSplicingData] This step may take some time --- 100.00%
[Info::processExacFiles] Time taken 4.13m
[Progress::getLengthInformation] This step may take some time --- 100.00%
[Info::processExacFiles] Time taken 3.32m
[Progress::flattenSplicingData] --- 100.00%
[Info::processExacFiles] Time taken 37.33s
[Progress::createSequences] This step may take some time -- 100.00%
[Info::processExacFiles] Time taken 4.00m
[Info::processExacFiles] Embedding Fraction 0.4 | Embedding Sequence Count 114040 | Classification sequence count 171062
[Info::processExacFiles] Creating switch sequences for Embedding
[Progress::createSwitchSequences] --- 100.00%
Saving metainfo output/preprocess/chrX_ExAC_embedding_0.4_sequences.bin: 100.00%
Saving sequences output/preprocess/chrX_ExAC_embedding_0.4_sequences.bin: 100.00%
[Info::processExacFiles] Time taken 1.65m
[Info::processExacFiles] Creating switch sequences for classification
[Progress::createSwitchSequences] --- 100.00%
[Info::removeDuplicates] Total Seqences 171062
[Progress::removeDuplicates] -- 100.00%
[Info::removeDuplicates] Remaining Sequences 149455. Reduction 87.37%
Saving metainfo output/preprocess/CHRX/chrX_ExAC_classificaiton_sequences.bin: 100.00%
Saving sequences output/preprocess/CHRX/chrX_ExAC_classificaiton_sequences.bin: 100.00%
[Info::processExacFiles] Time taken 2.09m
```

Output Files:

|Filename | Filesize | Format | Info |
| -----:| -----:| -----:| -----:|
output/preprocess/CHRX/chrX_ExAC_spliced_mutation_data.csv|14.4MB|Text;CSV|Variants repeated as per splicing with other meta information.|
output/preprocess/CHRX/chrX_ExAC_sequence_ref.seq|759.1MB|Text|Reference genome spliced sequence.|
output/preprocess/CHRX/chrX_ExAC_sequence_alt.seq|759.1MB|Text|Reference genome spliced sequence with alternate allele.|
output/preprocess/CHRX/chrX_ExAC_sequence_ind.txt|1.9MB|Text|SNVs indexes in splice file|
output/preprocess/CHRX/chrX_ExAC_embedding_0.4_count.txt|5KB|Text|Switch counts in the data kept for learning embeddings.|
output/preprocess/CHRX/chrX_ExAc_embedding_0.4_sequences.bin|202.5MB|Binary;Custom|Switch sequences in binarizef form to learn embeddings.|
output/preprocess/CHRX/chrX_ExAC_embedding_index.txt|754.6KB|Text|Embedding sequence index relative to sequence(ref/alt) file.|
output/preprocess/CHRX/chrX_ExAC_classificaiton_sequences.bin|274.9MB|Binary;Custom|Switch sequences in binarized form.|
output/preprocess/CHRX/chrX_ExAC_classification_index.txt|986.7KB|Text|Sequence indexes in splice file (After dropping duplicates).|

To process __COSMIC__ file:

```bash
python src/preprocess.py -C ${CHROMOSOME} -c ${COSMIC_VCF_FILE} -u ${UCSC_KG_FILE} -U ${UCSC_KGxREF_FILE} -R ${REFERENCE_GENOME_FILE} -s ${SWITCH_FILE} -E ${EMBEDDING_FRACTION} -r ${RANDOM_STATE} -o ${PREPROCESS_OUTPUTPATH}
```

Expected output:

```
[Info::Preprocess] Extracting Genome
[Info::Preprocess] Time taken 16.52s
[Progress::processUCSCFiles] --- 100.00%
[Info::Preprocess] Time taken 3.09s
[Message::loadCOSMICVCF] Reading file data/CosmicCodingMuts_37_cosmic_build_89.vcf.gz [This step may take few minutes]
[Info::processCOSMICFiles] Time taken 12.89s
[Progress::createAlternateSplicingData] This step may take some time --- 100.00%
[Info::processCOSMICFiles] Time taken 8.55m
[Progress::getLengthInformation] This step may take some time --- 100.00%
[Info::processCOSMICFiles] Time taken 7.98m
[Progress::flattenSplicingData] --- 100.00%
[Info::processCOSMICFiles] Time taken 1.44m
[Progress::createSequences] This step may take some time -- 100.00%
[Info::processCOSMICFiles] Time taken 8.65m
[Progress::createSwitchSequences] --- 100.00%
[Info::removeDuplicates] Total Seqences 590171
[Progress::removeDuplicates] -- 100.00%
[Info::removeDuplicates] Remaining Sequences 281973. Reduction 47.78%
Saving metainfo output/preprocess/CHRX/chrX_COSMIC_classificaiton_sequences.bin: 100.00%
Saving sequences output/preprocess/CHRX/chrX_COSMIC_classificaiton_sequences.bin: 100.00%
[Info::processCOSMICFiles] Time taken 6.53m
```

Output Files:

|Filename | Filesize | Format | Info |
| -----:| -----:| -----:| -----:|
output/preprocess/CHRX/chrX_COSMIC_spliced_mutation_data.csv|75.3MB|Text;CSV|Variants repeated as per splicing with other meta information.|
output/preprocess/CHRX/chrX_COSMIC_sequence_ref.seq|1.8GB|Text|Reference genome spliced sequence.|
output/preprocess/CHRX/chrX_COSMIC_sequence_alt.seq|1.8GB|Text|Reference genome spliced sequence with alternate allele.|
output/preprocess/CHRX/chrX_COSMIC_sequence_ind.txt|4MB|Text|SNVs indexes in splice file|
output/preprocess/CHRX/chrX_COSMIC_classificaiton_sequences.bin|542.5MB|Binary;Custom|Switch sequences in binarized form.|
output/preprocess/CHRX/chrX_COSMIC_index.txt|1.9MB|Text|Sequence indexes in splice file (After dropping duplicates).|



To process __dbSNP__ file:

```bash
python src/preprocess.py -C ${CHROMOSOME} -d ${dbSNP_VCF_FILE} -u ${UCSC_KG_FILE} -U ${UCSC_KGxREF_FILE} -R ${REFERENCE_GENOME_FILE} -s ${SWITCH_FILE} -E ${EMBEDDING_FRACTION} -r ${RANDOM_STATE} -o ${PREPROCESS_OUTPUTPATH}
```

Expected output:

```
[Info::Preprocess] Extracting Genome
[Info::Preprocess] Time taken 17.52s
[Progress::processUCSCFiles] --- 100.00%
[Info::Preprocess] Time taken 3.12s
[Info::filterDBSNPVCF] Filtering chromosome X from VCF
[Info::filterDBSNPVCF] Extracting informatoin for chromosome X
[Info::filterDBSNPVCF] Filtering informatoin for chromosome X: Total filtered samples 24247152
[Info::processDBSNPFiles] Time taken 55.14m
[Progress::createAlternateSplicingData] This step may take some time --- 100.00%
[Info::processDBSNPFiles] Time taken 6.82h
[Progress::getLengthInformation] This step may take some time --- 100.00%
[Info::processDBSNPFiles] Time taken 8.47m
[Progress::flattenSplicingData] --- 100.00%
[Info::processDBSNPFiles] Time taken 1.61m
[Progress::createSequences] This step may take some time -- 100.00%
[Info::processDBSNPFiles] Time taken 10.46m
[Progress::createSwitchSequences] --- 100.00%
[Info::removeDuplicates] Total Seqences 798480
[Progress::removeDuplicates] -- 100.00%
[Info::removeDuplicates] Remaining Sequences 664740. Reduction 83.25%
Saving metainfo output/preprocess/chrX_dbSNP_classificaiton_sequences.bin: 100.00%
Saving sequences output/preprocess/chrX_dbSNP_classificaiton_sequences.bin: 100.00%
[Info::processDBSNPFiles] Time taken 9.13m

```

Output Files:


|Filename | Filesize | Format | Info |
| -----:| -----:| -----:| -----:|
output/preprocess/CHRX/temp_dbSNP_chrX.vcf|3.5GB|Text;VCF|CHROMOSOME variants from the provided vcf. Temporary file deleted after filtering.|
output/preprocess/CHRX/chrX_dbSNP_spliced_mutation_data.csv|59.7MB|Text;CSV|Variants repeated as per splicing with other meta information.|
output/preprocess/CHRX/chrX_dbSNP_sequence_ref.seq|2.1GB|Text|Reference genome spliced sequence.|
output/preprocess/CHRX/chrX_dbSNP_sequence_alt.seq|2.1GB|Text|Reference genome spliced sequence with alternate allele.|
output/preprocess/CHRX/chrX_dbSNP_sequence_ind.txt|5.5MB|Text|SNVs indexes in splice file|
output/preprocess/CHRX/chrX_dbSNP_classificaiton_sequences.bin|1.2GB|Binary;Custom|Switch sequences in binarized form.|
output/preprocess/CHRX/chrX_dbSNP_index.txt|4.6MB|Text|Sequence indexes in splice file (After dropping duplicates).|



To process __Met__ file:

```bash
python src/preprocess.py -C ${CHROMOSOME} -m ${mET_VCF_FILE} -u ${UCSC_KG_FILE} -U ${UCSC_KGxREF_FILE} -R ${REFERENCE_GENOME_FILE} -s ${SWITCH_FILE} -E ${EMBEDDING_FRACTION} -r ${RANDOM_STATE} -o ${PREPROCESS_OUTPUTPATH}
```

Expected output:

```
[Info::Preprocess] Extracting Genome
[Info::Preprocess] Time taken 16.89s
[Progress::processUCSCFiles] --- 100.00%
[Info::Preprocess] Time taken 3.17s
[Info::filterMETData] 12527 after splitting coordinate and keeping only SNVs
[Info::processMETFiles] Time taken 13.18s
[Progress::createAlternateSplicingData] This step may take some time --- 100.00%
[Info::processMETFiles] Time taken 1.60s
[Progress::getLengthInformation] This step may take some time --- 100.00%
[Info::processMETFiles] Time taken 2.14s
[Progress::flattenSplicingData] --- 100.00%
[Info::processMETFiles] Time taken 0.39s
[Progress::createSequences] This step may take some time -- 100.00%
[Info::processMETFiles] Time taken 2.49s
[Progress::createSwitchSequences] --- 100.00%
[Info::removeDuplicates] Total Seqences 2617
[Progress::removeDuplicates] -- 100.00%
[Info::removeDuplicates] Remaining Sequences 2193. Reduction 83.80%
Saving metainfo output/preprocess/CHRX/chrX_mET_classificaiton_sequences.bin: 100.00%
Saving sequences output/preprocess/CHRX/chrX_mET_classificaiton_sequences.bin: 100.00%
[Info::processMETFiles] Time taken 2.48s

```

Output Files:

|Filename | Filesitze | Format | Info |
| -----:| -----:| -----:| -----:|
output/preprocess/CHRX/chrX_mET_spliced_mutation_data.csv|362.1KB|Text;CSV|Variants repeated as per splicing with other meta information.|
output/preprocess/CHRX/chrX_mET_sequence_ref.seq|9.8MB|Text|Reference genome spliced sequence.|
output/preprocess/CHRX/chrX_mET_sequence_alt.seq|9.8MB|Text|Reference genome spliced sequence with alternate allele.|
output/preprocess/CHRX/chrX_mET_sequence_ind.txt|12.0KB|Text|SNVs indexes in splice file|
output/preprocess/CHRX/chrX_mET_classificaiton_sequences.bin|5.6MB|Binary;Custom|Switch sequences in binarized form.|
output/preprocess/CHRX/chrX_mET_index.txt|10KB|Text|Sequence indexes in splice file (After dropping duplicates).|


All the files can be processed in single command

```bash
python src/preprocess.py -C ${CHROMOSOME} -e {ExAC_VCF_FILE} -c ${COSMIC_VCF_FILE} -d ${dbSNP_VCF_FILE} -m ${mET_VCF_FILE} -u ${UCSC_KG_FILE} -U ${UCSC_KGxREF_FILE} -R ${REFERENCE_GENOME_FILE} -s ${SWITCH_FILE} -E ${EMBEDDING_FRACTION} -r ${RANDOM_STATE} -o ${PREPROCESS_OUTPUTPATH}
```

Output Files:

|Filename | Filesize | Format | Info |
| -----:| -----:| -----:| -----:|
|output/preprocess/CHRX/chrX*|11.8GB|Multiple|...|

<a name="cskipgram"></a>
<h4>2. Create skipgram </h4>

```bash
python src/create_skipgram.py 
usage: create_skipgram.py [-h] [-w WINDOW_SIZE] [-n NEGATIVE_SAMPLING] -c
                          SWITCH_COUNT_FILE [-r RANDOM_STATE] -S
                          SWITCH_SEQUENCE_FILE [-o OUTPUT_PATH] [-s]
```

|Parameter | Required or Optional| Datatype | Default Value | Description |
| -----:| -----:| -----:|-----:|-----:|
|-w|Optional|`int`|3|Window size|
|-n|Optional|`float`|0.2|Rate of negative sampling|
|-c|Required|`str`|None|Path of switch file|
|-r|Optional|`int`|0|Random state|
|-S|Required|`str`|None|Path of switch sequences file|
|-o|Optional|`str`|Current directory|Output path|
|-s|Optional|`bool`|`False`|If provided, then word-context-label file is saved|

```bash
WINDOW_SIZE=3
NEGATIVE_SAMPLING=0.2
SWITCH_COUNT_FILE=${PREPROCESS_OUTPUTPATH}/chrX_ExAC_embedding_${EMBEDDING_FRACTION}_count.txt
SWITCH_SEQUENCE_FILE=${PREPROCESS_OUTPUTPATH}/chrX_ExAC_embedding_${EMBEDDING_FRACTION}_sequences.bin
CREATED_SKIPGRAM_SAVE_PATH=output/skipgram/CHR${CHROMOSOME}

mkdir -p ${CREATED_SKIPGRAM_SAVE_PATH}

python src/create_skipgram.py -w ${WINDOW_SIZE} -n ${NEGATIVE_SAMPLING} -c ${SWITCH_COUNT_FILE} -r ${RANDOM_STATE} -S ${SWITCH_SEQUENCE_FILE} -o ${CREATED_SKIPGRAM_SAVE_PATH} -s
```

Expected output:

```
Loading file output/preprocess/CHRX/chrX_ExAC_embedding_0.4_sequences.bin
Arranging data 100.00%
[Progress::createSkipgram] --- This may take some time 100.00%
[Info::createSkipgram] Time taken 7.32m
[Info::saveSkipgramCounts] Saving counts to file output/skipgram/CHRX/tupleCount_ws3_ns0.2.csv
[Info::saveSkipgram] Total unique tuples/Total tuples: 197832/219656262
[Progress::saveSkipgram] Generating full list --- 100.00%
[Info::saveSkipgram] Shuffling array
[Info::saveSkipgram] Time taken 3.32m
[Info::saveSkipgram] Saving skipgrams to
[Progress::saveSkipgram] Saving group 1/1-(0,219656262)
        output/skipgram/CHRX/idx/skipgram_ws3_ns0.2_data_0_219656262.idx        output/skipgram/CHRX/idx/skipgram_ws3_ns0.2_label_0_219656262.idx
[Info::saveSkipgram] Saved file info at path output/skipgram/CHRX/idx/skipgram_ws3_ns0.2_fnames.csv
[Info::__main__] Total Time taken 14.09m
```

Output Files:

|Filename | Filesize | Format | Info |
| -----:| -----:| -----:| -----:|
|output/skipgram/CHRX/tupleCount_ws3_ns0.2.csv|2.4MB|Text;CSV|File containing tuple of `word`, `context`, `label` and `count`. Used for generating skipgram training file.|
|output/skipgram/CHRX/idx/skipgram_ws3_ns0.2_data_0_219656262.idx|878.6MB|Binary;IDX|Skipgram for training|
|output/skipgram/CHRX/idx/skipgram_ws3_ns0.2_label_0_219656262.idx|219.7MB|Binary;IDX|Target label for skipgram for training|
|output/skipgram/CHRX/idx/skipgram_ws3_ns0.2_fnames.csv|138bytes|Text;CSV|First row of file contains __vocabulary size__. Second row onwards name of training (data & labels) files.|


Note:

1. If skipgram training files (`-s`) are dump, these files are divided into chunks of 4GB files with minimum filesize of 1GB (while dividing). In case, if all data can be stored into less than 4GB files. Only one file will be written. Maximum file size is 5GB, i.e. if last chunk is less than 1 GB, it will be merged with penultimate chunck.

2. If `-s` option is not provided, script will write only tuple metadata. Metadata file contains 4 columns, `word`, `context`, `label` and `count`. `label=1` represents `word` and `context` were in window and `label=0` otherwise. `count` is the frequecny of `word`, `context` and `label` tuple. This information can be used to create full skipgram file.
<br/><br/>This option comes handy in those cases when we have very large chromosome. Possibly, Information present in this file can be used to downsample the the count to learn skipgram.

<a name="tskipgram"></a>
<h4>3. Train skipgram </h4>

```bash
python src/train_skipgram.py 
usage: train_skipgram.py [-h] -t TUPLE_COUNT_INFO_FILE [-e EMBEDDING_SIZE]
                         [-p PEEK_INTERVAL] [-S START_EPOCH] [-E END_EPOCH]
                         [-b BATCH_SIZE] [-o OUTPUT_PATH]
```

|Parameter | Required or Optional| Datatype | Default Value | Description |
| -----:| -----:| -----:|-----:|-----:|
|-t|Required|`str`|None|Path of file for tuple information|
|-e|Optional|`int`|300|Size of embedding vector|
|-p|Optional|`int`|10|Step size for saving|
|-S|Optional|`int`|0|Start epoch for training|
|-E|Optional|`int`|200|End epoch for training|
|-b|Optional|`int`|8192|Training batch size|
|-o|Optional|`str`|Current directory|Output path|

```bash
TUPLE_COUNT_INFO_FILE=${CREATED_SKIPGRAM_SAVE_PATH}/idx/skipgram_ws3_ns0.2_fnames.csv
EMBEDDING_SIZE=300
EMB_PEEK_INTERVAL=10
EMB_START_EPOCH=1
EMB_END_EPOCH=200
EMB_BATCH_SIZE=8192
EMB_OUTPUTPATH=output/skipgram/CHR${CHROMOSOME}/embedding_${EMBEDDING_SIZE}

mkdir -p ${EMB_OUTPUTPATH}

CUDA_VISIBLE_DEVICES=0 python src/train_skipgram.py -t ${TUPLE_COUNT_INFO_FILE} -e ${EMBEDDING_SIZE} -p ${EMB_PEEK_INTERVAL} -S ${EMB_START_EPOCH} -E ${EMB_END_EPOCH} -b ${EMB_BATCH_SIZE} -o ${EMB_OUTPUTPATH}
```

Expected Output:
```
[Progress::loadSkipgram] Loading group 1/1
__________________________________________________________________________________________________
Layer (type)                    Output Shape         Param #     Connected to                     
==================================================================================================
input_1 (InputLayer)            [(None, 1)]          0                                            
__________________________________________________________________________________________________
input_2 (InputLayer)            [(None, 1)]          0                                            
__________________________________________________________________________________________________
embedding (Embedding)           (None, 1, 300)       192000      input_1[0][0]                    
                                                                 input_2[0][0]                    
__________________________________________________________________________________________________
reshape (Reshape)               (None, 300, 1)       0           embedding[0][0]                  
__________________________________________________________________________________________________
reshape_1 (Reshape)             (None, 300, 1)       0           embedding[1][0]                  
__________________________________________________________________________________________________
dot (Dot)                       (None, 1, 1)         0           reshape[0][0]                    
                                                                 reshape_1[0][0]                  
__________________________________________________________________________________________________
reshape_2 (Reshape)             (None, 1)            0           dot[0][0]                        
__________________________________________________________________________________________________
dense (Dense)                   (None, 1)            2           reshape_2[0][0]                  
==================================================================================================
Total params: 192,002
Trainable params: 192,002
Non-trainable params: 0
__________________________________________________________________________________________________
None
Epoch 1/200
Train on 219656262 samples
219656262/219656262 [==============================] - 302s 1us/sample - loss: 0.0958
Epoch 2/200
Train on 219656262 samples
51617792/219656262 [======>.......................] - ETA: 4:14 - loss: 0.0945
.
.
.
```

Output Files:

|Filename | Filesize | Format | Info |
| -----:| -----:| -----:| -----:|
|output/skipgram/CHRX/embedding_300/midway/weights_emb300_itr_*.pkl|768.2KB|Binary;Pickle|Intermediate embedding weights, Can be used to investigate embeddings in the middle of training. These files are generated on the `EMB_PEEK_INTERVAL`|
|output/skipgram/CHRX/embedding_300/finish/weights_emb300_itr_200.pkl|768.2KB|Binary;Pickle|Trained embeddings to be fed into the classifier|

![](images/embedding.png?raw=true)

2D UMAP representation of 300 sized vectors learned for switches.

1. Train embeddings in text format are present in the file `models/skipgram_weights_emb300_itr_200.txt`. (Order of embeddings are same as switch file `data/snp_switch_dict.csv`)

2. UMAP projections of these figures are present in the file `models/embedding_umap_projections.txt`

<a name="trtesplit"></a>
<h4>4. Train-Test split </h4>

```bash
python src/create_train_test_split.py 
usage: create_train_test_split.py [-h] -e EXAC_SPLICED_VARIANT_FILE -E
                     768.2KB|Binary;Pickle|768.2KB|Binary;Pickle|             EXAC_SEQUENCES_FILE -i EXAC_INDEX_FILE -c
                                  COSMIC_SPLICED_VARIANT_FILE -C
                                  COSMIC_SEQUENCES_FILE -I COSMIC_INDEX_FILE
                                  [-f TRAIN_SPLIT_FRACTION]
                                  [-l MAX_SEQUENCE_LENGTH]
                                  [-m MIN_GENE_FREQUENCEY] [-r RANDOM_STATE]
                                  -a AMINO_ACID_MAP_FILE -s SWITCH_FILE [-R]
                                  [-u] [-o OUTPUT_PATH]
```

|Parameter | Required or Optional| Datatype | Default Value | Description |
| -----:| -----:| -----:|-----:|-----:|
|-e|Required|`str`|None|Path of ExAC VCF with alternative splicing|
|-E|Required|`str`|None|Path of ExAC sequences file|
|-i|Required|`str`|None|Path of indexes for EXAC sequences|
|-c|Required|`str`|None|Path of COSMIC VCF with alternative splicing|
|-C|Required|`str`|None|Path of COMSIC sequences file|
|-I|Required|`str`|None|Path of indexes for COSMIC sequneces|
|-f|Optional|`float`|0.7|if 0 < v < 1 Ratio of train-test split, if v > 2 folds, v==1 is not valid|
|-l|Optional|`int`|1500|Max length of sequences considered|
|-m|Optional|`int`|200|Min variations in a gene|
|-r|Optional|`int`|0|Random state|
|-a|Required|`str`|None|Path of Amino acid map|
|-s|Required|`str`|None|Path of switch file|
|-R|Optional|`bool`|False|If provided, removes synonymous mutations|
|-u|Optional|`bool`|False|If provided, removes label noise|
|-o|Optional|`str`|Current directory|Output path|

```bash
ExAC_SPLICED_VARIANT_FILE=${PREPROCESS_OUTPUTPATH}/chrX_ExAC_spliced_variants_data.csv
ExAC_SEQUENCE_FILE=${PREPROCESS_OUTPUTPATH}/chrX_ExAC_classificaiton_sequences.bin
ExAC_INDEX_FILE=${PREPROCESS_OUTPUTPATH}/chrX_ExAC_classificaiton_index.txt
COSMIC_SPLICED_VARIANT_FILE=${PREPROCESS_OUTPUTPATH}/chrX_COSMIC_spliced_mutation_data.csv
COSMIC_SEQUENCE_FILE=${PREPROCESS_OUTPUTPATH}/chrX_COSMIC_classificaiton_sequences.bin
COSMIC_INDEX_FILE=${PREPROCESS_OUTPUTPATH}/chrX_COSMIC_index.txt
TRAIN_SPLIT_FRACTION=0.7
MAX_SEQUENCE_LENGTH=1500
MIN_GENE_FREQUENCEY=200
AMINO_ACID_MAP_FILE=data/aminoMap.csv
TRAIN_TEST_SPLIT_OUTPUT_FILE=output/train_test_split/CHR${CHROMOSOME}

mkdir -p ${TRAIN_TEST_SPLIT_OUTPUT_FILE}

python src/create_train_test_split.py -e ${ExAC_SPLICED_VARIANT_FILE} -E ${ExAC_SEQUENCE_FILE} -i ${ExAC_INDEX_FILE} -c ${COSMIC_SPLICED_VARIANT_FILE} -C ${COSMIC_SEQUENCE_FILE} -I ${COSMIC_INDEX_FILE} -f ${TRAIN_SPLIT_FRACTION} -l ${MAX_SEQUENCE_LENGTH} -m ${MIN_GENE_FREQUENCEY} -r ${RANDOM_STATE} -a ${AMINO_ACID_MAP_FILE} -s ${SWITCH_FILE} -o ${TRAIN_TEST_SPLIT_OUTPUT_FILE} -R -u
```

Expected Output:

```
[Info::createTrainTestSplit] Loading files
Loading file output/preprocess/CHRX/chrX_ExAC_classificaiton_sequences.bin
Arranging data 100.00%
Loading file output/preprocess/CHRX/chrX_COSMIC_classificaiton_sequences.bin
Arranging data 100.00%
[Info::createTrainTestSplit] Remove synonymous from train data
[Progress::removeSynonymousMutations] -- 100.00%
[Info::removeSynonymousMutations] Total percentage missense and nonsense mutations 73.19%
[Info::createTrainTestSplit] Time taken 3.99m
[Info::createTrainTestSplit] Remove synonymous from test data
[Progress::removeSynonymousMutations] -- 100.00%
[Info::removeSynonymousMutations] Total percentage missense and nonsense mutations 72.62%
[Info::createTrainTestSplit] Time taken 1.78m
[Info::createTrainTestSplit] Remove label noise from train data
[Progress::removeLabelNoise] --- 100.00%
[Info::removeLabelNoise] Keep percentage after removing label noise 95.98%
[Info::createTrainTestSplit] Time taken 9.24s
[Info::createTrainTestSplit] Remove label noise from test data
[Progress::removeLabelNoise] --- 100.00%
[Info::removeLabelNoise] Keep percentage after removing label noise 95.51%
[Info::createTrainTestSplit] Time taken 4.60s
[Info::createTrainTestSplit] Remove sequences longer than 1500 from train data
[Info::createTrainTestSplit] Time taken 0.09s
[Info::createTrainTestSplit] Remove sequences longer than 1500 from test data
[Info::createTrainTestSplit] Time taken 0.04s
[Info::createTrainTestSplit] Creating Train groups
[Info::createTrainTestSplit] Time taken 0.06s
[Info::createTrainTestSplit] Creating Test groups
[Info::createTrainTestSplit] Time taken 0.01s
Saving metainfo output/train_test_split/CHRX/train_data.bin: 100.00%
Saving sequences output/train_test_split/CHRX/train_data.bin: 100.00%
Saving metainfo output/train_test_split/CHRX/test_data.bin: 100.00%
Saving sequences output/train_test_split/CHRX/test_data.bin: 100.00%
[Info::createTrainTestSplit] Train Gene 232	Train splice variants 918	Train samples 146680	Train imbalance 0.7033
[Info::createTrainTestSplit] Test Gene 101	Test splice variants 400	Test samples 60518	 test imbalance 0.6670
```

Output Files:

|Filename | Filesize | Format | Info |
| -----:| -----:| -----:| -----:|
|output/train_test_split/CHRX/train_data.bin|225.9MB|Binary;Custom|Switch sequences for training classifier|
|output/train_test_split/CHRX/train_info.csv|15.3MB|Text;CSV|Splicing/Gene/Label info for training data|
|output/train_test_split/CHRX/train_group.pkl|1.2MB|Binary;Pickle|Group information for creating training batches| 
|output/train_test_split/CHRX/test_data.bin|86.4MB|Binary;Custom|Switch sequences for testing/validation of classifier|
|output/train_test_split/CHRX/test_info.csv|6.1MB|Test;CSV|Splicing/Gene/Label info for testing data|
|output/train_test_split/CHRX/test_group.pkl|499.8KB|Binary;Pickle|Group information for creating testing batches|

Note:

1. Because of the randomization in creating embedding splits and train-test splits, output of the code might be different from what used in original study.
2. To insure reproducibility, filtered files have also been provided. (rows marked as `included=1` are used in study.)

|Filename | Filesize | Format | Info | Gene (in-study/Total) | splicing (in-study/Total) | Samples (in-study/Total)
| -----:| -----:| -----:| -----:| -----:| -----:| -----:|
|output/train_geneinfo_chrX.csv.gz|783.4KB|Text;CSV(gzip)|Train dataset used in study|234/240|930/981|143924/182002|
|output/test_geneinfo_chrX.csv.gz|276.3KB|Text;CSV(gzip)|Test dataset used in study|98/103|375/388|55263/64397|

<a name="train"></a>
<h4>5. Train deep network </h4>

```bash
python src/train_classifier.py
usage: train_classifier.py [-h] -e EMBEDDING_FILE -t TRAIN_SEQUENCES_FILE -l
                           TRAIN_LABEL_FILE -v VALIDATION_SEQUENCES_FILE -L
                           VALIDATION_LABEL_FILE -T TRAIN_GROUP_FILE -V
                           VALIDATION_GROUP_FILE [-f FILE_TYPE]
                           [-S START_EPOCH] [-E END_EPOCH] [-p PEEK_INTERVAL]
                           [-P VALIDATION_INTERVAL] [-b BATCH_SIZE]
                           [-w WEIGHT_REGULARIZER] [-a ATTENTION_REGULARIZER]
                           [-o OUTPUT_PATH]
```

|Parameter | Required or Optional| Datatype | Default Value | Description |
| -----:| -----:| -----:|-----:|-----:|
|-e|Required|`str`|None|Path of learned embeddings|
|-t|Required|`str`|None|Path of training sequences|
|-l|Required|`str`|None|path of training labels|
|-v|Required|`str`|None|Path of validation sequences|
|-L|Required|`str`|None|Path of validation labels|
|-T|Required|`str`|None|Path of training groups|
|-V|Required|`str`|None|Path of validation groups|
|-f|Optional|`str`|'pkl' `or` 'txt'|File type of embedding file| 
|-S|Optional|`int`|0|Start epoch for training|
|-E|Optional|`int`|200|End epoch for training|
|-p|Optional|`int`|1|Training step size for saving|
|-P|Optional|`int`|1|Validation step size for saving|
|-b|Optional|`int`|25|Training batch size|
|-w|Optional|`float`|1e-5|Bi-LSTM weight regularization|
|-a|Optional|`float`|1e-6|Attention layer weight regularization|
|-o|Optional|`str`|Current directory|Output path|

The weights of trained network are available at: <a href="https://drive.google.com/file/d/1QLzHQqzZ8Zzoq8BdMPMNJB6PkKqwu0Lx/view?usp=sharing">link</a>

```bash
EMBEDDING_FILE_PATH=${EMB_OUTPUTPATH}/finish/weights_emb${EMBEDDING_SIZE}_itr_${EMB_END_EPOCH}.pkl
TRAIN_SEQUENCE_FILE=${TRAIN_TEST_SPLIT_OUTPUT_FILE}/train_data.bin
TRAIN_LABEL_FILE=${TRAIN_TEST_SPLIT_OUTPUT_FILE}/train_info.csv
VALIDATION_SEQUENCE_FILE=${TRAIN_TEST_SPLIT_OUTPUT_FILE}/test_data.bin
VALIDATION_LABEL_FILE=${TRAIN_TEST_SPLIT_OUTPUT_FILE}/test_info.csv
TRAIN_GROUP_FILE=${TRAIN_TEST_SPLIT_OUTPUT_FILE}/train_group.pkl
VALIDATION_GROUP_FILE=${TRAIN_TEST_SPLIT_OUTPUT_FILE}/test_group.pkl
TRAINING_EMBEDDING_FILE_TYPE='pkl'
TRAIN_START_EPOCH=1
TRAIN_END_EPOCH=200
TRAIN_PEEK_INTERVAL=1
TRAIN_VALIDATION_INTERVAL=1
TRAIN_BATCH_SIZE=25
WEIGHT_REG=1e-5
ATN_REG=1e-6
TRAIN_OUTPUT_PATH=outupt/train/CHR${CHROMOSOME}

CUDA_VISIBLE_DEVICES=0 python src/train_classifier.py -e ${EMBEDDING_FILE_PATH} -t ${TRAIN_SEQUENCE_FILE} -l ${TRAIN_LABEL_FILE} -v ${VALIDATION_SEQUENCE_FILE} -L ${VALIDATION_LABEL_FILE} -T ${TRAIN_GROUP_FILE} -V ${VALIDATION_GROUP_FILE} -f ${TRAINING_EMBEDDING_FILE_TYPE} -S ${TRAIN_START_EPOCH} -E ${TRAIN_END_EPOCH} -p ${TRAIN_PEEK_INTERVAL} -p ${TRAIN_VALIDATION_INTERVAL} -b ${TRAIN_BATCH_SIZE} -w ${WEIGHT_REG} -a ${ATN_REG} -o ${TRAIN_OUTPUT_PATH}

```

Expected output:

```
Loading file output/train_test_split/CHRX/train_data.bin
Arranging data 100.00%
Loading file output/train_test_split/CHRX/test_data.bin
Arranging data 100.00%
[Progress::createBatch] --- 100.00%
[Progress::createBatch] --- 100.00%
Layer (type)                 Output Shape              Param #   
=================================================================
embedding (Embedding)        (None, None, 300)         192000    
_________________________________________________________________
bidirectional (Bidirectional (None, None, 600)         1444800   
_________________________________________________________________
batch_normalization (BatchNo (None, None, 600)         2400      
_________________________________________________________________
bidirectional_1 (Bidirection (None, None, 600)         2164800   
_________________________________________________________________
batch_normalization_1 (Batch (None, None, 600)         2400      
_________________________________________________________________
time_distributed (TimeDistri (None, None, 100)         60100     
_________________________________________________________________
batch_normalization_2 (Batch (None, None, 100)         400       
_________________________________________________________________
atnlayer (AttentionWithConte (None, 100)               10200     
_________________________________________________________________
dense_1 (Dense)              (None, 1)                 101       
=================================================================
Total params: 3,877,201
Trainable params: 3,682,601
Non-trainable params: 194,600
_________________________________________________________________
None
Train epochs 1/200
3/3353 [..............................] - ETA: 5:00:34 - loss: 0.7860 - acc: 0.6319
.
.
.
```

Output Files:

|Filename | Filesize | Format | Info |
| -----:| -----:| -----:| -----:|
|outupt/train/CHRX/model_epoch*/weights.pkl|15.5MB|Binary;Pickle|Model weights at certain epoch, * will be decided by `TRAIN_PEEK_INTERVAL`|
|outupt/train/CHRX/model_epoch*/train_preds.csv|-|Text;CSV|Train Prediction at * epoch|
|outupt/train/CHRX/model_epoch*/train_metric.txt|-|Text|Different metrics for train prediction at * epoch|
|outupt/train/CHRX/model_epoch*/valid_preds.csv|-|Text;CSV|Validation prediction at * epoch, * is decided by `TRAIN_VALIDATION_INTERVAL`|
|outupt/train/CHRX/model_epoch*/valid_metric.txt|-|Text|Different metrics for validation prediction at '*' epoch|

<a name="create"></a>
<h4>6. Create prediction data </h4>

```bash
python src/create_prediction_data.py 
usage: create_prediction_data.py [-h] -f SPLICE_VARIANT_FILE -S SEQUENCE_FILE
                                 [-i INDEX_FILE] [-b BASE_LABEL] -t
                                 TRAIN_SEQUENCE_FILE -a AMINO_ACID_MAP_FILE -s
                                 SWITCH_FILE [-R] [-l MAX_SEQUENCE_LENGTH]
                                 [-o OUTPUT_PATH]
```

|Parameter | Required or Optional| Datatype | Default Value | Description |
| -----:| -----:| -----:|-----:|-----:|
|-f|Required|`str`|None|Path of splice variant file|
|-S|Required|`str`|None|Path of sequence file|
|-i|Optional|`str`|None|Path of Index file|
|-b|Optional|`int`|None|labels to be assigned to all samples in|
|-t|Required|`str`|None|Path of training sequence file|
|-a|Required|`str`|None|Path of Amino acid map|
|-s|Required|`str`|None|Path of switch file|
|-R|Optional|`bool`|False|If provided, remove synonymous mutations|
|-l|Optional|`int`|1500|Max length of sequence considered|
|-o|Optional|`str`|Current directory|Output path|

Note: 
1. If index file is not provided (`-i`) an error will be raised if number of lines in splice variant file (`-f`) and sequence file is not same (`-S`).
2. if base labels (`-b`) not provided, it is expected that splice variant file (`-f`) will have `label` columns.


To create __dbSNP__ data:

```bash
dbSNP_SPLICED_VARIANT_FILE=${PREPROCESS_OUTPUTPATH}/chrX_dbSNP_spliced_mutation_data.csv
dbSNP_SEQUENCE_FILE=${PREPROCESS_OUTPUTPATH}/chrX_dbSNP_classificaiton_sequences.bin
dbSNP_INDEX_FILE=${PREPROCESS_OUTPUTPATH}/chrX_dbSNP_index.txt
dbSNP_BASE_LABEL=0
dbSNP_OUTPUT_PATH=output/prediction_data/CHR${CHROMOSOME}/dbSNP/

mkdir -p ${dbSNP_OUTPUT_PATH}

python src/create_prediction_data.py -f ${dbSNP_SPLICED_VARIANT_FILE} -S ${dbSNP_SEQUENCE_FILE} -i ${dbSNP_INDEX_FILE} -b ${dbSNP_BASE_LABEL} -t ${TRAIN_SEQUENCE_FILE} -a ${AMINO_ACID_MAP_FILE} -s ${SWITCH_FILE} -l ${MAX_SEQUENCE_LENGTH} -o ${dbSNP_OUTPUT_PATH} -R
```

Expected output:

```
Loading file output/preprocess/CHRX/chrX_dbSNP_classificaiton_sequences.bin
Arranging data 100.00%
Loading file output/train_test_split/CHRX/train_data.bin
Arranging data 100.00%
[Info::createPredictionData] Removing synonymous mutations
[Progress::removeSynonymousMutations] -- 100.00%
[Info::removeSynonymousMutations] Total percentage missense and nonsense mutations 66.32%
[Info::createPredictionData] Time taken 9.20m
[Info::createPredictionData] Removing sequences larget than 1500
[Info::createPredictionData] Time taken 0.36s
[Info::createPredictionData] Removing overlap with training data
[Info::removeTrainOverlap] Keep Fraction 84.38%
[Info::createPredictionData] Time taken 16.91s
[Info::createPredictionData] Creating groups
[Info::createPredictionData] Time taken 0.14s
Saving metainfo output/prediction_data/CHRX/dbSNP/data.bin: 100.00%
Saving sequences output/prediction_data/CHRX/dbSNP/data.bin: 100.00%
[Info::createPredictionData] Total remaining samples 316582
```

Output Files:

|Filename | Filesize | Format | Info |
| -----:| -----:| -----:| -----:|
|output/prediction_data/CHRX/dbSNP/data.bin|399.5MB|Binary;Custom|Switch sequences for predictions|
|output/prediction_data/CHRX/dbSNP/info.csv|27.4MB|Text;CSV|Splicing/Gene/Label information for sequences|
|output/prediction_data/CHRX/dbSNP/groups.pkl|2.6MB|Binary;Pickle|Group information for batch creation during testing|

To create __mET__ data:

```bash
mET_SPLICED_VARIANT_FILE=${PREPROCESS_OUTPUTPATH}/chrX_mET_spliced_mutation_data.csv
mET_SEQUENCE_FILE=${PREPROCESS_OUTPUTPATH}/chrX_mET_classificaiton_sequences.bin
mET_INDEX_FILE=${PREPROCESS_OUTPUTPATH}/chrX_mET_index.txt
mET_BASE_LABEL=1
mET_OUTPUT_PATH=output/prediction_data/CHR${CHROMOSOME}/mET/

mkdir -p ${mET_OUTPUT_PATH}

python src/create_prediction_data.py -f ${mET_SPLICED_VARIANT_FILE} -S ${mET_SEQUENCE_FILE} -i ${mET_INDEX_FILE} -b ${mET_BASE_LABEL} -t ${TRAIN_SEQUENCE_FILE} -a ${AMINO_ACID_MAP_FILE} -s ${SWITCH_FILE} -l ${MAX_SEQUENCE_LENGTH} -o ${mET_OUTPUT_PATH} -R
```

Expected output:

```
Loading file output/preprocess/CHRX/chrX_mET_classificaiton_sequences.bin
Arranging data 100.00%
Loading file output/train_test_split/CHRX/train_data.bin
Arranging data 100.00%
[Info::createPredictionData] Removing synonymous mutations
[Progress::removeSynonymousMutations] -- 100.00%
[Info::removeSynonymousMutations] Total percentage missense and nonsense mutations 99.82%
[Info::createPredictionData] Time taken 2.44s
[Info::createPredictionData] Removing sequences larget than 1500
[Info::createPredictionData] Time taken 0.00s
[Info::createPredictionData] Removing overlap with training data
[Info::removeTrainOverlap] Keep Fraction 80.64%
[Info::createPredictionData] Time taken 1.58s
[Info::createPredictionData] Creating groups
[Info::createPredictionData] Time taken 0.00s
Saving metainfo output/prediction_data/CHRX/mET/data.bin: 100.00%
Saving sequences output/prediction_data/CHRX/mET/data.bin: 100.00%
[Info::createPredictionData] Total remaining samples 1216
```

Output Files:

|Filename | Filesize | Format | Info |
| -----:| -----:| -----:| -----:|
|output/prediction_data/CHRX/mET/data.bin|2.2MB|Binary;Custom|Switch sequences for predictions|
|output/prediction_data/CHRX/mET/info.csv|184.5KB|Text;CSV|Splicing/Gene/Label information for sequences|
|output/prediction_data/CHRX/mET/groups.pkl|12.9KB|Binary;Pickle|Group information for batch creation during testing|

Note:

1. Because of the randomization in creating embedding splits and train-test splits, output of the code might be different from what used in original study.
2. To insure reproducibility, filtered files have also been provided. (After all filtering)

|Filename | Filesize | Format | Info | Gene | splicing | Samples)|
| -----:| -----:| -----:| -----:| -----:| -----:| -----:|
|output/info_dbsnp.csv.gz|6.7MB|Text;CSV(gzip)|Processed dbsnp variants info with splicing|772|1820|315040|
|output/info_mET.csv.gz|20.1kb|Text;CSV(gzip)|Processed dbsnp variants info with splicing|19|71|1220|

<a name="prediction"></a>
<h4>7. Prediction </h4>

```bash
python src/predict.py
usage: predict.py [-h] -e EMBEDDING_FILE -m MODEL_FILE -t SEQUENCES_FILE
                  [-f FILE_TYPE] -l LABEL_FILE -T GROUP_FILE [-b BATCH_SIZE]
                  [-w WEIGHT_REGULARIZER] [-a ATTENTION_REGULARIZER]
                  [-o OUTPUT_PATH]
```

|Parameter | Required or Optional| Datatype | Default Value | Description |
| -----:| -----:| -----:|-----:|-----:|
|-e|Required|`str`|None|Path of embedding file|
|-m|Required|`str`|None|Path of model file|
|-t|Required|`str`|None|Path of sequence file|
|-l|Required|`str`|None|Path of label file (splice variant file)|
|-T|Required|`str`|None|Path of group file|
|-b|Optional|`int`|25|Testing batch size|
|-w|Optional|`float`|1e-5|BI-LSTM weight regularization|
|-a|Optional|`float`|1e-6|Attention layer weight regulatization|
|-o|Optional|`str`|Current directory|Output path|

```bash
MODEL_FILE=${TRAIN_OUTPUT_PATH}/model_epoch200/weights.pkl
mET_SEQUENCE_FILE=${mET_OUTPUT_PATH}/data.bin
mET_LABEL_FILE=${mET_OUTPUT_PATH}/info.csv
mET_GROUP_FILE=${mET_OUTPUT_PATH}/groups.pkl
EMBEDDING_FILE_TYPE=pkl
PREDICT_BATCH_SIZE=25
PREDICT_OUTPUT_PATH=output/predictions/CHR${CHROMOSOME}/mET

mkdir -p ${PREDICT_OUTPUT_PATH}

CUDA_VISIBLE_DEVICES=0 python src/predict.py -e ${EMBEDDING_FILE_PATH} -m ${MODEL_FILE} -t ${mET_SEQUENCE_FILE} -f ${EMBEDDING_FILE_TYPE} -l ${mET_LABEL_FILE} -T ${mET_GROUP_FILE} -b ${PREDICT_BATCH_SIZE} -w ${WEIGHT_REG} -a ${ATN_REG} -o ${PREDICT_OUTPUT_PATH}
```

Output Files

|Filename | Filesize | Format | Info |
| -----:| -----:| -----:| -----:|
|output/predictions/CHRX/mET/prediction_metric.txt|-|Text|Different metrics for test prediction|
|output/predictions/CHRX/mET/predictions.csv|-|Text;CSV|Test prediction|

Similar steps for `dbSNP` or any other dataset.

<a name="pretrain"></a>
<h4>8. Predict with pre-trained model </h4>

Steps:

1. Download trained weights from <a href="https://drive.google.com/file/d/1QLzHQqzZ8Zzoq8BdMPMNJB6PkKqwu0Lx/view?usp=sharing">link</a> and place them in models folder and run

```bash
gunzip weights_model.pkl.gz
```

2. Run the blow snippet.

```bash
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
```

This script show how to get predictions for `mET` data. Similarly, Other data can be passed to get their predictions.
