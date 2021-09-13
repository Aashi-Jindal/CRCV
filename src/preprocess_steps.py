import numpy as np
import pandas as pd
from Bio import SeqIO, Seq
import utils as mu
import gzip
import mio as io
import subprocess
import os

def splitRefGenome(ref_file, chromosome):

    for rec in SeqIO.parse(ref_file, 'fasta'):
        #print(f'[Split Reference Genome]Chromosome {rec.id}')
        if rec.id == chromosome:
            return str(rec.seq)

def processUCSCFiles(kg_file, kg_ref_file, chromosome):

    exon_loc = pd.read_csv(kg_file, delimiter='\t', low_memory=False)
    exon_ref = pd.read_csv(kg_ref_file, delimiter='\t', low_memory=False)

    req_cols = ['#name', 'geneSymbol', 'chrom', 'strand', 'txStart', 'txEnd', 'cdsStart', 'cdsEnd', 'exonCount', 'proteinID', 'alignID']

    desc = exon_ref['description']
    desc_split = np.array([a.split(',')[-1] for a in desc])
    mRNA_locs = np.where(desc_split == ' mRNA.')[0]
    exon_ref_sel = exon_ref.iloc[mRNA_locs].reset_index(drop=True)
    exon_loc_chr = exon_loc[exon_loc['chrom'] == f'chr{chromosome}'].reset_index(drop=True)
    exon_mrna = exon_ref_sel.merge(exon_loc_chr, left_on='#kgID', right_on='#name', how='inner')

    values = []
    for i in range(exon_mrna.shape[0]):
        print(f'[Progress::processUCSCFiles] --- {100*(i+1)/exon_mrna.shape[0]:0.2f}%', end='\r')
        vals = exon_mrna[req_cols].iloc[i].values
        _strt = exon_mrna['exonStarts'].iloc[i].split(',')[:-1]
        _end = exon_mrna['exonEnds'].iloc[i].split(',')[:-1]
        if len(_strt) == len(_end):
            for _s in range(len(_strt)):
                values.append(np.append(vals, [_strt[_s], _end[_s]]))
        else:
            mrna_name = exon_mrna['#kgID'].iloc[i].values[0]
            print(f'[ERROR::processUCSCFiles] start length != end length of exon at {mrna_name}')

    print('')
    processed_info = pd.DataFrame(values, columns=req_cols+['exonStart', 'exonEnd'])
    processed_info['exonStart'] = processed_info['exonStart'].astype(np.int)
    processed_info['exonEnd'] = processed_info['exonEnd'].astype(np.int)
    return processed_info

def filterExACVCF(fname, chromosome, csq_ind_from_end=13, exon_ind=8, comment='##'):

    fileformat = mu.getFileType(fname)

    if fileformat == 'gz':
        with gzip.open(fname, 'rb') as fid:
            comment_lines = mu.getCommentLinesFromGzip(fid, comment)
    elif fileformat == 'txt':
        with open(fname, 'r') as fid:
            comment_lines = mu.getCommentLinesFromText(fid, comment)
    else:
        raise ValueError(f'File format is not supported, Supported formats are gzip and txt. Provided file {fname}')

    print(f'[Message::filterExACVCF] Reading file {fname} [This step may take few minutes]')
    vcf_data = pd.read_table(fname, skiprows=comment_lines, low_memory=False)

    print(f'[Message::filterExACVCF] Filtering ...')
    vcf_data_filter = vcf_data[vcf_data['FILTER'] == 'PASS'].reset_index(drop=True)
    vcf_data_chr = vcf_data_filter[vcf_data_filter['#CHROM'] == chromosome].reset_index(drop=True)
    vcf_filter_info = vcf_data_chr['INFO'].values

    extract_info_csq = [a.split(';')[-csq_ind_from_end] for a in vcf_filter_info]
    len_info_csq = np.array([len(a.split('|')) for a in extract_info_csq])
    keep_info_csq_ind = np.where(len_info_csq > 1)[0]

    exon_vals = np.array([extract_info_csq[_sel_ind].split('|')[exon_ind] for _sel_ind in keep_info_csq_ind])
    filter_locs = keep_info_csq_ind[np.where(exon_vals != '')[0]]
    header = ['#CHROM', 'POS', 'ID', 'REF', 'ALT']
    filtered_data = vcf_data_chr.iloc[filter_locs][header].reset_index(drop=True)

    return filtered_data

def loadCOSMICVCF(fname, chromosome, comment='##'):

    fileformat = mu.getFileType(fname)

    if fileformat == 'gz':
        with gzip.open(fname, 'rb') as fid:
            comment_lines = mu.getCommentLinesFromGzip(fid, comment)
    elif fileformat == 'txt':
        with open(fname, 'r') as fid:
            comment_lines = mu.getCommentLinesFromText(fid, comment)
    else:
        raise ValueError(f'File format is not supported, Supported formats are gzip and txt. Provided file {fname}')

    print(f'[Message::loadCOSMICVCF] Reading file {fname} [This step may take few minutes]')
    vcf = pd.read_table(fname, skiprows=comment_lines, low_memory=False, delimiter="\t")
    vcf_sel = vcf[vcf['#CHROM'] == chromosome].reset_index(drop=True)

    return vcf_sel

def filterDBSNPVCF(fname, chromosome, temp_out_path):

    #Inner function
    def getval(s1, s2):

        for v in s1:
            if s2 in v:
                s = v.split('=')
                return s[-1]

    #zcat GCF_000001405.25.gz | grep NC_000024 > chrY.vcf
    fileformat = mu.getFileType(fname)
    temp_out_file = os.path.join(temp_out_path, f'temp_dbSNP_chr{chromosome}.vcf')

    if fileformat == 'gz':
        command = 'zcat'
    else:
        command = 'cat'

    #https://stackoverflow.com/questions/13332268/how-to-use-subprocess-command-with-pipes
    print(f'[Info::filterDBSNPVCF] Filtering chromosome {chromosome} from VCF')
    with open(temp_out_file, 'w') as fid:
        ps = subprocess.Popen((command, fname), stdout=subprocess.PIPE)
        subprocess.call(('grep', mu.ncbi_chr_id[chromosome]), stdin=ps.stdout, stdout=fid)

    header = ['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'pub', 'common', 'clnvi', 'clnorigin', 'clnsig', 'clndisdb', 'clndn', 'clnrevstat', 'clnacc']

    print(f'[Info::filterDBSNPVCF] Extracting informatoin for chromosome {chromosome}')
    data = []
    with open(temp_out_file, 'r') as fid:
        while True:
            line = fid.readline()

            if not line:
                break

            line = line.strip('\n')
            c = [chromosome]
            split_line = line.split('\t')
            info = split_line[-1]
            info_split = info.split(';')

            info_keep = []

            if 'VC=' in info:
                vc = getval(info_split, 'VC=')
                if vc != 'SNV':
                    continue
            else:
                print(line)

            for info_part in ['PUB', 'COMMON', 'CLNVI=', 'CLNORIGIN=', 'CLNSIG=', 'CLNDISDB=', 'CLNDN=', 'CLNREVSTAT=', 'CLNACC=']:

                if info_part in info:
                    if info_part in ['PUB', 'COMMON']:
                        info_keep.append(str(1))
                    else:
                        info_keep.append(getval(info_split, info_part))
                else:
                    info_keep.append('')


            c.append(split_line[1])
            c.append(split_line[2])
            c.append(split_line[3])
            c.append(split_line[4])
            c += info_keep

            data.append(c)

    content = pd.DataFrame(data, columns=header)
    all_vals = np.array([f'{i}' for i in range(256)], dtype=np.str)
    keep_vals = np.array(['0', '1', '2', '3', '255'], dtype=np.str)
    toremove_vals = np.setdiff1d(all_vals, keep_vals)

    keep_content = content[content['clnsig'].isna()]

    to_process = content[np.logical_not(content['clnsig'].isna())]
    keep_indexes = []
    for i in range(to_process.shape[0]):

        v = to_process.iloc[i]['clnsig']
        sv = v.split(',')
        ssv = []
        for _sv in sv:
            __ssv = _sv.split('|')
            ssv += __ssv
        ins = np.intersect1d(ssv, toremove_vals)
        if len(ins) == 0:
            keep_indexes.append(i)

    further_keep_content = to_process.iloc[keep_indexes]

    filtered_data = pd.concat([keep_content, further_keep_content], axis=0).reset_index(drop=True)

    print(f'[Info::filterDBSNPVCF] Filtering informatoin for chromosome {chromosome}: Total filtered samples {filtered_data.shape[0]}')

    filtered_data['POS'] = filtered_data['POS'].astype(np.int)

    os.unlink(temp_out_file)
    return filtered_data

def filterMETData(fname, chromosome):

    data = pd.read_excel(fname)

    index_coord = np.where(np.array(data.columns) == 'coordinate')[0][0]
    index_vc = np.where(np.array(data.columns) == 'category')[0][0]

    values = []
    for i in range(data.shape[0]):
        coord = data.iloc[i]['coordinate']
        if not data.iloc[i].isna()['coordinate']:
            coord_l = coord.split(',')
            if len(coord) > 1:
                vc_l = data.iloc[i][index_vc].split(',')
                for j in range(len(coord_l)):
                    _vc = vc_l[j]
                    if _vc == 'SNV':
                        values.append(data.iloc[i][:index_coord].tolist() + [coord_l[j]] + [_vc] + data.iloc[i][index_vc+1:].tolist())
            else:
                if data.iloc[i][index_vc] == 'SNV':
                    values.append(data.iloc[i].tolist())

    print(f'[Info::filterMETData] {len(values)} after splitting coordinate and keeping only SNVs')

    splitted_values = []

    for i in range(len(values)):
        coord = values[i][index_coord]
        _chr, _mut = coord.split(':')
        _p1, _alt = _mut.split('>')
        _pos, _ref = _p1[:-len(_alt)], _p1[-len(_alt)]

        splitted_values.append(values[i] + [_chr, _pos, _ref, _alt])

    _df = pd.DataFrame(splitted_values, columns = data.columns.tolist() + ['#CHROM', 'POS', 'REF', 'ALT'])
    _df = _df[_df['#CHROM'] == chromosome].reset_index(drop=True)
    _df['POS'] = _df['POS'].astype(np.int)
    return _df

def createAlternateSplicingData(filtered_vcf, mrna_locs):

    values = []
    for i in range(filtered_vcf.shape[0]):
        print(f'[Progress::createAlternateSplicingData] This step may take some time --- {100*(i+1)/filtered_vcf.shape[0]:0.2f}%', end='\r')
        if filtered_vcf['ALT'].iloc[i][0] == 'N': #This is a very weak condition, as Alteration may be more the one allele long. However, this will not pose further problem as only SNPs will be kept.
            #print(filtered_vcf['ALT'].iloc[i])
            continue

        _pos = filtered_vcf['POS'].iloc[i]

        cdsStart_ind = np.where(mrna_locs['cdsStart'] < _pos)[0]
        cdsEnd_ind = np.where(mrna_locs['cdsEnd'] >= _pos)[0]
        cdscommon_ind = np.intersect1d(cdsStart_ind, cdsEnd_ind)

        if cdscommon_ind.shape[0] > 0:
            cdsRegion = mrna_locs.iloc[cdscommon_ind]
            exonStart_ind = np.where(cdsRegion['exonStart'] < _pos)[0]
            exonEnd_ind = np.where(cdsRegion['exonEnd'] >= _pos)[0]
            exoncommon_ind = np.intersect1d(exonStart_ind, exonEnd_ind)

            if exoncommon_ind.shape[0] > 0:
                exonRegion = cdsRegion.iloc[exoncommon_ind]
                name_strand = np.unique(exonRegion[['#name', 'strand', 'geneSymbol']].values.astype(np.str), axis=0)

                v = np.array([filtered_vcf.iloc[i].values]* name_strand.shape[0])
                values.append(np.append(v, name_strand, axis=1))
    print('')

    values = np.vstack(values)
    df = pd.DataFrame(values, columns=filtered_vcf.columns.tolist()+['#name', 'ucsc_strand', 'Gene'])
    return df

def getLengthInformation(splicing_data, mrna_locs):

    values = []

    groups = mrna_locs.groupby(by='#name')
    for i in range(splicing_data.shape[0]):
        print(f'[Progress::getLengthInformation] This step may take some time --- {100*(i+1)/splicing_data.shape[0]:0.2f}%', end='\r')

        name = splicing_data['#name'].iloc[i]
        group = groups.get_group(name)
        cdsStrt = group['cdsStart'].iloc[0]
        cdsEnd = group['cdsEnd'].iloc[0]

        total_len = 0
        for j in range(group.shape[0]):

            _str = group['exonStart'].iloc[j]
            _end = group['exonEnd'].iloc[j]

            if _end >= cdsStrt and _str <= cdsEnd:
                if _str < cdsStrt:
                    _str = cdsStrt
                if _end > cdsEnd:
                    _end = cdsEnd

                total_len += _end - _str

        if total_len%3 != 0:
            print(f'{name} not multiple of 3')
        else:
            values.append(np.append(splicing_data.iloc[i].values, [total_len]))
    print('')
    len_df = pd.DataFrame(values, columns=list(splicing_data.columns)+['length'])
    return len_df

def flattenSplicingData(splicing_data, output_file):

    values = []
    alt_column = np.where(splicing_data.columns == 'ALT')[0][0]
    for i in range(splicing_data.shape[0]):
        print(f'[Progress::flattenSplicingData] --- {100*(i+1)/splicing_data.shape[0]:0.2f}%', end='\r')
        _m = splicing_data.iloc[i]
        _ref = _m['REF']
        _alt = _m['ALT'].split(',')
        _m_l = list(_m.values)

        if len(_alt) > 1:
            for j in range(len(_alt)):
                values.append(_m_l[:alt_column] + [_alt[j]] + _m_l[alt_column+1:])
        else:
            values.append(_m_l)

    print('')
    df = pd.DataFrame(values, columns=splicing_data.columns)
    df.to_csv(output_file)
    return df

def createSequences(splicing_data, mrna_locs, reference_genome, output_fstub):

    output_file_ref = f'{output_fstub}_ref.seq'
    output_file_alt = f'{output_fstub}_alt.seq'
    output_file_ind = f'{output_fstub}_ind.txt'
    index_processed = []

    groups = mrna_locs.groupby('#name')

    with open(output_file_ref, 'w') as ref_fid:
        with open(output_file_alt, 'w') as alt_fid:

            for i in range(splicing_data.shape[0]):
                print(f'[Progress::createSequences] This step may take some time -- {100*(i+1)/splicing_data.shape[0]:0.2f}%', end='\r')

                _m = splicing_data.iloc[i]
                _pos = _m['POS']
                _ref = _m['REF']
                _alt = _m['ALT']
                _len = _m['length']
                _name = _m['#name']
                _strand = _m['ucsc_strand']

                group = groups.get_group(_name)
                cdsStrt = group['cdsStart'].iloc[0]
                cdsEnd = group['cdsEnd'].iloc[0]
                _flag = False

                if len(_ref) > 1 and len(_alt) > 1: #SNP, del
                    if len(_ref) == len(_alt):
                        if _ref[1:] == _alt[1:]: #then SNP
                            _ref = _ref[0]
                            _alt = _ref[0]
                            _flag = True
                    else: #del
                        pass
                elif len(_ref) > 1: #del only
                    pass
                elif len(_ref) == 1 and len(_alt) > 1: #ins only
                    pass
                elif len(_ref) == 1 and len(_alt) == 1: #SNPs only
                    _flag = True
                else:
                    print(f'[Error::createSequences] missed combination REF {_ref} -> ALT {_alt} at POS {_m} (index {i})')


                if _flag:
                    _sub_str = ''
                    first = True

                    for j in range(group.shape[0]):

                        _str = group['exonStart'].iloc[j]
                        _end = group['exonEnd'].iloc[j]

                        if _end >= cdsStrt and _str <= cdsEnd:
                            if _str < cdsStrt:
                                _str = cdsStrt
                            if _end > cdsEnd:
                                _end = cdsEnd

                            if first:
                                if _pos > _str and _pos <= _end:
                                    _pos_sub_str = len(_sub_str) + _pos - 1 - _str
                                    first = False

                            _sub_str += reference_genome[_str:_end]

                    if len(_sub_str) == _len:
                        if _sub_str[_pos_sub_str] == _ref:
                            index_processed.append([i])
                            _sub_str_alt = _sub_str[:_pos_sub_str] + _alt + _sub_str[_pos_sub_str+1:]

                            if _strand == '-':
                                _sub_str = str(Seq.Seq(_sub_str).reverse_complement())
                                _sub_str_alt = str(Seq.Seq(_sub_str_alt).reverse_complement())

                            ref_fid.write(f'{_sub_str}\n')
                            alt_fid.write(f'{_sub_str_alt}\n')

                        else:
                            print(f'[Error::createSequences] position mismatch, {_pos_sub_str}, {_sub_str[_pos_sub_str]}, {_ref}')
                    else:
                        print(f'[Error::createSequences] mismatch in string found and id _sub_str len {len(_sub_str)} and {_len}')

    print('')
    np.savetxt(output_file_ind, index_processed, fmt='%d')
    return np.array(index_processed, dtype=np.int)

def createSwitchSequences(sequence_fstub, indexes, out_fstub, switch_dict, remove_duplicates=False, dump_count=False):

    ref_file = f'{sequence_fstub}_ref.seq'
    alt_file = f'{sequence_fstub}_alt.seq'

    sequences = []

    if dump_count:
        counts = np.zeros(len(switch_dict), dtype=np.int)

    with open(ref_file, 'r') as ref_fid:
        with open(alt_file, 'r') as alt_fid:

            curr_line = 0
            for i, ind in enumerate(indexes):
                print(f'[Progress::createSwitchSequences] --- {100*(i+1)/len(indexes):0.2f}%', end='\r')

                while(curr_line <= ind):
                    ref = ref_fid.readline().strip('\n')
                    alt = alt_fid.readline().strip('\n')
                    curr_line += 1

                _bin_index = []

                for k in range(0, len(ref), 3):
                    switch = ref[k:k+3]+alt[k:k+3]
                    _ind = switch_dict[switch]
                    _bin_index.append(_ind)
                    if dump_count:
                        counts[_ind] += 1
                sequences.append(_bin_index)

    print('')
    if remove_duplicates:
        sequences, kept_indexes = removeDuplicates(sequences, indexes)

    bin_fname = f'{out_fstub}_sequences.bin'
    io.saveBin(bin_fname, sequences)

    if dump_count:
        count_fname = f'{out_fstub}_count.txt'
        io.totxtWithIndex(count_fname, counts, header=['switch_id','count'])

    if remove_duplicates:
        return kept_indexes

def removeDuplicates(sequences, indexes):

    print(f'[Info::removeDuplicates] Total Seqences {indexes.shape[0]}')

    lengths = np.array([len(s) for s in sequences])
    unique_lengths = np.unique(lengths)

    kept_indexes = []
    kept_sequences = []
    for u, ul in enumerate(unique_lengths):
        print(f'[Progress::removeDuplicates] -- {100*(u+1)/len(unique_lengths):0.2f}%', end='\r')

        ind = np.where(lengths == ul)[0]
        arr = np.zeros((ind.shape[0], ul), dtype=np.int)
        for i, j in enumerate(ind):
            arr[i] = sequences[j]

        unique_ar, r_indexes = np.unique(arr, return_index=True, axis=0)
        kept_sequences += unique_ar.tolist()
        kept_indexes += ind[r_indexes].tolist()

    print(f'\n[Info::removeDuplicates] Remaining Sequences {len(kept_indexes)}. Reduction {100 * len(kept_indexes)/indexes.shape[0]:0.2f}%')
    return kept_sequences, kept_indexes
