import sys
import os
import numpy as np
import sklearn.metrics as skme

def mkdir(path):
    if not os.path.exists(path):
        os.mkdir(path)

ncbi_chr_id = {'Y':'NC_000024.9', 'X':'NC_000023.10'}

#https://www.garykessler.net/library/file_sigs.html
#https://stackoverflow.com/questions/13044562/python-mechanism-to-identify-compressed-file-type-and-uncompress
magic_dict = {
    b"\x1f\x8b\x08": "gz",
    b"\x42\x5a\x68": "bz2",
    b"\x50\x4b\x03\x04": "zip"
    }

max_len = max(len(x) for x in magic_dict)

def getFileType(fname):
    with open(fname, 'rb') as f:
        file_start = f.read(max_len)
    for magic, filetype in magic_dict.items():
        if file_start.startswith(magic):
            return filetype
    return "txt"

def getCommentLinesFromGzip(file_identifier, comment):

    count = 0
    while(1):
        l = file_identifier.readline().decode()[:len(comment)]
        if l == comment:
            count += 1
        else:
            break

    return count

def getCommentLinesFromText(file_identifier, comment):

    count = 0
    while(1):
        l = file_identifier.readline()[:len(comment)]
        if l == comment:
            count += 1
        else:
            break

    return count

def removeSynonymousMutations(data, info_df, aa_dict, switch_dict):

    selected_data = []
    selected_ind = []
    n_sequences = len(data)

    for i, a in enumerate(data):
        print(f'[Progress::removeSynonymousMutations] -- {100*(i + 1)/n_sequences:0.2f}%', end='\r')
        for b in a:
            s = switch_dict[b]
            ref = s[:3]
            alt = s[3:]
            if ref != alt:
                if aa_dict[ref] != aa_dict[alt]:
                    selected_data.append(a)
                    selected_ind.append(i)
                    break
    print('')
    selected_info_df = info_df.iloc[selected_ind].reset_index(drop=True)

    print(f'[Info::removeSynonymousMutations] Total percentage missense and nonsense mutations {100*len(selected_data)/len(data):0.2f}%')
    return selected_data, selected_info_df

def removeLognerSequence(data, info_df, max_length):

    lengths = info_df['lengths'].values.astype(np.int)
    keep_inds = np.where(lengths < max_length)[0]

    selected_df = []
    for i in keep_inds:
        selected_df.append(data[i])

    selected_info_df = info_df.iloc[keep_inds].reset_index(drop=True)

    return selected_df, selected_info_df

def removeLabelNoise(data, info_df):

    total_count = len(data)

    index_set_train = []
    lengths = info_df['lengths'].values.astype(np.int)
    labels = info_df['label'].values.astype(np.int)
    unique_lengths = np.unique(lengths)

    for u in unique_lengths:
        ind = np.where(lengths == u)[0]
        index_set_train.append(ind)

    total_groups = len(index_set_train)
    keep_ind = []
    for i, ind in enumerate(index_set_train):

        print(f'[Progress::removeLabelNoise] --- {100*(i+1)/total_groups:0.2f}%', end='\r')

        local_db = []
        for j in ind:
            local_db.append(data[j])
        local_db = np.array(local_db)
        uq_local_db, return_index, inverse, return_counts = np.unique(local_db, axis=0, return_index=True, return_inverse=True, return_counts=True)

        u_ind = ind[return_index[np.where(return_counts == 1)[0]]].tolist()
        nu_ind_ = return_index[np.where(return_counts > 1)[0]]

        for nu in nu_ind_:
            locs = ind[np.where(inverse == nu)[0]]
            c = np.sum(labels[locs])
            if c == locs.shape[0] or c == 0:
                u_ind.append(ind[nu])

        keep_ind += u_ind
    print('')

    keep_count = len(keep_ind)

    selected_seqs = []
    for i in keep_ind:
        selected_seqs.append(data[i])
    selected_info_df = info_df.iloc[keep_ind].reset_index(drop=True)

    print(f'[Info::removeLabelNoise] Keep percentage after removing label noise {100*keep_count/total_count:0.2f}%')

    return selected_seqs, selected_info_df


def removeTrainOverlap(seqs, info_df, tr_seqs):

    tr_len = np.array([len(s) for s in tr_seqs])
    d_len = info_df['lengths'].values.astype(np.int)

    unique_len = np.unique(d_len)

    selected_data = []
    selected_inds = []

    keep_count = 0
    for u in unique_len:
        tr_ind = np.where(tr_len == u)[0]
        s_ind = np.where(d_len == u)[0]

        if len(tr_ind) > 0:
            ltdb = []
            for t in tr_ind:
                ltdb.append(tr_seqs[t])

            for t in s_ind:
                ltdb.append(seqs[t])

            ltdb = np.array(ltdb, dtype=np.int)
            uq_ltdb, return_index, return_counts = np.unique(ltdb, axis=0, return_index=True, return_counts=True)
            rc_m = return_index[np.where(return_counts == 1)[0]]
            rc_l = rc_m[np.where(rc_m >= len(tr_ind))[0]]

            if len(rc_l) > 0:
                lind = s_ind[rc_l - len(tr_ind)].tolist()
                selected_inds += lind
                keep_count += len(lind)
        else:
            selected_inds += s_ind.tolist()
            keep_count += len(s_ind)

    selected_data = [seqs[i] for i in selected_inds]
    selected_info_df = info_df.iloc[selected_inds].reset_index(drop=True)

    print(f'[Info::removeTrainOverlap] Keep Fraction {100*keep_count/len(seqs):0.2f}%')

    return selected_data, selected_info_df

def createGroups(info_df):
    groups = []

    lengths = info_df['lengths'].values.astype(np.int)
    uq_lengths = np.unique(lengths)

    for u in uq_lengths:
        ind = np.where(lengths == u)[0]
        groups.append(ind)

    return groups

def createBatch(data, label, group, max_seq_length, min_batch_size):
    group_by_data = []

    for e, g in enumerate(group):
        print(f'[Progress::createBatch] --- {100*(e+1)/len(group):0.2f}%', end='\r')

        t = []
        l = []

        max_batch_size = (max_seq_length * min_batch_size) // (len(data[g[0]]))

        #print(len(g), max_batch_size)
        if len(g) > max_batch_size:
            for i in range(len(g)):
                t.append(data[g[i]].astype(np.float32))
                l.append([label[g[i]]])

                if len(t) == max_batch_size or i == len(g) - 1:
                    t = np.array(t, dtype=np.float32)
                    l = np.array(l, dtype=np.float32)
                    group_by_data.append((t, l))

                    t = []
                    l = []
        else:

            for i in g:
                t.append(data[i].astype(np.float32))
                l.append([label[i]])
            #print([len(_) for _ in t])
            #print('before t')
            t = np.array(t, dtype=np.float32)
            #print('after t')
            l = np.array(l, dtype=np.float32)

            group_by_data.append((t, l))

    print('')
    return group_by_data

def createBatchFromDF(data, info_df, group, max_seq_length, min_batch_size):

    group_by_data = []

    for e, g in enumerate(group):
        print(f'[Progress::createBatchFromDF] -- {100*(e+1)/len(group):0.2f}%', end='\r')
        t = []
        l = []

        max_batch_size = (max_seq_length * min_batch_size) //  (len(data[g[0]]))

        if len(g) > max_batch_size:
            for i in range(len(g)):
                t.append(data[g[i]].astype(np.float32))
                l.append(info_df.iloc[g[i]])

                if len(t) == max_batch_size or i == len(g) - 1:
                    t = np.array(t, dtype=np.float32)
                    group_by_data.append((t, l))

                    t = []
                    l = []
        else:

            for i in g:
                t.append(data[i].astype(np.float32))
                l.append(info_df.iloc[i])

            t = np.array(t, dtype=np.float32)

            group_by_data.append((t, l))

    print('')
    return group_by_data

def get_generator(group_by_data, test=False):
    count = -1
    max_count = len(group_by_data)
    while(True):
        count += 1
        if(count >= max_count):
            count = 0
        x = group_by_data[count][0]
        y = group_by_data[count][1]
        if test:
            yield x
        else:
            yield x, y

def get_generator_shuffle(group_by_data):

    max_count = len(group_by_data)
    perm = np.random.permutation(max_count)
    count = -1
    while(True):
        count += 1
        if(count >= max_count):
            count = 0

        x = group_by_data[perm[count]][0]
        y = group_by_data[perm[count]][1]

        yield x, y

def computeMetrics(true_labels, predicted_labels, metric_out_file, th=0.5):

    if not isinstance(true_labels, np.ndarray):
        _t = np.array(true_labels)
    else:
        _t = true_labels.copy()

    if not isinstance(predicted_labels, np.ndarray):
        _p = np.array(predicted_labels)
    else:
        _p = predicted_labels.copy()

    mse = np.sqrt(np.mean(np.square(_t - _p)))

    a = skme.classification_report(_t, (_p>=th).astype(np.int), labels=[0, 1], target_names=['normal', 'cancer'], digits=4)
    print(a)

    print('MSE {}'.format(mse))

    with open(metric_out_file, 'w') as fid:
        fid.write(a)

        fid.write('\nMSE : {}\n'.format(mse))

def processTime(time):
    s = 's'
    if time > 60:
        time = time / 60
        s = 'm'
    if time > 60:
        time = time / 60
        s = 'h'
    return f'{time:0.2f}{s}'
