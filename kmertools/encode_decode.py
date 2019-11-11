import multiprocessing as mp
import time
from collections import Counter


def encoder(sequence):
    encd = {
        'A': 1,
        'C': 2,
        'G': 3,
        'T': 4,
        'N': 0
    }
    result = []
    for c in list(sequence):
        result.append(encd[c])
    return result


def decoder(sequence):
    dcd = {
        1: 'A',
        2: 'C',
        3: 'G',
        4: 'T',
        0: 'N'
    }
    result = ""
    for i in sequence:
        result += dcd[int(i)]
    return result


def encode_sequence(sequence):
    pool = mp.Pool()
    nprocs = mp.cpu_count()
    args = split_seq(sequence, nprocs)
    results = [funccall.get() for funccall in [pool.map_async(encoder, args)]]
    pool.close()
    encoded_seq = []
    for seq_list in results[0]:
        encoded_seq.extend(seq_list)
    return encoded_seq


def split_seq(sequence, nprocs, overlap=None):
    chunk_size = int(len(sequence) / nprocs) + 1
    args = []
    start = 0
    end = chunk_size
    for proc in range(nprocs):
        if overlap is not None:
            args.append(sequence[start:(end + overlap - 1)])
        else:
            args.append(sequence[start:end])
        start = end
        end += chunk_size
        if end > len(sequence):
            end = len(sequence)
    return args


def decode_sequence(sequence):
    pool = mp.Pool()
    nprocs = mp.cpu_count()
    chunk_size = int(len(sequence) / nprocs) + 1
    args = []
    start = 0
    end = chunk_size
    for proc in range(nprocs):
        args.append(sequence[start:end])
        start = end
        end += chunk_size
        if end > len(sequence):
            end = len(sequence)
    results = [funccall.get() for funccall in [pool.map_async(encoder, args)]]
    pool.close()
    decoded_seq = ""
    for s in results[0]:
        decoded_seq += s
    return decoded_seq


def convert_seq(sequence, encode=True):
    if encode:
        return encode_sequence(sequence)
    else:
        return decode_sequence(sequence)



