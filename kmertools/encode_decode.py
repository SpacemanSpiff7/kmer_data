import multiprocessing as mp
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
    args = split_seq(nprocs, sequence)
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


def kmer_search(sequence, kmer_length):
    counts = Counter()
    for i in range(len(sequence) - (kmer_length - 1)):  # This takes the (1-based) reference sequence for chromosome 22
        next_seq = sequence[i:(i + kmer_length)]
        if not ('N' in next_seq or 'n' in next_seq):
            counts[next_seq] += 1
    return counts


def get_kmer_count(sequence, kmer_length):
    args = split_seq(sequence, mp.cpu_count(), overlap=kmer_length)
    args = [[seq, kmer_length] for seq in args]
    pool = mp.Pool(mp.cpu_count())
    results = [res.get() for res in [pool.starmap_async(kmer_search, args)]]
    pool.close()
    counts = Counter()
    for result in results[0]:
        for k, v in result.items():
            counts[k] += v
    return counts
