import psutil
from multiprocessing import Pool


def parallel_variant(vcf_path, nthreads=0):
    if nthreads == 0:
        nthreads = psutil.cpu_count(logical=False)
    pool = Pool(nthreads)
