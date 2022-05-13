import concurrent.futures
import numpy as np

from datetime import timedelta
from functools import partial
from multiprocessing import Process, Queue
from multiprocessing.dummy import Pool
from sys import argv
from time import time
from ROOT import TFile


S_RUN = 'run26' if len(argv) < 2 else argv[1]  # S_RUN is the first argument of the script, 'run26' if not given
TREE_NAME = 'tree'
N_ADC = 192  # the number of ADCs
N_THREADS = 1 if len(argv) < 3 else int(argv[2])  # N_THREADS is the second argument, 1 if not given
DEFAULT_THRESHOLD = 100  # a default noise threshold for the amplitude on an ADC
DEFAULT_N_CUTOFF = 800  # a default number of hits detected on an ADC that would indicate that most of the hits at that amplitude are noise
NOISE_SIGNAL_RATIO = 5  # the ratio of something deemed to be noise compared to a regular signal
ROLLING_AVG_OFFSET = 15  # how far ahead of the current point the algorithm looks to determine the signal strength
ROLLING_AVG_HALF_WIDTH = 5  # how many points ahead of and behind the avg center point to use to calculate the rolling avg

q = Queue()  # stores the return data from the count function

# open file and get tree
file = TFile(f'./data/{S_RUN}.root')
t = file.Get(TREE_NAME)
n_entries = t.GetEntriesFast()
file.Close()

"""
A function to follow the status of program execution, gives: elapsed time, estimate of remaining time, number of events evaluated up to that point
"""
def print_time(start_time, step, total):
    elapsed_time = round(time() - start_time)
    percentage = step/total
    remaining_time = 0
    if percentage != 0:
        remaining_time = round(elapsed_time * (1-percentage)/percentage)
        remaining_time = timedelta(0, remaining_time)
    elapsed_time = timedelta(0, elapsed_time)
    print(f"{percentage*100:.1f}%, elapsed {elapsed_time}, remaining {remaining_time}, event {step} out of {total}")



"""
Given the start and end point of the data, returns the amplitudes of all the recorded events per ADC
"""
def count(start, end):
    # loads the raw data root file in once per thread
    # some systems may experience RAM issues when trying to run too many threads
    root_file = TFile(f'./data/{S_RUN}.root')
    tree = root_file.Get(TREE_NAME)

    arr = [ [] for i in range(N_ADC) ]
    for ind in range(start, end):
        if start == 0 and ind % 10000 == 0:
            print_time(start_time, ind*N_THREADS, n_entries)

        tree.GetEntry(ind)

        for j in range(tree.wm):
            arr[tree.wadc[j] - 1].append(tree.wampl[j])

    root_file.Close()

    # stores the data in a global variable which is later used to merge the data
    # it would be better practice for q to be an argument of the count function instead of a global variable
    q.put(arr)


"""
Takes a tuple of the number of the ADC and an array of all the amplitudes measured on that ADC
and uses a simple algorithm to calculate the appropriate noise threshold for that ADC
"""
def get_threshold(pack):
    i, out_data = pack
    ampl_count = {}
    for data_point in out_data:
        if data_point not in ampl_count.keys():
            ampl_count[data_point] = 0

        ampl_count[data_point] += 1

    ### for testing purposes only - slows the script down a lot ###
    # with open(f'adc/{S_RUN}_{i+1}.txt', 'w') as f:
    #     for key in sorted(ampl_count.keys()):
    #         f.write(f'{key}: {ampl_count[key]}\n')

    #     f.close()

    ks, vs = list(ampl_count.keys()), list(ampl_count.values())
    k_max = ks[vs.index(max(vs))]

    avg = lambda l: sum(l)/len(l)

    if ampl_count[k_max] >= DEFAULT_N_CUTOFF:
        for k in range(k_max, len(ks)-10):
            k_start = k + ROLLING_AVG_OFFSET - ROLLING_AVG_HALF_WIDTH
            k_end = k + ROLLING_AVG_OFFSET + ROLLING_AVG_HALF_WIDTH
            a = [ ampl_count[key] for key in range(k_start, k_end) ]
            if ampl_count[k] < DEFAULT_N_CUTOFF and ampl_count[k]/avg(a) <= NOISE_SIGNAL_RATIO:
                return k
    return DEFAULT_THRESHOLD


if __name__ == '__main__':
    # execute the count function on all the data, split into N_THREADS threads
    print(f'Counting amplitudes on {N_THREADS} threads...')
    start_time = time()

    ps = []
    for i in range(N_THREADS):
        time1 = time()
        print(f'Starting thread {i+1}')

        # every process does 1//N_THREADS of the calculations
        p = Process(target=count, args=(n_entries*i//N_THREADS, n_entries*(i+1)//N_THREADS))
        p.start()
        ps.append(p)

        print(f'Thread {i+1} started in {timedelta(0, time() - time1)}')

    ampls = [ [] for i in range(N_ADC) ]
    # get all the return values from the count functions and merge the data
    for _ in range(N_THREADS):
        ret = q.get()
        for i in range(N_ADC):
            ampls[i].extend(ret[i])

    # wait for all the processes to finish
    for p in ps:
        p.join()

    print(f'Amplitudes counted successfully in {timedelta(0, time() - start_time)}')

    # calculate the noise thresholds based on the amplitude data

    thresholds = []
    print('Calculating thresholds...')
    time1 = time()

    pool = Pool(N_THREADS)
    thresholds = pool.map(get_threshold, enumerate(ampls))

    print(f'Thresholds calculated successfully in {timedelta(0, time() - time1)}')

    # outputs the calculated thresholds to a file to be used for noise filtering
    print('Writing thresholds...')
    time1 = time()
    with open(f'./thresholds/noise_thresholds_{S_RUN}', 'w') as f:
        for i, t in enumerate(thresholds):
            f.write(f'{i+1}: {t}\n')

    print(f'Thresholds written successfully in {timedelta(0, time() - time1)}')
