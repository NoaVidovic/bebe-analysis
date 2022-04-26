# -*- coding: utf-8 -*-
"""
@author: Deni

This is a script for removing noise form the data. Primarily used at the start of the analyis to cut the amount of (faulty) data. 
Prior to use we need to define noise thresholds for the channels, either a universal threshold or a specific one for each channel. 
"""

""" import predefined functions """
from ROOT import TFile, TTree
from numpy import zeros, loadtxt
from time import time
from datetime import timedelta
from multiprocessing import Process, Queue
from sys import argv



""" 
Define program-wide info
    !ADJUST!: S_RUN, TREE_NAME, N_ADC, ALL(True,False), CHANNEL_LIST, UNIVERSAL_CUTOFF(True,False)
"""
S_RUN = 'run26' if len(argv) < 2 else argv[1] # name of the run, string form
TREE_NAME = 'tree'  # usually we used names "tree" or "T"
N_ADC = 192  # number of ADC channels (strips)
N_THREADS = 4
ALL = True  # do we remove noise in ALL of the channels (TRUE) or just some (FALSE)
CHANNEL_LIST = []  # channels for which we do NOT remove noise
UNIVERSAL_CUTOFF = False  # do we have a single amplitude cutoff for all channels (TRUE), or do we have a specific cutoff for each channel (FALSE)


if ALL:
    S_SUFFIX = 'killnoiseALL'  # suffix to add to the output files 
else:
    S_SUFFIX = 'killnoiseSOME'  # suffix to add to the output files 
    
    
    
"""
A function to follow the status of program execution, gives: elapsed time, estimate of remaining time, number of events evaluated up to that point
"""
def print_time(start_time, step, total):
    elapsed_time = round(time()-start_time)
    percentage = (i*1.)/n_entries
    remaining_time = 0
    if percentage != 0:
        remaining_time = round(elapsed_time * (1-percentage)/percentage)
        remaining_time = timedelta(0, remaining_time)
    elapsed_time = timedelta(0, elapsed_time)
    print("{:.1f}%, elapsed {}, remaining {}, event {} out of {} ".format(percentage*100, elapsed_time, remaining_time, i, n_entries))


""" 
Define cutoffs: either specific for each channel or universal
    !ADJUST!: file input of specific cutoffs, or a single value for the universal cutoff
"""
if not UNIVERSAL_CUTOFF:  # input from file noise_thresholds_BeBe.txt
    noise_cutoff = loadtxt('noise_thresholds_{}'.format(S_RUN),
                           delimiter=': ',
                           dtype={'names': ('adc', 'weight'),
                           'formats': ('i2', 'i2')})
if UNIVERSAL_CUTOFF:
    noise_cutoff = 160



"""
A function to determine if a recorded hit in a channel is noise or not
"""
def is_noise(adc, amplitude, cutoff):
    
    if UNIVERSAL_CUTOFF:
        if ALL or (not ALL and adc not in CHANNEL_LIST):
            if amplitude < cutoff:
                return True
        elif adc in CHANNEL_LIST:
            return True

    else:
        if ALL or (not ALL and adc not in CHANNEL_LIST):
            if amplitude < cutoff[adc-1][1]:
                return True
        elif adc in CHANNEL_LIST:
            return True
        
    return False
    
    
        
"""open the file and get the tree"""
myfile = TFile(f'{S_RUN}.root')  # open
# myfile.ls() #check the file contents
t = myfile.Get(TREE_NAME)   # get the tree
# t.Print() #check the tree contents    
    
    

"""creating the output root file; defining name, output tree, output tree variables, output tree structure
"""
out_file = TFile(f'{S_RUN}_{S_SUFFIX}.root', "recreate")  # output root file
out_t = TTree('tree', f'{S_RUN}_{S_SUFFIX}')  # output tree


# output tree variables
wm = zeros(1, dtype='short')  # multiplicity of an event - number of hits, value of type Short
wadc = zeros(N_ADC, dtype='short')  # channel number, field of values of type Short
wampl = zeros(N_ADC, dtype='short')  # amplitude, field of values of type Short

# output tree structure; identical to the structure of the input tree
out_t.Branch('wm', wm, 'wm/S')  # multiplicity branch
out_t.Branch('wadc', wadc, 'wadc[wm]/S')  # channel number branch, field of lenght "wm" 
out_t.Branch('wampl', wampl, 'wampl[wm]/S')  # amplitude branch, field of lenght "wm" 


""" go trough the run data """
n_entries = t.GetEntriesFast()  # get the  number of entries in the tree

start_time = time()  # start the timer

for i in range(n_entries):  # for: 0 to n_entries-1

    if i % 10000 == 0:  # print the timer info every 10000 entries
        print_time(start_time, i, n_entries)
    
    t.GetEntry(i)  # get an entry(event)

    wm[0] = 0  # setting the number of hits which are not noise to 0
    
    for j in range(t.wm):  # going through all the hits in an entry(event)
    
        if not is_noise(t.wadc[j], t.wampl[j], noise_cutoff):  # if the hit is not noise, we add it to the output file/tree
            wadc[wm] = t.wadc[j]  # output channel number same as input
            wampl[wm] = t.wampl[j]  # output amplitude same as input
            wm[0] += 1  # raising the number of hits (multiplicity) which are not noise
    
    out_t.Fill()  # write an entry with the hits which are not noise into the output tree

""" writing the output file and closing input and output files """
out_file.Write()
out_file.Close()
myfile.Close()
