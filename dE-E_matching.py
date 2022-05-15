# -*- coding: utf-8 -*-
"""
@author: Deni

This is a script for geometrically matching the ADC channels (strips) of different detectors. Primarily front and dE strips, because they are vertical and constitute significant change in angle of detection.  
We are using it prior to determining what hits constitute a particle. 
The idea is to find those pairs that have the highest geometrical overlap and define the area around them. We do this simply by counting and comparing the number of hits that occur simultaneously in those pairs during the same event.

"""

""" import predefined functions """
from ROOT import TFile, TH2S
from time import time
from datetime import timedelta
from sys import argv


""" 
Define program-wide info
    !ADJUST!: N_RUN, TREE_NAME, front_low, front_high, back_low, back_high, dE_low, dE_high
        DEBUG(True,False)
"""
S_RUN = 'run26' if len(argv) < 2 else argv[1]  # S_RUN is the first argument of the scripts, 'run26' if not given
S_RUN = f'{S_RUN}_killnoiseALL'
front = 1 if len(argv) < 3 else int(argv[2])  # S_RUN is the first argument of the scripts, 'run26' if not given
S_SUFFIX = 'dE-E_matching' # suffix to add to the output files and histograms
TREE_NAME = 'tree' #usually we used names "tree" or "T"

DEBUG = True #are we checking the output


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


#lower and upper limit of included strips, upper limit not included in output (for using just one strip, lower limit is "i" and upper "i+1") 
#when we match with the front for the same detector, we usually use all back strips of that detector while matching front and dE strips
for x in range(17, 129, 32):
    if front in range(x, x+16):
        print('Error: detector number is a back detector.')
        exit(0)

if front >= 129:
    print('Error: detector number is a dE detector.')
    exit(0)

back_low = (front//16 + 1) * 16 + 1
back_high = back_low + 16
lst_back = range(back_low, back_high)


"""  
definition of histograms: "histo_0" (ampl,ampl), "histo_1" (ampl,ampl)  
    !ADJUST!: put histogram ranges (X,Y_low,X,Y_high) of interest and choose binning/resolution (Xres,Yres)
"""
X_low, X_high, X_res = 0, 2600, 2600 #front or back strips ranges and resolution
Y_low, Y_high, Y_res = 0, 3000, 3000 #dE strips ranges and resolution

histo_0 = TH2S('(front, dE)', f'ampl-ampl histogram: {S_SUFFIX} {S_RUN}', X_res, X_low, X_high, Y_res, Y_low, Y_high) # 2D histogram records amplitude of a simultaneous hit in front and dE
histo_1 = TH2S('(back, dE)', f'ampl-ampl histogram: {S_SUFFIX} {S_RUN}', X_res, X_low, X_high, Y_res, Y_low, Y_high) # 2D histogram records amplitude of a simultaneous hit in front and dE

""" 
function for saving the histogram in .root form
"""
def write_hists(hist1, hist2):
    f = TFile(f'./matching_data/EF[{front}],dE[{dE}],Eb[{back_low},{back_high-1}]_{S_SUFFIX}_{S_RUN}.root', 'recreate')
    hist1.Write()
    hist2.Write()
    f.Close()


""" EXCECUTION """ 

"""open the file and get the tree"""
myfile = TFile(f'./killnoiseALL/{S_RUN}.root') #open

if DEBUG:
    myfile.ls() #check the file contents

t = myfile.Get(TREE_NAME)  #get the tree

if DEBUG:
    t.Print() #check the tree contents

""" go trough the run data """
n_entries = t.GetEntriesFast() # get the  number of entries in the tree

if DEBUG:
    print("number of entries: ", n_entries)

start_time = time() #start the timer

dE_detectors = { x: 0 for x in range(129, 193) }
for dE in range(129, 193):
    for i in range(n_entries): # for: 0 to n_entries-1
        if i % 50000 == 0: # print the timer info every 50000 entries
            print_time(start_time, (dE-129) * n_entries + i, n_entries*64)
        
        t.GetEntry(i) # get an entry(event)
     
        #make list variables for amplitudes of hits in different detectors
        lst_ampl_front = []
        lst_ampl_back = []
        lst_ampl_dE = []
        
        
        for j in range(t.wm): # go trough all the hits in an entry(event)
            # if t.wadc[j] in lst_front: #fill the list for amplitudes detected in the front strips
            #     lst_ampl_front.append(t.wampl[j])
            if t.wadc[j] == front: #fill the list for amplitudes detected in the front strips
                lst_ampl_front.append(t.wampl[j])

            if t.wadc[j] in lst_back: #fill the list for amplitudes detected in the back strips
                lst_ampl_back.append(t.wampl[j])

            # if t.wadc[j] in lst_dE: #fill the list for amplitudes detected in the dE strips
            #     lst_ampl_dE.append(t.wampl[j])
            if t.wadc[j] == dE: #fill the list for amplitudes detected in the dE strips
                lst_ampl_dE.append(t.wampl[j])
                
        # sort the amplitude by value        
        lst_ampl_front.sort()
        lst_ampl_back.sort() 
        lst_ampl_dE.sort(reverse=True) # reverse the order of dE values since the particles with most total energy should loose the least amound in the dE detector
        # if DEBUG:
            # print(lst_ampl_front, lst_ampl_back, lst_ampl_dE)
        
        # skip further steps if one of the lists is empty, meaning there is no complete information to form particles
        if not lst_ampl_front or not lst_ampl_back or not lst_ampl_dE:
            continue
        
        # skip writing the data if the lists are of different lenghts, meaning there is noise in certain strips or the particle missed some detector
        #there is a possibility that noise "replaces" the missing hit so as to pass through this check, but the proportion of these occurences should be the same for all possible strip matches
        if (len(lst_ampl_dE) != len (lst_ampl_front)) or (len(lst_ampl_dE) > len (lst_ampl_back)): 
            continue
        
        else:  #write the dE-E amplitude pairs into the corresponding histograms
            # if len(lst_ampl_dE) > 2:
                # print(len(lst_ampl_dE))
            for k in range(0, len(lst_ampl_dE)):
                # print('\t', k)
                # print(f'\tfront: {lst_ampl_front[k]}, {lst_ampl_dE[k]}')
                # print(f'\tback: {lst_ampl_front[k]}, {lst_ampl_dE[k]}')
                if lst_ampl_front[k] >= 200:
                    dE_detectors[dE] += 1
                    histo_0.Fill(lst_ampl_front[k], lst_ampl_dE[k]) 
                    histo_1.Fill(lst_ampl_back[k], lst_ampl_dE[k])


""" writing the output file and closing input and output files """
print('Writing amplitude pairs into a file')

keys = sorted(dE_detectors.keys(), key=lambda x: dE_detectors[x], reverse=True)
s = f'{front}\t{keys[0]}'
for i in range(1, 4):
    if dE_detectors[keys[i]] == 0 or dE_detectors[keys[i-1]]/dE_detectors[keys[i]] > 6:
        s += '\n'
        break
    else:
        s += f'\t{keys[i]}'

with open('matching', 'a') as f:
    f.write(s)

print("Writing histograms into a file")
write_hists(histo_0, histo_1)
print("DONE")
myfile.Close()
