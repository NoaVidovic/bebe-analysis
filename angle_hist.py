# -*- coding: utf-8 -*-
"""
Created on Tue May 24 06:54:49 2022

@author: pmf1

This is a script used to check the quality of angle calibration. 

Detection angles of elastically scatered 9Be+9Be should give 90 degrees when added since the two particles have the same mass.
Cut on 9Be particles should not be necessary since the elastic scattering should be dominant over other reaction channels.
Combinations of detectors that should have 90 degrees when added: A i D, A i C


"""

from datetime import timedelta
from sys import argv
from time import time
from ROOT import TFile, TH2S


DEBUG = False


run_start = 'run18' if len(argv) < 2 else argv[1]
run_end = 'run30' if len(argv) < 3 else argv[2]

run_start = int(run_start[3:])
run_end = int(run_end[3:])

S_RUN_LIST = [ f'run{n_run}_TMIN2Be_particles(E=E,F1-4)_ptype(4He,6He,6Li,7Li,8Li,9Be,10Be)_CUT' for n_run in range(run_start, run_end+1) ]
                
CALIB = 'TMIN2Be' # string to mark the used calibration
REACTION = '9Be+9Be'  # projectile and target
PARTICLES = '9Be9Be'  # particles that the script filters for

OUT_RUN = f'./angle_check/angle_hist_{REACTION}_run{run_start}-{run_end}_TMIN2Be_particles(E=E,F1-4)_{PARTICLES}.root'




""" Defining histograms""" 
res_t, low_t, high_t = 550, 10, 65  # theta1 and theta2
res_s, low_s, high_s = 600, 60, 120  # theta1 + theta2

h1s = TH2S('theta1-theta1+theta2', "theta1 - theta1+theta2 in 9Be elastic collisions, detectors A and D", res_s, low_s, high_s, res_t, low_t, high_t)
h2s = TH2S('theta2-theta1+theta2', "theta2 - theta1+theta2 in 9Be elastic collisions, detectors A and D", res_s, low_s, high_s, res_t, low_t, high_t)
h12 = TH2S('theta1-theta2', "theta1 - theta2 in 9Be elastic collisions, detectors A and D", res_t, low_t, high_t, res_t, low_t, high_t)



"""FUNCTIONS """

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
    
    
""" EXECUTION"""


for s_run in S_RUN_LIST:
    # open the file
    myfile = TFile(f'./ptype/{s_run}.root')
    #myfile.ls()
    t = myfile.Get('tree')
    print(f"reading file {s_run}.root...")
    
    
    n_entries = t.GetEntriesFast() # get the  number of entries in the tree
    if DEBUG:
        print("number of entries: ", n_entries)
        n_entries = 30
        
    start_time = time() #start the timer 


    for i in range(n_entries): 
        if i % 50000 == 0: # print the timer info every 50000 entries
            print_time(start_time,i,n_entries)
            
        t.GetEntry(i) # get an entry
        
        if t.cnc == 2: #fill with coincident detections
            if t.ptype[0] == 409 and t.ptype[1] == 409 and t.detector[0] == 1 and t.detector[1] == 4:
                h1s.Fill(t.theta[0]+t.theta[1], t.theta[0])
                h2s.Fill(t.theta[0]+t.theta[1], t.theta[1])
                h12.Fill(t.theta[0], t.theta[1])
                    
f = TFile(OUT_RUN, 'recreate')  
h1s.Write()
h2s.Write()
h12.Write()
f.Close() 
