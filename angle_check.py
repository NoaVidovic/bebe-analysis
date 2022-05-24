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
from ROOT import TFile, TH1F


DEBUG = False


run_start = 'run18' if len(argv) < 2 else argv[1]
run_end = 'run30' if len(argv) < 3 else argv[2]

run_start = int(run_start[3:])
run_end = int(run_end[3:])

S_RUN_LIST = [ f'run{n_run}_TMIN2Be_particles(E=E,F1-4)' for n_run in range(run_start, run_end+1) ]
                
CALIB = 'TMIN2Be' #string to mark the used calibration
REACTION = '9Be+9Be'  #projectile and target

OUT_RUN = f'./angle_check/angle_check_{REACTION}_run{run_start}-{run_end}_TMIN2Be_particles(E=E,F1-4).root'



""" defining detectors used """
det1 = 1
det2 = 3
det3 = 4



""" Defining histograms""" 
h12 = TH1F('theta1+theta2_AB','theta1+theta2 in two particle detections, detectors A and B', 1000, 30, 130)
h13 = TH1F('theta1+theta2_AC','theta1+theta2 in two particle detections, detectors A and C', 1000, 30, 130)
h14 = TH1F('theta1+theta2_AD','theta1+theta2 in two particle detections, detectors A and D', 1000, 30, 130)
h23 = TH1F('theta1+theta2_BC','theta1+theta2 in two particle detections, detectors B and C', 1000, 30, 130)
h24 = TH1F('theta1+theta2_BD','theta1+theta2 in two particle detections, detectors B and D', 1000, 30, 130)
h34 = TH1F('theta1+theta2_CD','theta1+theta2 in two particle detections, detectors C and D', 1000, 30, 130)



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
    myfile = TFile(f'./angles/{s_run}.root')
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
            if t.detector[0] == 1:
                if t.detector[1] == 2:
                    h12.Fill(t.theta[0] + t.theta[1])
                elif t.detector[1] == 3:
                    h13.Fill(t.theta[0] + t.theta[1])
                elif t.detector[1] == 4:
                    h14.Fill(t.theta[0] + t.theta[1])
            elif t.detector[0] == 2:
                if t.detector[1] == 3:
                    h23.Fill(t.theta[0] + t.theta[1])
                elif t.detector[1] == 4:
                    h24.Fill(t.theta[0] + t.theta[1])
            elif t.detector[0] == 3:
                if t.detector[1] == 4:
                    h34.Fill(t.theta[0] + t.theta[1])
                    
f = TFile(OUT_RUN, 'recreate')  
h12.Write()
h13.Write()
h14.Write()
h23.Write()
h24.Write()
h34.Write()
f.Close() 
