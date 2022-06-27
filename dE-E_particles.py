# -*- coding: utf-8 -*-
"""
Created on Tue Jul 16 15:57:49 2019

@author: Deni

This is a script that draws all the desired dE-E graphs after it has been determined which detections constitute particles.
It can draw either all histograms, meaning everything that goes through a certain adc goes into the same graph, or it can draw histogram for specific sets of strips (adcs)
ALL HISTS has a memory overload problem, probably because it attempts to hold too many histograms and histogram entries in memory at the same time. 

"""

""" import predefined functions """
from datetime import timedelta
from sys import argv
from time import time
from ROOT import TFile, TH2S


""" 
Define program-wide info
    !ADJUST!: S_RUN_LIST, S_RUN
        HAS_ENERGY_E (True,False), HAS_ENERGY (True,False), ALFA(True,False), GOLD(True,False), PARTICLE_ENTRY(True,False), EVENT_ENTRY(True,False), 
        REMOVE_CONFLICTS (True,False), REMOVE_INTERSTRIP (True,False), 
        SPECIFIC_HISTS(True,False), STRIP_PAIRS_1(True,False), STRIP_PAIRS_2(True,False)
        STRIP_FRONT, STRIP_dE, strip_front
        DEBUG(True,False)
"""

DEBUG = False


def arg(i, default):
    return default if len(argv) < i+1 else argv[i]

strip_front = int(arg(1, 1))
ALL_PAIRS = bool(arg(2, True))
run_start = arg(3, 'run18')
run_end = arg(4, 'run30')

run_start = int(run_start[3:])
run_end = int(run_end[3:])

S_RUN_LIST = [ f'run{n_run}_TMIN2Be_particles(E=E,F1-4)' for n_run in range(run_start, run_end+1) ]
                
CALIB = 'TMIN2Be' #string to mark the used calibration

OUT_RUN = f'./hist_particles/run{run_start}-{run_end}_TMIN2Be_particles(E=E,F1-4)'
TREE_NAME = 'tree' #usually we used names "tree" or "T"

# do we have energies of the hits, and are only E detectors calibrated (HAS_ENERGY_E = True) or are all detectors calibrated (HAS_ENERGY = True)
HAS_ENERGY_E = False
HAS_ENERGY = True

#are we looking at alfa calibration runs (ALFA = True), or gold scattering calibration runs (GOLD = True), If both are false then we are looking at standard "non-calibration" runs
ALFA = False
GOLD = False

#Define the output tree entry structure
PARTICLE_ENTRY = False  # will an entry in the outgoing file be a single particle or an entire event?

#are we also drawing dE-E graphs without potential interstrip events (REMOVE_INTERSTRIP=True) and/or conflicting events (REMOVE_CONFLICTS=True)
REMOVE_INTERSTRIP = True
REMOVE_CONFLICTS = True

""" Are we drawing all histograms (SPECIFIC_HISTS=False) or some specific subset (SPECIFIC_HISTS=True) that we define by compiling a list of strips or strip pairs (STRIP_PAIRS1,2 = True) that we want to include
 """ 
 
N_ADC = 192 # total number of strips(adcs)

S_PREFIX = 'dE-E(hists)' # a string to mark the output file if we are drawing all histograms

SPECIFIC_HISTS = True

if SPECIFIC_HISTS:
    # STRIP_PAIRS_1 = False #pairs that have the best dE-E match
    # STRIP_PAIRS_2 = True #pairs that have the second best dE-E match
    # ALL_PAIRS = False

    list_pairs = {}
    with open('matching', 'r') as f:
        for l in f.readlines()[1:]:
            l = l.split()
            front = int(l[0])

            if ALL_PAIRS:
                dE = [ int(x) for x in l[1:] ]
            else:
                dE = [ int(l[1]) ]

            list_pairs[front] = dE

    ##lists of used front and dE strips
    strip_front= [strip_front,]
    
    #string mark of used front and dE strips
    STRIP_FRONT = str(strip_front[0])
    STRIP_dE = 'a' if ALL_PAIRS else 's'
                            

"""FUNCTIONS """

""" a function to write histograms into a file when SPECIFIC_HISTS=False """
def write_hists(hists, hists1):
    f = TFile(f'{OUT_RUN}_{S_PREFIX}.root', 'recreate')
    for i in range(N_ADC):
        hists[i].Write()
        if HAS_ENERGY_E: 
            hists1[i].Write()
    f.Close()
    

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
    


""" Defining histograms"""
#defining limits and bin size for "non-calibration" runs
Xd_ampl, Xg_ampl, Xres_ampl = 0, 2600, 2600
Yd_ampl, Yg_ampl, Yres_ampl = 0, 3000, 3000

Xd_energy, Xg_energy, Xres_energy = 0, 50, 10000
Yd_energy, Yg_energy, Yres_energy = 0, 28, 5600

#redefining limits and bin size if we are dealing with calibration runs, alphas or gold scattering
if ALFA:
    Xd_ampl, Xg_ampl, Xres_ampl = 0, 600, 600
    Yd_ampl, Yg_ampl, Yres_ampl = 0, 600, 600
    
    Xd_energy, Xg_energy, Xres_energy = 0, 10, 1000
    Yd_energy, Yg_energy, Yres_energy = 0, 10, 1000
if GOLD: 
    Xd_ampl, Xg_ampl, Xres_ampl = 2000, 3500, 1500
    Yd_ampl, Yg_ampl, Yres_ampl = 0, 800, 800

    Xd_energy, Xg_energy, Xres_energy = 40, 55, 1500
    Yd_energy, Yg_energy, Yres_energy = 0, 10, 1000

#defining the specific set of histograms    
if SPECIFIC_HISTS:    
    h1 = TH2S('dE-E_ampl', "deltaE-E amplitudes", Xres_ampl, Xd_ampl, Xg_ampl, Yres_ampl, Yd_ampl, Yg_ampl) # amplitude dE-E graph
    h1.SetOption("COLZ") # contrast scale
    
    h2 = TH2S('dE-E_energy', "deltaE-E energies", Xres_energy, Xd_energy, Xg_energy, Yres_energy, Yd_energy, Yg_energy) # energy dE-E graf 
    #h2.SetOption("COLZ") # contrast scale
    
    h3 = TH2S('dE-E_energy-ampl', "deltaE-E amplitude vs. energy", Xres_energy, Xd_energy, Xg_energy, Yres_ampl, Yd_ampl, Yg_ampl) # energy E vs amlitude dE 
    #h3.SetOption("COLZ") # contrast scale
    
    if REMOVE_INTERSTRIP:
        h4 = TH2S('dE-E_ampl_NOinterstrip', "deltaE-E amplitudes without potential interstrip events", Xres_ampl, Xd_ampl, Xg_ampl, Yres_ampl, Yd_ampl, Yg_ampl) # amplitude dE-E graph 
        h4.SetOption("COLZ") # contrast scale
        h5 = TH2S('dE-E_energy_NOinterstrip', "deltaE-E energies without potential interstrip events", Xres_energy, Xd_energy, Xg_energy, Yres_energy, Yd_energy, Yg_energy) # energy dE-E graph 
        #h5.SetOption("COLZ") # contrast scale
    if REMOVE_CONFLICTS:
        h6 = TH2S('dE-E_ampl_NOconflicts', "deltaE-E amplitudes without conflicting events", Xres_ampl, Xd_ampl, Xg_ampl, Yres_ampl, Yd_ampl, Yg_ampl) # amplitude dE-E graph 
        h6.SetOption("COLZ") # contrast scale
        h7 = TH2S('dE-E_energy_NOconflicts', "deltaE-E energies without conflicting events", Xres_energy, Xd_energy, Xg_energy, Yres_energy, Yd_energy, Yg_energy) # energy dE-E graph
        #h7.SetOption("COLZ") # contrast scale
#defining all histograms
else:
    histos = [] #field for all emplitude or energy histograms
    histos1 = [] #field for all "mixed" amplitude (dE) energy (E) graphs polje koje će sadržavati histograme miješane amplitude dE i energije E
    for i in range(N_ADC): # for from 0 to N_ADC-1
        name = f"dE-E_ADC={i+1}"  # histogram name
        title = f"Histogram {OUT_RUN}, ADC {i+1}, dE-E"  # histogram title
        
        name1 = f"dE-E_ADC={i+1}_ampl-energy"  # histogram name
        title1 = f"Histogram {OUT_RUN}, ADC {i+1}, dE-E, ampl-energy"  # histogram title
        
        if not HAS_ENERGY:
            h = TH2S(name,title, Xres_ampl, Xd_ampl, Xg_ampl, Yres_ampl, Yd_ampl, Yg_ampl) 
        else:
            h = TH2S(name,title, Xres_energy, Xd_energy, Xg_energy, Yres_energy, Yd_energy, Yg_energy)
        #h.SetOption("COLZ") # contrast scale
        
        #histograms for "mixed" amplitude (dE) energy (E) graphs 
        h1 = TH2S(name1,title1, Xres_energy, Xd_energy, Xg_energy, Yres_ampl, Yd_ampl, Yg_ampl)
        #h1.SetOption("COLZ") #contrast scale
        
        #adding to the field
        histos1.append(h1)
        histos.append(h)
        

""" EXCECUTION """ 
#counters of total number of particles and the number of particles with dE energies below zero (a sign of questionable calibration)
counter_dE_negative = 0
counter_dE_negative_stripPair = 0
counter_particles = 0
counter_particles_stripPair = 0
if REMOVE_INTERSTRIP:
    counter_dE_negative_NOinterstrip = 0
    counter_dE_negative_stripPair_NOinterstrip = 0
    counter_particles_NOinterstrip = 0
    counter_particles_stripPair_NOinterstrip = 0
if REMOVE_CONFLICTS:
    counter_dE_negative_NOconflicts = 0
    counter_dE_negative_stripPair_NOconflicts = 0
    counter_particles_NOconflicts = 0
    counter_particles_stripPair_NOconflicts = 0

#counter of total programme excecution time
start_time = time() #start the excecution timer 

for s_run in S_RUN_LIST:
    myfile = TFile(f'./angles/{s_run}.root') #open
    if DEBUG:
        myfile.ls() #check the file contents

    t = myfile.Get(TREE_NAME)  #get the tree
    if DEBUG:
        t.Print() #check the tree contents

    print(f"reading file {s_run}.root...")
    
    """ go trough the run data """
    n_entries = t.GetEntriesFast() # get the  number of entries in the tree
    if DEBUG:
        print("number of entries: ", n_entries)
        n_entries = 30
        
    start_time_run = time() #start the run timer 

    for i in range(n_entries):
        if i % 50000 == 0: # print(the timer info every 50000 entries
            print_time(start_time_run, i, n_entries)
        
        t.GetEntry(i) # get an entry
        
        if PARTICLE_ENTRY: #if particles are entries
            #counting particles
            counter_particles += 1 
            if REMOVE_INTERSTRIP and t.interstrip == 0: counter_particles_NOinterstrip += 1
            if REMOVE_CONFLICTS and t.conflicts == 0: counter_particles_NOconflicts += 1

            
            if HAS_ENERGY:
                #counting particles with dE energies below zero
                if t.nrg[2] <0: 
                    counter_dE_negative+=1
                    if REMOVE_INTERSTRIP and t.interstrip == 0: counter_dE_negative_NOinterstrip += 1                     
                    if REMOVE_CONFLICTS and t.conflicts == 0: counter_dE_negative_NOconflicts += 1
   
            if SPECIFIC_HISTS:
                if  t.adc[0] in strip_front and t.adc[2] in list_pairs[t.adc[0]]:
                    #counting particels and filling the appropriate histograms
                    counter_particles_stripPair += 1
                    h1.Fill(t.ampl[0],t.ampl[2]) #amplitude dE-E
                    
                    if REMOVE_INTERSTRIP and t.interstrip == 0:
                        counter_particles_stripPair_NOinterstrip += 1
                        h4.Fill(t.ampl[0],t.ampl[2])
                    if REMOVE_CONFLICTS and t.conflicts == 0:
                        counter_particles_stripPair_NOconflicts += 1
                        h6.Fill(t.ampl[0],t.ampl[2])    
                        
                    if HAS_ENERGY:
                        
                        if t.nrg[2] <0:  #counting particles with dE energies below zero
                            counter_dE_negative_stripPair +=1
                            if REMOVE_INTERSTRIP and t.interstrip == 0: counter_dE_negative_stripPair_NOinterstrip += 1 
                            if REMOVE_CONFLICTS and t.conflicts == 0: counter_dE_negative_stripPair_NOconflicts += 1
                               
                        h2.Fill(t.nrg[0],t.nrg[2]) #energies dE-E
                        if REMOVE_INTERSTRIP and t.interstrip == 0: h5.Fill(t.nrg[0],t.nrg[2])
                        if REMOVE_CONFLICTS and t.conflicts == 0: h7.Fill(t.nrg[0],t.nrg[2])
                            
                    if HAS_ENERGY_E:# energy E vs amplitude dE
                        h3.Fill(t.nrg[0],t.ampl[2])
            else: 
                # fill front vs dE energy and amplitude histograms from the perspective of front strips
                if not HAS_ENERGY:
                    histos[t.adc[0]-1].Fill(t.ampl[0],t.ampl[2])
                else:
                    histos[t.adc[0]-1].Fill(t.nrg[0],t.nrg[2])
                if HAS_ENERGY_E: 
                    histos1[t.adc[0]-1].Fill(t.nrg[0],t.ampl[2])
            
                # fill back vs deltaE energy and amplitude histograms
                if not HAS_ENERGY:
                    histos[t.adc[1]-1].Fill(t.ampl[1],t.ampl[2])
                else:
                    histos[t.adc[1]-1].Fill(t.nrg[1],t.nrg[2])
                if HAS_ENERGY_E: 
                    histos1[t.adc[1]-1].Fill(t.nrg[1],t.ampl[2])
            
                # fill front vs deltaE histograms from the perspective of dE strips
                if not HAS_ENERGY:
                    histos[t.adc[2]-1].Fill(t.ampl[0],t.ampl[2])
                else:
                    histos[t.adc[2]-1].Fill(t.nrg[0],t.nrg[2])
                if HAS_ENERGY_E: 
                    histos1[t.adc[2]-1].Fill(t.nrg[0],t.ampl[2])
        else:  # if events are entries
            #counting particles
            counter_particles += t.cnc
            if REMOVE_INTERSTRIP and t.interstrip == 0: counter_particles_NOinterstrip += t.cnc
            if REMOVE_CONFLICTS and t.conflicts == 0: counter_particles_NOconflicts += t.cnc
                
            for i in range(t.cnc): #going through all the particles in an event
                #counters to go through front, back and dE data positions 
                front_counter = i*3 
                back_counter = i*3+1
                dE_counter = i*3+2
                
                if HAS_ENERGY:
                    if t.nrg[dE_counter] <0: #counting particles with dE energies below zero
                        counter_dE_negative+=1
                        if REMOVE_INTERSTRIP and t.interstrip == 0: counter_dE_negative_NOinterstrip += 1                     
                        if REMOVE_CONFLICTS and t.conflicts == 0: counter_dE_negative_NOconflicts += 1
                
                if SPECIFIC_HISTS:
                    if t.adc[front_counter] in strip_front and t.adc[front_counter] in list_pairs.keys() and t.adc[dE_counter] in list_pairs[t.adc[front_counter]]:
                        #counting particels and filling the appropriate histograms
                        counter_particles_stripPair += 1
                        h1.Fill(t.ampl[front_counter],t.ampl[dE_counter]) #amplitude dE-E
                        
                        if REMOVE_INTERSTRIP and t.interstrip == 0:
                            counter_particles_stripPair_NOinterstrip += 1
                            h4.Fill(t.ampl[front_counter],t.ampl[dE_counter])
                        if REMOVE_CONFLICTS and t.conflicts == 0:
                            counter_particles_stripPair_NOconflicts += 1
                            h6.Fill(t.ampl[front_counter],t.ampl[dE_counter]) 
                            
                        if HAS_ENERGY:
                            if t.nrg[dE_counter] <0: #counting particles with dE energies below zero
                                counter_dE_negative_stripPair +=1
                                if REMOVE_INTERSTRIP and t.interstrip == 0: counter_dE_negative_stripPair_NOinterstrip += 1 
                                if REMOVE_CONFLICTS and t.conflicts == 0: counter_dE_negative_stripPair_NOconflicts += 1
                                
                            h2.Fill(t.nrg[front_counter],t.nrg[dE_counter]) #energies dE-E
                            if REMOVE_INTERSTRIP and t.interstrip == 0: h5.Fill(t.nrg[front_counter],t.nrg[dE_counter])
                            if REMOVE_CONFLICTS and t.conflicts == 0: h7.Fill(t.nrg[front_counter],t.nrg[dE_counter])
                            
                        if HAS_ENERGY_E: # energy E vs amplitude dE
                            h3.Fill(t.nrg[front_counter],t.ampl[dE_counter])
                else: 
                    # fill front vs dE energy and amplitude histograms from the perspective of front strips
                    if not HAS_ENERGY:
                        histos[t.adc[front_counter]-1].Fill(t.ampl[front_counter],t.ampl[dE_counter])
                    else:
                        histos[t.adc[front_counter]-1].Fill(t.nrg[front_counter],t.nrg[dE_counter])
                    if HAS_ENERGY_E: 
                        histos1[t.adc[front_counter]-1].Fill(t.nrg[front_counter],t.ampl[dE_counter])
                
                    # fill back vs deltaE energy and amplitude histograms
                    if not HAS_ENERGY:
                        histos[t.adc[back_counter]-1].Fill(t.ampl[back_counter],t.ampl[dE_counter])
                    else:
                        histos[t.adc[back_counter]-1].Fill(t.nrg[back_counter],t.nrg[dE_counter])
                    if HAS_ENERGY_E: 
                        histos1[t.adc[back_counter]-1].Fill(t.nrg[back_counter],t.ampl[dE_counter])
                
                    # fill front vs deltaE histograms from the perspective of dE strips
                    if not HAS_ENERGY:
                        histos[t.adc[dE_counter]-1].Fill(t.ampl[front_counter],t.ampl[dE_counter])
                    else:
                        histos[t.adc[dE_counter]-1].Fill(t.nrg[front_counter],t.nrg[dE_counter])
                    if HAS_ENERGY_E: 
                        histos1[t.adc[dE_counter]-1].Fill(t.nrg[front_counter],t.ampl[dE_counter])
                    

                            
#when there are multiple runs, we calculate the total elapsed time
total_time = round(time()-start_time)
total_time = timedelta(0,total_time)                                            

print("Writing histograms into a file")
if SPECIFIC_HISTS:
    g = TFile(f"{OUT_RUN}_{STRIP_FRONT}{STRIP_dE}.root", "recreate")
    
    h1.Write()
    if REMOVE_INTERSTRIP:
        h4.Write()
    if REMOVE_CONFLICTS:
        h6.Write()
    
    if HAS_ENERGY: 
        h2.Write()
        if REMOVE_INTERSTRIP:
            h5.Write()
        if REMOVE_CONFLICTS:
            h7.Write()
        
    if HAS_ENERGY_E:
        h3.Write()
        
    g.Close()
else:
    write_hists(histos, histos1)  
    
print(OUT_RUN)

if SPECIFIC_HISTS:
    print(f"Chosen strip pairs ({STRIP_FRONT}, {STRIP_dE})")

print("Total number of particles in all detectors:", counter_particles)
print("Total number of particles with dE energy below zero in all detectors:", counter_dE_negative)

if SPECIFIC_HISTS:
    print("The number of particles in the chosen strip pairs:", counter_particles_stripPair)
    print("The number of particles with energies below zero in chosen strip pairs:", counter_dE_negative_stripPair)

if REMOVE_INTERSTRIP:
    print("WITHOUT POTENTIAL INSTERSTRIP EVENTS:")
    print("\tTotal number of particles in all detectors:", counter_particles_NOinterstrip)
    print("\tTotal number of particles with dE energy below zero in all detectors:", counter_dE_negative_NOinterstrip)
    if SPECIFIC_HISTS:
        print("\tThe number of particles in the chosen strip pairs:", counter_particles_stripPair_NOinterstrip)
        print("\tThe number of particles with energies below zero in chosen strip pairs:", counter_dE_negative_stripPair_NOinterstrip)

if REMOVE_CONFLICTS:
    print("WITHOUT CONFLICTING EVENTS:")
    print("\tTotal number of particles in all detectors:", counter_particles_NOconflicts)
    print("\tTotal number of particles with dE energy below zero in all detectors:", counter_dE_negative_NOconflicts)
    if SPECIFIC_HISTS:
        print("\tThe number of particles in the chosen strip pairs:", counter_particles_stripPair_NOconflicts)
        print("\tThe number of particles with energies below zero in chosen strip pairs:", counter_dE_negative_stripPair_NOconflicts)

print("The total time of excecution", total_time)
print("END")
