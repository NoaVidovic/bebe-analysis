# -*- coding: utf-8 -*-
"""
Created on Thu Jun  2 02:38:33 2022

@author: Deni

This is a script that takes the graphical cuts made on dE-E graphs and uses them to name the detected particles. 

"""

""" import predefined functions """
from datetime import timedelta
from numpy import zeros
from sys import argv
from time import time
from ROOT import TFile, TTree, TCutG


""" 
Define program-wide info
    !ADJUST!: N_RUN, CALIB, PARTICLES, PTYPE
        HAS_ENERGY (True,False), PARTICLE_ENTRY(True,False), EVENT_ENTRY(True,False), INTERSTRIP_MARKING (True,False), CONFLICT_MARKING (True,False) 
        particles that we have cuts for; He4(True,False), He6(True,False) ...
        DEBUG (True,False)
"""

DEBUG = False # do we want an output check of certain parts of the code

S_RUN = 'run26' if len(argv) < 2 else argv[1]  # S_RUN is the first argument of the scripts, 'run26' if not given

CALIB = 'RUTH90' # the name of the calibration;
PARTICLES = 'particles(E=E,F1-4)'#info about the particle finding process
S_RUN = f'{S_RUN}_{CALIB}_{PARTICLES}'
TREE_NAME = 'tree' #usually we used names "tree" or "T"


HAS_ENERGY = True # do we already have the energies of the hits, or do we only have amplitudes

MIN_mult = 3 # minimum number of hits in the event needed to define a physical event, dE, FRONT and BACK info needed to constitue a single particle


#Define the output tree entry structure
PARTICLE_ENTRY = False  # will an entry in the outgoing file be a single particle or an entire event?

#maximum number of particles that we write into a single root entry
if PARTICLE_ENTRY: 
    MAX_PARTICLES_ENTRY = 1 
else:
    MAX_PARTICLES_ENTRY = 30 # any number that is higher than the expected maximum number of particles in an experimental event, 10 is enough when the conflict filter is on, more if it is not
    
CONFLICT_MARKING = True # do we have conflicted particles who share the same adc strip marked
INTERSTRIP_MARKING = True  # do we have potentital interstrip particles marked

#all the particles that we are naming i.e. that we have the graphical cuts for
He4 = True
He6 = True
Li6 = True
Li7 = True
Li8 = True
Be9 = True
Be10 = True

particles = []
markings = {}
included_particles = ''

if He4:
    particles.append('4He')
    markings['4He'] = 204
if He6:
    particles.append('6He')
    markings['6He'] = 206
if Li6:
    particles.append('6Li')
    markings['6Li'] = 306
if Li7:
    particles.append('7Li')
    markings['7Li'] = 307
if Li8:
    particles.append('8Li')
    markings['8Li'] = 308
if Be9:
    particles.append('9Be')
    markings['9Be'] = 409
if Be10:
    particles.append('10Be')
    markings['10Be'] = 410
    
#do we also have cuts for the second best geometrical matches of dE and E
MATCHES = 's'
    
# do we have cuts made on energy or amplitude graphs
cutAMPL = True  # True if cut on AMPL, False if cut on NRG

# the name of the output file 
PTYPE = f'ptype({",".join(particles)})'
OUT_RUN = f'{S_RUN}_{PTYPE}_CUT'


""" CUTS manipulation 
direct variable definition for all cuts as a first option --> different naming and dictionary structure possible (example below)
ideally we would just "get" cuts inside the particle loops, but this creates a memory overload problem --> try to clean the memory after each loop, an idea given in that part of the code

""" 
# 1. open the cut file
f = TFile('./hist_particles/cuts.root') #master file with cuts for all particles and strip pairs

# 2. get the cuts
front_list = sum([ list(range(x, x+16)) for x in range(1, 100, 32) ], [])

# outer dictionary: first key is particle type, pointing to another dictionary
# inner dictionary: key is front strip, pointing to the appropriate cut
cuts = { p: { i: f.Get(f'{i}{MATCHES}_{p}') for i in front_list } for p in particles }
for p in particles:
    for i in front_list:
        if type(cuts[p][i]) is not TCutG:
            print(f'{i}{MATCHES}_{p}: MISSING')
f.Close()

#pairs of front and dE strips with best and 2nd best geometrical match             
list_pairs = {}
with open('matching', 'r') as f:
    for l in f.readlines()[1:]:
        l = l.split()
        front = int(l[0])

        if MATCHES == 'a':
            dE = [ int(x) for x in l[1:] ]
        elif MATCHES == 's':
            dE = [ int(l[1]) ]
        else:
            print('INVALID MATCHES')
            exit(0)

        list_pairs[front] = dE


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
    
    
    
""" EXCECUTION """ 
       
"""open the file and get the tree"""
myfile = TFile(f'./angles/{S_RUN}.root') #open
if DEBUG:
    myfile.ls() #check the file contents

t = myfile.Get(TREE_NAME)  #get the tree
if DEBUG:
    t.Print() #check the tree contents
    
"""creating the output root file; defining name, output tree, output tree variables, output tree structure """
out_file = TFile(f'./ptype/{OUT_RUN}.root', 'recreate') # output root file
out_t = TTree('tree', f'{OUT_RUN}.root') # output tree


# output tree variables
# variables that exist in the starting tree
mult = zeros(1, dtype='short') # multiplicity of an event - number of hits, value of type Short
adc = zeros(MIN_mult * MAX_PARTICLES_ENTRY, dtype='short') # channel number, field of values of type Short
ampl = zeros(MIN_mult * MAX_PARTICLES_ENTRY, dtype='short') # amplitude, field of values of type Short
if HAS_ENERGY:
    nrg = zeros(MIN_mult * MAX_PARTICLES_ENTRY, dtype='float') # energy, field of values of type Float
event = zeros(1, dtype= 'int')# the number of events in the starting run, value of type Int
detector = zeros(MAX_PARTICLES_ENTRY, dtype='short')# which detector the particle went into; 1,2,3 or 4, value of type Short
cnc = zeros(1, dtype='short') # the number of particles in an event (coincidences), value of type Short
if CONFLICT_MARKING:
    conflicts = zeros(1, dtype='short') # the number of conflicting particles in an event, value of type Short
if INTERSTRIP_MARKING:
    interstrip = zeros(1, dtype='short') # the number of potential interstrip particles in an event, value of type Short         
r = zeros(MAX_PARTICLES_ENTRY, dtype='float') # distance from the target to the pixel where the particle is detected, field of values of type Float
theta = zeros(MAX_PARTICLES_ENTRY, dtype='float') # angle from the target to the pixel where the particle is detected, field of values of type Float
phi = zeros(MAX_PARTICLES_ENTRY, dtype='float') # angle from the target to the pixel where the particle is detected, field of values of type Float
#new variables  
ptype = zeros(MAX_PARTICLES_ENTRY, dtype='short')# marking of the particle type, value of type Short

out_t.Branch('event', event, 'event/I') # event number branch
out_t.Branch('cnc', cnc, 'cnc/S') # coincidences branch    
out_t.Branch('mult',mult,'mult/S') # multiplicity branch  
if PARTICLE_ENTRY:  
    out_t.Branch('adc',adc,'adc[3]/S') # channel number branch, field of lenght 3
    out_t.Branch('ampl',ampl,'ampl[3]/S') # amplitude branch, field of lenght 3
    if  HAS_ENERGY:
        out_t.Branch('nrg', nrg, 'nrg[3]/D') # energy branch, field of lenght 3
    out_t.Branch('detector', detector, 'detector/S') # detector branch
    out_t.Branch('r', r, 'r/D') # particle distance branch
    out_t.Branch('theta', theta, 'theta/D') # particle angle branch
    out_t.Branch('phi', phi, 'phi/D') # particle angle branch
    out_t.Branch('ptype', ptype, 'ptype/S') # particle type branch
else:
    out_t.Branch('adc',adc,'adc[mult]/S') # channel number branch, field of lenght "mult"
    out_t.Branch('ampl',ampl,'ampl[mult]/S') # amplitude branch, field of lenght "mult"
    if  HAS_ENERGY:
        out_t.Branch('nrg', nrg, 'nrg[mult]/D') # energy branch, field of lenght "mult"
    out_t.Branch('detector', detector, 'detector[cnc]/S') # detector branch, field of lenght "cnc"  
    out_t.Branch('r', r, 'r[cnc]/D') # particle distance branch
    out_t.Branch('theta', theta, 'theta[cnc]/D') # particle angle branch
    out_t.Branch('phi', phi, 'phi[cnc]/D') # particle angle branch
    out_t.Branch('ptype', ptype, 'ptype[cnc]/S') # particle type branch, field of lenght "cnc"  

if CONFLICT_MARKING:
    out_t.Branch('conflicts', conflicts, 'conflicts/S') #conflict branch

if INTERSTRIP_MARKING:
    out_t.Branch('interstrip', interstrip, 'interstrip/S') #interstrip branch    
    
    
""" go trough the run data """
n_entries = t.GetEntriesFast() # get the  number of entries in the tree
if DEBUG:
    print("number of entries: ", n_entries)
    n_entries = 30
    
start_time = time() #start the timer 



for i in range(n_entries):
    if i % 50000 == 0: # print the timer info every 50000 entries
        print_time(start_time,i,n_entries)
        
    t.GetEntry(i) # get an entry
    
    # variables that are attributed to all of the particles in an event: event, mult, cnc, conflicts, interstrip
    event[0] = t.event
    mult[0] = t.mult
    cnc[0] = t.cnc
    if CONFLICT_MARKING:
        conflicts[0] = t.conflicts
    if INTERSTRIP_MARKING:
        interstrip[0] = t.interstrip
        
        # variables that are different for different particles in an event: detector, adc, ampl, nrg, r, theta, phi, ptype   
    if PARTICLE_ENTRY:
        front_counter = 0 #counter that goes though field positions that contain information about the front strip
        dE_counter = 2 #counter that goes through field positions that contain information about the dE strip
        
        detector[0] = t.detector        
        r[0] = t.r
        theta[0] = t.theta
        phi[0] = t.phi
        #copying adc, ampl and energy info from the previous run          
        for i in range(3):
            adc[i] = t.adc[i]
            ampl[i] = t.ampl[i]
            if HAS_ENERGY:
                nrg[i] = t.nrg[i]
                
        #adding the particle type info
        ptype[0] = 0     
        #if particles belong to the matched strip pairs
        if t.adc[dE_counter] in list_pairs[t.adc[front_counter]]:
            if cutAMPL:
                x, y = t.ampl[front_counter], t.ampl[dE_counter] 
            else:
                x, y = t.nrg[front_counter], t.nrg[dE_counter]

            for p in particles:
                c = cuts[p][t.adc[front_counter]]
                if type(c) is TCutG:
                    if c.IsInside(x, y):
                        ptype[0] = markings[p]
                        break
                else:
                    pass
                    # print(p, t.adc[front_counter])
                
                
        out_t.Fill() # filling the out tree; particle is the output entry
    else:
        #going through all the particles in an event
        for i in range(cnc[0]):
            front_counter = i*3 #counter that goes though field positions that contain information about the front strip
            dE_counter = i*3+2 #counter that goes through field positions that contain information about the dE strip
            front = t.adc[front_counter] #strip(adc) front number
            dE = t.adc[dE_counter] #strip(adc) dE number
            
            detector[i] = t.detector[i]
            r[i] = t.r[i]
            theta[i] = t.theta[i]
            phi[i] = t.phi[i]
            
            #adding the particle type info
            ptype[i] = 0     
            
            #do particles belong to the matched strip pairs
            if t.adc[dE_counter] in list_pairs[t.adc[front_counter]]:
                if cutAMPL: 
                    x, y = t.ampl[front_counter], t.ampl[dE_counter] 
                else: 
                    x, y = t.nrg[front_counter], t.nrg[dE_counter]

                for p in particles:
                    c = cuts[p][t.adc[front_counter]]
                    if type(c) is TCutG:
                        if c.IsInside(x, y):
                            ptype[i] = markings[p]
                            break
                    else:
                        pass
                        # print(p, t.adc[front_counter])
        
        #copying adc, ampl and energy info from the previous run            
        for j in range(mult[0]):
            adc[j] = t.adc[j]
            ampl[j] = t.ampl[j]
            if HAS_ENERGY:
                nrg[j] = t.nrg[j]   

                        
        out_t.Fill() # filling the out tree; event is the output entry as in the raw root files
        
out_file.Write()
out_file.Close()
myfile.Close()
