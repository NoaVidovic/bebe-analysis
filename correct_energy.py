# -*- coding: utf-8 -*-
"""
Created on Mon Jun  6 03:28:39 2022

@author: Deni

This is a script that adds up all energy contributions of a particle; energies detected at different detector parts and the energy losses in the target and the dead layer.
This forms the total particle energy in the reaction which then enables the production of excitation spectra. 
It is necessary to determine particle types, angles, and detected energies prior to this step. 
"""


""" import predefined functions """
import numpy as np

from datetime import timedelta
from numpy import zeros
from sys import argv
from time import time
from ROOT import TFile,TTree

""" 
Define program-wide info
    !ADJUST!: N_RUN, CALIB, PARTICLES, CUT, TARGET
        PARTICLE_ENTRY(True,False), INTERSTRIP_MARKING (True,False), CONFLICT_MARKING (True,False) 
        DEBUG (True,False)
"""

DEBUG = False # do we want an output check of certain parts of the code

S_RUN = 'run26' if len(argv) < 2 else argv[1]  # S_RUN is the first argument of the scripts, 'run26' if not given

CALIB = 'RUTH90'  # the name of the calibration
PARTICLES = 'particles(E=E,F1-4)'  # info about the particle finding process
PTYPE = 'ptype(4He,6He,6Li,7Li,8Li,9Be,10Be)'  # the particles for which the type is determined
S_RUN = f'{S_RUN}_{CALIB}_{PARTICLES}_{PTYPE}_CUT'
TREE_NAME = 'tree' #usually we used names "tree" or "T"


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

# the name of the output file 
SUFFIX_ENERGY = 'cnrg'
OUT_RUN = f'{S_RUN}_cnrg'


"""
Target and assumed aluminium dead layer info.
The dead layer width and half of the target width.
Correction coefficients that come from fiting a function in the form of the below defined "fit_function" to energy loss data coming from SRIM software.
Since we are reconstructing energies from the final point of the particle path (E detector) to the start od the reaction, we appropriately modify the SRIM data (exit point energy).
"""
#debljina mrtvog sloja i korak debljine pri računu korekcija
DLAYER = 135 #dead layer width in ug/cm2
CORRECTION_STEP = 5 #the step of correction calculation in ug/cm2

TARGET = '9Be' #the target used
HALF_TARGET = {'9Be': 239.5} #half of the target width in ug/cm2 


## -------------------- Be target, correction coefficients
CORRECT_TARGET_HIGH = {
    204: {'element': 'He4(2-55MeV)', 'array': [  0.34959465,  13.44789442,   1.46373525,   0.07630944,  -8.28262284]},
    206: {'element': 'He6(2-55MeV)', 'array': [  0.49668402,  16.12320631,   2.11048743,   0.06656297, -11.84410601]},
    306: {'element': 'Li6(3-55MeV)', 'array': [  4.24865824e+02,   7.13420165e+01,  -4.23834328e+02, -2.18177963e-05,  -9.67396726e+01]},
    307: {'element': 'Li7(3-55MeV)', 'array': [  1.15461707,  53.73517537,   4.84491891,   0.0775495 , -72.94684789]},  
    308: {'element': 'Li8(3-55MeV)', 'array': [  1.50529724,  46.43293467,   7.1348864 ,   0.07612878, -65.82766616]},       
    409: {'element': 'Be9(3-55MeV)', 'array': [  2.92034060e+00,   4.57408450e+01,   1.47508907e+01, 5.17633004e-02,  -7.42552211e+01]},
    410: {'element': 'Be10(3-55MeV)', 'array':[  3.21548235e+00,   4.30014007e+01,   1.57264651e+01, 4.94028785e-02,  -7.55787471e+01]},
    }
CORRECT_TARGET_LOW = {
    204: {'element': 'He4(0.3-2MeV)', 'array': [  4.28132349,  -1.46755886,  10.29661932,   0.66750658,   0.08509652]},
    206: {'element': 'He6(0.3-2MeV)', 'array': [  4.21186404,  -2.19243545,  10.34552502,   0.43943333,   0.1873723 ]},
    306: {'element': 'Li6(0.3-3MeV)', 'array': [ 10.50338577,  11.75361073, -61.50967052,   2.67638404,  -1.0391496 ]},
    307: {'element': 'Li7(0.3-3MeV)', 'array': [ 10.94622547,  12.52639622, -60.56257845,   2.36612142,  -1.18774464]},
    308: {'element': 'Li8(0.3-3MeV)', 'array': [ 11.2957328 ,  13.30323643, -59.53632644,   2.11711134,  -1.35582625]},       
    409: {'element': 'Be9(0.3-3MeV)', 'array': [ 21.17065791,   5.55437781, -32.2261036 ,   1.6966735 ,  -0.67590671]},
    410: {'element': 'Be10(0.3-3MeV)', 'array': [ 21.89886047,   3.98153361, -27.59301738,   1.54787644,  -0.50455359]},
    }
    
LIMITS_TARGET = {
    204: {'element': 'He4', 'low':  0.3, 'high': 2.0},
    206: {'element': 'He6', 'low':  0.3, 'high': 2.0},
    306: {'element': 'Li6', 'low':  0.3, 'high': 3.0},
    307: {'element': 'Li7', 'low':  0.3, 'high': 3.0},
    308: {'element': 'Li8', 'low':  0.3, 'high': 3},
    409: {'element': 'Be9', 'low':  0.3, 'high': 3.0},
    410: {'element': 'Be10', 'low':  0.3, 'high': 3.0},
    }


## -------------------- Aluminium Dlayer, correction coefficients
CORRECT_DLAYER_HIGH = {
    204: {'element': 'He4(3-55MeV)', 'array': [  0.27448595,  12.82525106,   0.81438843,   0.05580138, -10.92251642]},
    206: {'element': 'He6(3-55MeV)', 'array': [  0.4163073 ,  14.46910447,   1.36629708,   0.05492109, -13.36339859]},
    306: {'element': 'Li6(3-55MeV)', 'array': [  1.07028691,  29.20290404,   4.03405203,   0.06433485, -40.95640997]},
    307: {'element': 'Li7(3-55MeV)', 'array': [  1.30011295,  25.73037226,   5.20909122,   0.06301463, -38.93595838]},
    308: {'element': 'Li8(3-55MeV)', 'array': [  1.507279  ,  22.27482037,   6.21500898,   0.0609492 , -36.62119809]},       
    409: {'element': 'Be9(3-55MeV)', 'array': [  2.42477738,  16.39398731,  10.94076062,   0.03668005, -23.21906396]},
    410: {'element': 'Be10(3-55MeV)', 'array': [  2.48818972,  17.61240673,  10.97353712,   0.03352156, -27.38524116]},
    }
    
CORRECT_DLAYER_LOW = {
    204: {'element': 'He4(0.3-3MeV)', 'array': [ 3.02517364, -1.06032362,  7.42345954,  0.5635764 ,  0.05191343]},
    206: {'element': 'He6 (0.3-3MeV)', 'array': [ 3.25266764, -1.80469287,  7.55412956,  0.40370607,  0.16113877]},
    306: {'element': 'Li6(0.3-3MeV)', 'array': [ 13.00104136,  -2.02195511,  -0.82158575,  -0.4026955 ,   0.19375556]},
    307: {'element': 'Li7(0.3-3MeV)', 'array': [ 12.04955631,  -2.1020286 ,  -0.22140863,  -0.63512487,   0.20655154]},
    308: {'element': 'Li8(0.3-3MeV)', 'array': [ 11.6509913 ,  -2.21507943,  -0.05474745,  -0.92223369,   0.22611104]},       
    409: {'element': 'Be9(0.3-3MeV)', 'array': [ 13.95325397,   4.12291305, -25.19504877,   1.87065402,  -0.46153586]},
    410: {'element': 'Be10(0.3-3MeV)', 'array': [ 14.36647877,   3.345535  , -22.75844531,   1.7263639 ,  -0.3751475 ]},
    }
    
LIMITS_DLAYER = {
    204: {'element': 'He4', 'low':  0.3, 'high': 3.0},
    206: {'element': 'He6', 'low':  0.3, 'high': 3.0},
    306: {'element': 'Li6', 'low':  0.3, 'high': 3.0},
    307: {'element': 'Li7', 'low':  0.3, 'high': 3.0},
    308: {'element': 'Li8', 'low':  0.3, 'high': 3.0},
    409: {'element': 'Be9', 'low':  0.3, 'high': 3.0},
    410: {'element': 'Be10', 'low':  0.3, 'high': 3.0},
    }



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

"""
A function that was used to fit the SRIM data
"""    
def fit_function(x, a, b, c, d, j):
    return (a + b/x + c * np.exp(x*(-d)) + j/x**2)

"""
A function to calculate the energy correction for the particle in the second half of the target, after the reaction has occured
"""
def get_corrected_energy_target(energy, theta, ptype, half_target): 
    limit = LIMITS_TARGET[ptype]
    
    if energy < limit['low']: 
        print("energy of the particle in the target is below the correction range:", energy)
        return 0.
    else:
        distance = abs(half_target/np.cos(np.deg2rad(theta)))
        e_true = energy  
        
        while distance > CORRECTION_STEP:
            if e_true >= limit['high']:
                correction = CORRECT_TARGET_HIGH[ptype]
            else:  #the particles with energies below the lower limit are removed, so this is for the [low,high] interval
                correction = CORRECT_TARGET_LOW[ptype]
                
            e_true = e_true + fit_function(e_true,*correction['array'])/1000  #/1000 because the fit_function gives the output in keV and the detected energies are in MeV
            distance = distance - CORRECTION_STEP
                    
        return (e_true + fit_function(e_true,*correction['array'])/1000 * distance / CORRECTION_STEP)  

"""A function to calculate the energy correction for the particle in a dead layer of a detector"""
def get_corrected_energy_Dlayer(energy, ptype, dlayer):
    limit = LIMITS_DLAYER[ptype]
    
    if energy < limit['low']: 
        print("energy of the particle in a dead layer is below the correction range:", energy)
        return 0.
    else:
        distance = dlayer
        e_true = energy  
        
        while distance > CORRECTION_STEP:
            if e_true >= limit['high']:
                correction = CORRECT_DLAYER_HIGH[ptype]
            else:   #the particles with energies below the lower limit are removed, so this is for the [low,high] interval
                correction = CORRECT_DLAYER_LOW[ptype]
                
            e_true = e_true + fit_function(e_true,*correction['array'])/1000  #/1000 because the fit_function gives the output in keV and the detected energies are in MeV
            distance = distance - CORRECTION_STEP
                    
        return (e_true + fit_function(e_true,*correction['array'])/1000 * distance / CORRECTION_STEP)  
    
    
    
""" EXCECUTION """ 
       
"""open the file and get the tree"""
myfile = TFile(f'./ptype/{S_RUN}.root') #open
if DEBUG:
    myfile.ls() #check the file contents

t = myfile.Get(TREE_NAME)  #get the tree
if DEBUG:
    t.Print() #check the tree contents
    
"""creating the output root file; defining name, output tree, output tree variables, output tree structure """
out_file = TFile(f'./cnrg/{OUT_RUN}.root', 'recreate') # output root file
out_t = TTree('tree', f'{OUT_RUN}.root') # output tree


# output tree variables
# variables that exist in the starting tree
mult = zeros(1, dtype='short') # multiplicity of an event - number of hits, value of type Short
adc = zeros(MIN_mult * MAX_PARTICLES_ENTRY, dtype='short') # channel number, field of values of type Short
ampl = zeros(MIN_mult * MAX_PARTICLES_ENTRY, dtype='short') # amplitude, field of values of type Short
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
ptype = zeros(MAX_PARTICLES_ENTRY, dtype='short')# marking of the particle type, value of type Short
#new variables 
cnrg =  zeros(MAX_PARTICLES_ENTRY, dtype='float') # corrected energy, field of values of type Float

# output tree structure; adding new branches
out_t.Branch('event', event, 'event/I') # event number branch
out_t.Branch('cnc', cnc, 'cnc/S') # coincidences branch    
out_t.Branch('mult',mult,'mult/S') # multiplicity branch  
if PARTICLE_ENTRY:  
    out_t.Branch('adc',adc,'adc[3]/S') # channel number branch, field of lenght 3
    out_t.Branch('ampl',ampl,'ampl[3]/S') # amplitude branch, field of lenght 3
    out_t.Branch('nrg', nrg, 'nrg[3]/D') # energy branch, field of lenght 3
    out_t.Branch('cnrg', cnrg, 'cnrg/D') # corrected energy branch
    out_t.Branch('detector', detector, 'detector/S') # detector branch
    out_t.Branch('r', r, 'r/D') # particle distance branch
    out_t.Branch('theta', theta, 'theta/D') # particle angle branch
    out_t.Branch('phi', phi, 'phi/D') # particle angle branch
    out_t.Branch('ptype', ptype, 'ptype/S') # particle type branch
else:
    out_t.Branch('adc',adc,'adc[mult]/S') # channel number branch, field of lenght "mult"
    out_t.Branch('ampl',ampl,'ampl[mult]/S') # amplitude branch, field of lenght "mult"
    out_t.Branch('nrg', nrg, 'nrg[mult]/D') # energy branch, field of lenght "mult"
    out_t.Branch('cnrg', cnrg, 'cnrg[cnc]/D') # corrected energy branch, field of lenght "cnc"  
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

counter_below_range = 0

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
        #copying adc, ampl and energy info from the previous run          
        for i in range(3):
            adc[i] = t.adc[i]
            ampl[i] = t.ampl[i]
            nrg[i] = t.nrg[i]

        front_counter = 0 #counter that goes though field positions that contain information about the front strip
        dE_counter = 2 #counter that goes through field positions that contain information about the dE strip
        
        detector[0] = t.detector
        r[0] = t.r
        theta[0] = t.theta
        phi[0] = t.phi
        ptype[0] = t.ptype
        
        if ptype[0] == 0:
            cnrg[0] = 0
        else:
            cnrg_ = get_corrected_energy_Dlayer(t.nrg[front_counter], ptype[0], DLAYER)
            cnrg_ = get_corrected_energy_Dlayer(cnrg_, ptype[0], DLAYER)
            if cnrg_ <= 0.0:
                counter_below_range += 1
                continue

            cnrg_ = cnrg_ + t.nrg[dE_counter]
            cnrg_ = get_corrected_energy_Dlayer(cnrg_, ptype[0], DLAYER)

            cnrg[0] = get_corrected_energy_target(cnrg_, t.theta[0], ptype[0], HALF_TARGET)
            
            if cnrg[0] <= 0.0:
                counter_below_range += 1
                continue
            
        out_t.Fill() # filling the out tree; particle is the output entry
    else:
        #copying adc, ampl and energy info from the previous run            
        for i in range(mult[0]):
            adc[i]=t.adc[i]
            ampl[i]=t.ampl[i]
            nrg[i] = t.nrg[i] 

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
            ptype[i] = t.ptype[i]
            
            if ptype[i] == 0:
                cnrg[i] = 0
            else:
                cnrg_ = get_corrected_energy_Dlayer(t.nrg[front_counter], ptype[i], DLAYER)
                cnrg_ = get_corrected_energy_Dlayer(cnrg_, ptype[i], DLAYER)
                if cnrg_ <= 0.0:
                    counter_below_range += 1
                    continue

                cnrg_ = cnrg_ + t.nrg[dE_counter]
                cnrg_ = get_corrected_energy_Dlayer(cnrg_, ptype[i], DLAYER)

                cnrg[i] = get_corrected_energy_target(cnrg_, t.theta[i], ptype[i], HALF_TARGET[TARGET])
                
                if cnrg[i] == 0.0:
                    counter_below_range += 1
                    continue
            
        out_t.Fill() # filling the out tree; event is the output entry as in the raw root files


print("The number of particles with energies below the range of energy corrections:", counter_below_range)
        
out_file.Write()
out_file.Close()
myfile.Close()
