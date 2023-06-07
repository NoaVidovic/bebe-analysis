# -*- coding: utf-8 -*-
"""

@author: Deni

This is a script for turning detected amplitudes into energies.
Prior to use we need to finish the calibration i.e. get the calibration coefficients-
It is used during the dE calibration to turn amplitude into energies for E detectors in the used runs.
It is used in physical runs either before or after defining particles. In the final run versions this is done before defining particles so that we can use energy related particle filters.
"""

""" import predefined functions """
from ROOT import TFile, TTree
from numpy import genfromtxt, zeros
from time import time
from datetime import timedelta
from sys import argv


""" 
Define program-wide info
    !ADJUST!: N_RUN, TREE_NAME, N_ADC, COEFF_FILENAME, CALIB 
    PARTICLES(True,False), PARTICLE_ENTRY(True,False), EVENT_ENTRY(True,False), dE(True,False), CONFLICT_MARKING (True,False)
    DEBUG(True,False)
"""
S_RUN = 'run26' if len(argv) < 2 else argv[1]  # S_RUN is the first argument of the scripts, 'run26' if not given
S_RUN = f'{S_RUN}_killnoiseALL'
TREE_NAME = 'tree' #usually we used names "tree" or "T"

COEFF_FILENAME = 'calibALL_dEBe_RUTH.txt' #txt file with the calibration coefficients
CALIB = 'RUTH' # the name of the calibration

DEBUG = False #are we checking the output

# do we have defined particles in the starting run, is a particle or an event an entry in the starting tree
PARTICLES = False
PARTICLE_ENTRY = False
EVENT_ENTRY = False
CONFLICT_MARKING = False # do we have conflicted particles who share the same adc strip marked
INTERSTRIP_MARKING = False  # do we have potentital interstrip particles marked

dE = True  #is the dE detector calibrated

# number of ADC channels (strips), 128 if only E det is calibrated, 192 if dE is aswell
if dE:
    N_ADC = 192
else:
    N_ADC = 128

if dE:
    OUT_RUN = f'{S_RUN[:5]}_NRG({CALIB})' # the name of the output file 
else:
    OUT_RUN = f'{S_RUN[:5]}_NRG({CALIB})E'  # the name of the output file 

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
A function to calculate energies out of amplitudes and calibration coefficients
"""

def amplitude_2_energy(ampl, slope, intercept):
    return float(ampl) * slope + intercept
    
    
""" EXCECUTION """ 

"""open the file and get the tree"""
myfile = TFile(f'./killnoiseALL/{S_RUN}.root')  # open
# myfile.ls() #check the file contents
t = myfile.Get(TREE_NAME)   # get the tree
# t.Print() #check the tree contents    

"""creating the output root file; defining name, output tree, output tree variables, output tree structure """
out_file = TFile(f'./nrg/{OUT_RUN}.root', "recreate")  # output root file
out_t = TTree('tree', f'{OUT_RUN}.root')  # output tree

# output tree variables
# variables that exist in the starting tree
mult = zeros(1, dtype='short')  # multiplicity of an event - number of hits, value of type Short
adc = zeros(N_ADC, dtype='short')  # channel number, field of values of type Short
ampl = zeros(N_ADC, dtype='short')  # amplitude, field of values of type Short

if PARTICLES:
    event = zeros(1, dtype='int')  # the number of events in the starting run, value of type Int
    cnc = zeros(1, dtype='short')  # the number of particles in an event (coincidences), value of type Short
    detector = zeros(N_ADC, dtype='short')  # which detector the particle went into; 1,2,3 or 4, value of type Short
    if CONFLICT_MARKING:
        conflicts = zeros(1, dtype='short')  # the number of conflicting particles in an event, value of type Short
    if INTERSTRIP_MARKING:
        interstrip = zeros(1, dtype='short')  # the number of potential interstrip particles in an event, value of type Short
#new variable
nrg = zeros(N_ADC, dtype='float')  # energy, field of values of type Float


# output tree structure; adding energy branch to the structure of the input tree
if not PARTICLES:
    out_t.Branch('mult', mult, 'mult/S')  # multiplicity branch
    out_t.Branch('adc', adc, 'adc[mult]/S')  # channel number branch, field of lenght "mult"
    out_t.Branch('ampl', ampl, 'ampl[mult]/S')  # amplitude branch, field of lenght "mult"
    out_t.Branch('nrg', nrg, 'nrg[mult]/D')  # energy branch, field of lenght "mult"
else:
    out_t.Branch('event', event, 'event/I')  # event number branch
    out_t.Branch('cnc', cnc, 'cnc/S')  # coincidences branch
    out_t.Branch('mult', mult, 'mult/S')  # multiplicity branch

    if PARTICLE_ENTRY:  
        out_t.Branch('adc', adc, 'adc[3]/S')  # channel number branch, field of lenght 3
        out_t.Branch('ampl', ampl, 'ampl[3]/S')  # amplitude branch, field of lenght 3
        out_t.Branch('nrg', nrg, 'nrg[3]/D')  # energy branch, field of lenght 3
        out_t.Branch('detector', detector, 'detector/S')  # detector branch

    if EVENT_ENTRY: 
        out_t.Branch('adc', adc, 'adc[mult]/S')  # channel number branch, field of lenght "mult"
        out_t.Branch('ampl', ampl, 'ampl[mult]/S')  # amplitude branch, field of lenght "mult"
        out_t.Branch('nrg', nrg, 'nrg[mult]/D')  # energy branch, field of lenght "mult"
        out_t.Branch('detector', detector, 'detector[cnc]/S')  # detector branch, field of lenght "cnc"

    if CONFLICT_MARKING:
        out_t.Branch('conflicts', conflicts, 'conflicts/S')  #conflict branch

    if INTERSTRIP_MARKING:
        out_t.Branch('interstrip', interstrip, 'interstrip/S')  #interstrip branch

""" open the calibration file and get the calibration coefficients
turning a txt file into a field fit_coeff[0, ...., n-1], n= number of lines; one line of txt file is one subfield  """

fit_coeff = genfromtxt(COEFF_FILENAME, dtype=[('adc', '<i4'), ('slope', '<f4'), ('intercept', '<f4')])

if DEBUG:
    print(fit_coeff)


""" go trough the run data """
n_entries = t.GetEntriesFast()  # get the  number of entries in the tree
if DEBUG:
    print('number of entries: ', n_entries)
    # n_entries = 30

start_time = time()  #start the timer

for i in range(n_entries):  # for: 0 to n_entries-1 
    if i % 10000 == 0:  # print the timer info every 10000 entries
        print_time(start_time, i, n_entries)
        
    t.GetEntry(i)  # get an entry (event or a particle)
    
    """ calculating the output variables"""
    mult[0] = t.wm  # copying the number of hits in an event, change of variable name from "wm" to "mult"
    if PARTICLES:
        event[0] = t.event  # copying the cardinal number of the event
        cnc[0] = t.cnc  # copying the number of coincidences in an event

        if CONFLICT_MARKING:
            conflicts[0] = t.conflicts  # copying the number of conflicts in an event

        if INTERSTRIP_MARKING:
            interstrip[0] = t.interstrip  # copying the number of conflicts in an event
    else:
        for i in range(mult[0]):
            adc[i] = t.wadc[i]  # copying the channel number of hits in an event, change of variable name from "wadc" to "adc"
            ampl[i] = t.wampl[i]  # copying the amplitude of hits in an event, change of variable name from "wampl" to "ampl"

            if t.wadc[i] <= N_ADC:  # going trough the channels for which we have the calibration
                # calculating the energy
                coeff = fit_coeff[t.wadc[i]-1]  #the appropriate coefficients for the channel
                nrg[i] = amplitude_2_energy(t.wampl[i], coeff[1], coeff[2])

                if DEBUG:
                    print(f'Energy calculation: {nrg[i]}, Ampl: {ampl[i]}, Coeff: {fit_coeff[t.wadc[i]-1]}, Adc: {adc[i]}')
            else:
                nrg[i] = 0  # setting the energies for the uncalibrated channels to zero, this is a way of marking these channels while staying consistent with tree entries for following programs
            
    if PARTICLE_ENTRY:
        detector[0] = t.detector  #copying the detector number for the particles in the event
        for i in range(3):  # one entry is one particle, so there are 3 detector hits in an entry
            adc[i] = t.wadc[i]
            ampl[i] = t.wampl[i]
            if t.wadc[i] <= N_ADC:
                coeff = fit_coeff[t.wadc[i]-1]
                nrg[i] = amplitude_2_energy(t.wampl[i], coeff[1], coeff[2])
            else:
                nrg[i] = 0

    if EVENT_ENTRY:
        for i in range(cnc[0]):  #there is "cnc" particles in an event(entry)
            detector[i] = t.detector[i]
        for j in range(mult[0]):  #there is "mult" hits in an entry
            adc[j] = t.wadc[j]
            ampl[j] = t.wampl[j]
            if t.wadc[j] <= N_ADC:
                coeff = fit_coeff[t.wadc[j]-1]
                nrg[j] = amplitude_2_energy(t.wampl[j], coeff[1], coeff[2])
            else:
                nrg[j] = 0

    # write an entry into the output tree
    if DEBUG:  #output check
        num = mult[0]
        if PARTICLE_ENTRY:
            print('Entry: ', num, adc[:3], ampl[:3], nrg[:3])
        else:
            print('Entry: ', num, adc[:num], ampl[:num], nrg[:num])
    
    out_t.Fill()
         

""" writing the output file and closing input and output files """
out_file.Write()
out_file.Close()
myfile.Close()
