# -*- coding: utf-8 -*-
"""
Created on Thu Jun 16 13:43:00 2022

@author: Deni

This is a script that calculates excitation spectra of various nuclei from single particle detections. 
No need to separate "particle=entry" and "event=entry" since we are only using events which have a single detected particle (cnc=1).
No need to exclude possible interstrip and conflicting events since we are only using events which have a single detected particle (cnc=1).
"""

""" import predefined functions """
import numpy as np

from datetime import timedelta
from sys import argv
from time import time
from ROOT import TFile, TH2S, TH1F


""" 
Define program-wide info
    !ADJUST!: S_RUN_LIST, S_RUN
        (True,False): DEBUG, STRIP_HITS 
        choose projectile-target composition and the reaction channel
        choose particle filters and their options
            (True,False):DETECTORS,STRIP_RANGE,STRIP_PAIRS,STRIP_PAIRS_2
            DETECTORS_USED,detectors_used_s, front_low, front_high, dE_low, dE_high
"""
DEBUG = False

S_RUN_LIST = [ f'run{n_run}_TMIN2Be_particles(E=E,F1-4)_ptype(4He,6He,6Li,7Li,8Li,9Be,10Be)_CUT_cnrg' for n_run in range(18, 30) ]

ptype_DET = 409 if len(argv) < 2 else int(argv[1])
STRIP_USED = 47 if len(argv) < 3 else int(argv[2])

TREE_NAME = 'tree' #usually we used names "tree" or "T"

STRIP_HITS = False #do we count the number of hits in different strips
if STRIP_HITS:
    strip_hits = {}

""" REACTION DETAILS """
#choosing the projectile-target composition
Be9_Be9 = True

#choosing the reaction channel
if Be9_Be9:
    Ex8Be = False
    Ex9Be = True
    Ex10B = False
    Ex11B = False
    Ex12B = False
    Ex12C = False
    Ex14C = False

""" 
9Be+9Be projectile-target composition
REACTION CHANNELS:
    9Be(9Be,10Be)8Be, 8Be Excitation
    9Be(9Be,9Be)9Be, 9Be Excitation, elastic scattering
    9Be(9Be,8Li)10B, 10B Excitation
    9Be(9Be,7Li)11B, 11B Excitation
    9Be(9Be,6Li)12B, 12B Excitation
    9Be(9Be,6He)12C, 12C Excitation
    9Be(9Be,4He)14C, 14C Excitation
"""
if Be9_Be9:
    Energy_PROJECTILE_HALFtarget = 53.778 #projectile energy at half target where the reaction approximately occurs
    Mass_PROJECTILE = 9.0    #projectile mass
    Name_PROJECTILE = '9Be' #name of the projectile
    Name_TARGET = '9Be' #name of the target
    Energy_BEAM = '54' #projectile energy before entering the target

    #9Be(9Be,10Be)8Be reaction
    if ptype_DET == 410:
        Mass_DET = 10.0 #mass of the detected particle
        Name_DET = '10Be' #name of the detected particle
        
        Q0 = 5.147 #Q value of the reaction when all the particles are in their ground states; Q = E_det + E_undet - E_projectile
        
        Mass_UNDET = 8.0 #mass of the undetected particle
        Name_UNDET = '8Be'  #name of the undetected particle 

    #9Be(9Be,9Be)9Be reaction
    elif ptype_DET == 409:
        Mass_DET = 9.0 #mass of the detected particle
        Name_DET = '9Be' #name of the detected particle
        
        Q0 = 0.0 #Q value of the reaction when all the particles are in their ground states; Q = E_det + E_undet - E_projectile
        
        Mass_UNDET = 9.0 #mass of the undetected particle
        Name_UNDET = '9Be'  #name of the undetected particle
        
    #9Be(9Be,8Li)10B reaction       
    elif ptype_DET == 308:
        Mass_DET = 8.0 #mass of the detected particle
        Name_DET = '8Li' #name of the detected particle
        
        Q0 = -10.299 #Q value of the reaction when all the particles are in their ground states; Q = E_det + E_undet - E_projectile
        
        Mass_UNDET = 10.0 #mass of the undetected particle
        Name_UNDET = '10B'  #name of the undetected particle

    #9Be(9Be,7Li)11B reaction         
    elif ptype_DET == 307:
        Mass_DET = 7.0 #mass of the detected particle
        Name_DET = '7Li' #name of the detected particle
        
        Q0 = -0.877 #Q value of the reaction when all the particles are in their ground states; Q = E_det + E_undet - E_projectile
        
        Mass_UNDET = 11.0 #mass of the undetected particle
        Name_UNDET = '11B'  #name of the undetected particle        

    #9Be(9Be,6Li)12B reaction
    elif ptype_DET == 306:
        Mass_DET = 6.0 #mass of the detected particle
        Name_DET = '6Li' #name of the detected particle
        
        Q0 = -4.759 #Q value of the reaction when all the particles are in their ground states; Q = E_det + E_undet - E_projectile
        
        Mass_UNDET = 12.0 #mass of the undetected particle
        Name_UNDET = '12B'  #name of the undetected particle

    #9Be(9Be,6He)12C reaction
    elif ptype_DET == 206:
        Mass_DET = 6.0 #mass of the detected particle
        Name_DET = '6He' #name of the detected particle
        
        Q0 = 5.104   #Q value of the reaction when all the particles are in their ground states; Q = E_det + E_undet - E_projectile
        
        Mass_UNDET = 12.0 #mass of the undetected particle
        Name_UNDET = '12C'  #name of the undetected particle
        
    #9Be(9Be,4He)14C reaction
    elif ptype_DET == 204:
        Mass_DET = 4.0 #mass of the detected particle
        Name_DET = '4He' #name of the detected particle
        
        Q0 = 17.252  #Q value of the reaction when all the particles are in their ground states; Q = E_det + E_undet - E_projectile
        
        Mass_UNDET = 14.0 #mass of the undetected particle
        Name_UNDET = '14C'  #name of the undetected particle

    else:
        print('Invalid reaction channel')
        exit(0)

        
""" Choosing which particles to use
We can do this by filtering through one or more options depending on the wanted output
    - detectors; choose which to use
    - strip range; choose what range of strip to use
    - strip pairs; choose which strip pairs to use
 """
DETECTORS = True #do we filter by detectors
STRIP_RANGE = False #do we filter by strip range
MATCHES = 's'

filters_used = [] #field of string markings of all filters used, to add to output file name
output_suffix = '' #string marking of all filters used, to add to output file name
#choose which detectors to use
if DETECTORS: 
    DETECTORS_USED = [2, 3] #1,2,3,4 write all that we want to include
    detectors_used_s = 'detBC_' #string marking for the output
    
    filters_used.append(detectors_used_s)

#bottom and upper front and dE strip number that we let into further calculation                
if STRIP_RANGE:
    front_low, front_high = STRIP_USED, STRIP_USED
    dE_low, dE_high = 112, 192
    
    strip_range_s = f'E=[{front_low},{front_high}],dE=[{dE_low},{dE_high}]_'
    filters_used.append(strip_range_s)


#best matched strip pairs    
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

    if MATCHES == 'a':
        strip_pairs_s = 'allPairs_'
    elif MATCHES == 's':
        strip_pairs_s = 'bestPairs_'

    filters_used.append(strip_pairs_s)

strip_corr_pars = {}
with open('strip_9be_correction', 'r') as f:
    for l in f.readlines():
        l = l.split()

        strip = int(l[0])
        pars = [ float(x) for x in l[1:] ]

        strip_corr_pars[strip] = pars

# the name of the output file 
output_suffix = ''.join(filters_used)

out_name = f"single_corr_sigma_Ex{Name_UNDET}_from{Name_DET}_{output_suffix}run18-30"

        
""" Defining histograms"""
#defining limits and bin size
res_ampl, low_ampl, high_ampl  = 3000,0,3000 #1 amplitude resolution
res_energy_detected, low_energy_detected, high_energy_detected = 1100, 0, 55 #50keV resolution

res_energy_E, low_energy_E, high_energy_E  = 5000,0,50 #10keV resolution
res_energy_dE, low_energy_dE, high_energy_dE  = 3000,0,30 #10keV resolution
res_energy_Ex, low_energy_Ex, high_energy_Ex  = 1100,-5,50 #50keV resolution

res_theta, low_theta, high_theta  = 110, 10, 65 #0.5 degrees resolution

#defining histograms
h0 = TH1F('AMPLdet', f"Amplitude {Name_DET}", res_ampl, low_ampl, high_ampl) # amplitude of the detected particle in the front E detector
h = TH1F('Edet_nrg', f"Energies in E detector of {Name_DET}", res_energy_detected, low_energy_detected, high_energy_detected ) # energy of the detected particle in the front E detector
h1 = TH1F('Edet_cnrg', f"Total energies of {Name_DET}", res_energy_detected, low_energy_detected, high_energy_detected ) # total (corrected) energy of the detected particle
h1_c = TH1F('E\'det_cnrg', f"Total energies of {Name_DET} (corrected)", res_energy_detected, low_energy_detected, high_energy_detected ) # total (corrected) energy of the detected particle
h2 = TH2S('dE-E_ampl', f"deltaE-E amplitudes of {Name_DET}", res_ampl, low_ampl, high_ampl, res_ampl, low_ampl, high_ampl) # dE-E amplitued graph of detected particles
h2.SetOption("COLZ") # contrast scale
h3 = TH2S('dE-E', f"deltaE-E energies of {Name_DET}", res_energy_E, low_energy_E, high_energy_E, res_energy_dE, low_energy_dE, high_energy_dE) # dE-E energy graph of detected particles
h3.SetOption("COLZ") # contrast scale

h4 = TH1F('Ex', f'Ex({Name_UNDET}) from {Name_DET} single detections, Beam {Name_PROJECTILE} {Energy_BEAM} MeV, Target {Name_TARGET}', res_energy_Ex, low_energy_Ex, high_energy_Ex) #Excitation energy of the undetected particle
h5 = TH2S('Ex-theta', f"theta - Ex({Name_UNDET}) from {Name_DET} single detections, Beam {Name_PROJECTILE} {Energy_BEAM} MeV, Target Name_TARGET", res_energy_Ex, low_energy_Ex, high_energy_Ex, res_theta, low_theta, high_theta) #Excitation energy of the undetected particle vs angle theta of the detected particle
h5_c = TH2S('Ex-theta_colour', f"theta - Ex({Name_UNDET}) from {Name_DET} single detections, Beam {Name_PROJECTILE} {Energy_BEAM} MeV, Target Name_TARGET", res_energy_Ex, low_energy_Ex, high_energy_Ex, res_theta, low_theta, high_theta) # COLOUR version of h5
h5_c.SetOption("COLZ") # contrast scale

h6 = TH1F('Ex\'', f'Ex\'({Name_UNDET}) from {Name_DET} single detections, Beam {Name_PROJECTILE} {Energy_BEAM} MeV, Target {Name_TARGET}', res_energy_Ex, low_energy_Ex, high_energy_Ex) #Excitation energy of the undetected particle
h7 = TH2S('Ex\'-theta', f"theta - Ex\'({Name_UNDET}) from {Name_DET} single detections, Beam {Name_PROJECTILE} {Energy_BEAM} MeV, Target Name_TARGET", res_energy_Ex, low_energy_Ex, high_energy_Ex, res_theta, low_theta, high_theta) #Excitation energy of the undetected particle vs angle theta of the detected particle
h7_c = TH2S('Ex\'-theta_colour', f"theta - Ex\'({Name_UNDET}) from {Name_DET} single detections, Beam {Name_PROJECTILE} {Energy_BEAM} MeV, Target Name_TARGET", res_energy_Ex, low_energy_Ex, high_energy_Ex, res_theta, low_theta, high_theta) # COLOUR version of h5
h7_c.SetOption("COLZ") # contrast scale

h8 = TH1F('Ex\'\'', f'Ex\'\'({Name_UNDET}) from {Name_DET} single detections, Beam {Name_PROJECTILE} {Energy_BEAM} MeV, Target {Name_TARGET}', res_energy_Ex, low_energy_Ex, high_energy_Ex) #Excitation energy of the undetected particle
h9 = TH2S('Ex\'\'-theta', f"theta - Ex\'\'({Name_UNDET}) from {Name_DET} single detections, Beam {Name_PROJECTILE} {Energy_BEAM} MeV, Target Name_TARGET", res_energy_Ex, low_energy_Ex, high_energy_Ex, res_theta, low_theta, high_theta) #Excitation energy of the undetected particle vs angle theta of the detected particle
h9_c = TH2S('Ex\'\'-theta_colour', f"theta - Ex\'\'({Name_UNDET}) from {Name_DET} single detections, Beam {Name_PROJECTILE} {Energy_BEAM} MeV, Target Name_TARGET", res_energy_Ex, low_energy_Ex, high_energy_Ex, res_theta, low_theta, high_theta) # COLOUR version of h5
h9_c.SetOption("COLZ") # contrast scale

state0_nrg = 0        # MeV
state1_nrg = 2.43     # MeV
nrg_half_width = 0.5  # MeV

hE0t = TH1F('N-theta', f'theta - N for particles in the range of {state0_nrg-nrg_half_width}-{state0_nrg+nrg_half_width} MeV', res_theta, low_theta, high_theta)
hE1t = TH1F('N-theta', f'theta - N for particles in the range of {state1_nrg-nrg_half_width}-{state1_nrg+nrg_half_width} MeV', res_theta, low_theta, high_theta)


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
A function for excitation energy calculation; excitation energy of the undetected particle from the total energy of the detected particle
"""    

def calculate_Ex_single (Q0, m_beam, m_det, m_undet, beam_energy, E_det, theta_det):
    temp1 = ((m_undet - m_beam)/m_undet) * beam_energy
    temp2 = ((m_undet + m_det)/m_undet) * E_det
    temp3 = m_beam * m_det * beam_energy * E_det
    temp4 = (2.0/m_undet) * np.sqrt(temp3) * np.cos(theta_det)
              
    Ex = Q0 + temp1 - temp2 + temp4
        
    return Ex   

    
""" EXECUTION """ 
print(f"Excitation of {Name_UNDET} from the detection of {Name_DET}, Filters:{output_suffix}  Run: 18-30")
#counter of total programme excecution time
start_time = time() #start the excecution timer

#counters
counter_particles_single = 0 #counter for all single detections in the input runs
counter_particles_reaction = 0 #counter for all the particles that are of the wanted type
counter_particles_filtered = 0 #counter for all the particles that are of the wanted type that have gone through the chosen filters

for s_run in S_RUN_LIST:
    """open the file and get the tree"""
    myfile = TFile(f'./cnrg/{s_run}.root') #open
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
        if i % 50000 == 0: # print the timer info every 50000 entries
            print_time(start_time_run, i, n_entries)
        
        t.GetEntry(i) # get an entry

        #field positions of front, back and dE information
        front_number = 0
        dE_number = 2
        
        #don't include events that have more than one detected particle
        if t.cnc != 1:
            continue
        counter_particles_single += 1
        
        
        #don't include particles that are not a part of the chosen reaction
        if t.ptype[0] != ptype_DET:
            continue
        counter_particles_reaction += 1
        
        #don't include particles that are not detected in the detectors that we want to use
        if DETECTORS and t.detector[0] not in DETECTORS_USED: 
            continue
        
        #don't include particles that are not in the strip range that we want to use
        if STRIP_RANGE and not (dE_low <= t.adc[dE_number] <= dE_high and front_low <= t.adc[front_number] <= front_high):
            continue
        
        #don't include particles that are not in the strip pairs
        if t.adc[dE_number] not in list_pairs[t.adc[front_number]]:
            continue
              
        counter_particles_filtered += 1

        #calculate the excitation energy of the undetected particle
        Ex = calculate_Ex_single(Q0, Mass_PROJECTILE, Mass_DET, Mass_UNDET, Energy_PROJECTILE_HALFtarget, t.cnrg[0], np.deg2rad(t.theta[0]))
        a, b = strip_corr_pars[t.adc[front_number]]

        Ex_ = a*Ex + b
        E_ = (2-a)*t.cnrg[0] - b

        Ex_2 = calculate_Ex_single(Q0, Mass_PROJECTILE, Mass_DET, Mass_UNDET, Energy_PROJECTILE_HALFtarget, E_, np.deg2rad(t.theta[0]))
         
        h0.Fill(t.ampl[0]) # amplitude of the detected particle in the front E detector
        h.Fill(t.nrg[0]) # energy of the detected particle in the front E detector
        h1.Fill(t.cnrg[0]) # total (corrected) energy of the detected particle
        h1_c.Fill(E_)  # total (corrected corrected) energy of the detected particle
        h2.Fill(t.ampl[0],t.ampl[2]) # dE-E amplitued graph of detected particles
        h3.Fill(t.nrg[0],t.nrg[2])  # dE-E energy graph of detected particles
        h4.Fill(Ex) #Excitation energy of the undetected particle
        h5.Fill(Ex,t.theta[0]) #Excitation energy of the undetected particle vs angle theta of the detected particle
        h5_c.Fill(Ex,t.theta[0]) # COLOUR version of h5
        h6.Fill(Ex_) #Excitation energy of the undetected particle
        h7.Fill(Ex_, t.theta[0]) #Excitation energy of the undetected particle vs angle theta of the detected particle
        h7_c.Fill(Ex_, t.theta[0]) # COLOUR version of h5
        h8.Fill(Ex_2) #Excitation energy of the undetected particle
        h9.Fill(Ex_2, t.theta[0]) #Excitation energy of the undetected particle vs angle theta of the detected particle
        h9_c.Fill(Ex_2, t.theta[0]) # COLOUR version of h5
        
        if state0_nrg - nrg_half_width <= Ex_ <= state0_nrg + nrg_half_width:
            hE0t.Fill(t.theta[0])

        if state1_nrg - nrg_half_width <= Ex_ <= state1_nrg + nrg_half_width:
            hE1t.Fill(t.theta[0])

        
        if STRIP_HITS:
            for i in range(t.mult):
                try:
                    strip_hits[t.adc[i]] += 1
                except:#setting up the counter (dictionary member) for the strip at its first mention
                    strip_hits[t.adc[i]] = 1   
        
#when there are multiple runs, we calculate the total elapsed time
total_time = round(time()-start_time)
total_time = timedelta(0,total_time)  

print("Writing histograms into a file")
g = TFile(f"./singles/{out_name}.root", "recreate")       
h0.Write()
h.Write()
h1.Write()
h1_c.Write()
h2.Write()
h3.Write()
h4.Write()
h5.Write()
h5_c.Write()
h6.Write()
h7.Write()
h7_c.Write()
h8.Write()
h9.Write()
h9_c.Write()
hE0t.Write()
hE1t.Write()
g.Close()
print(f"Excitation of {Name_UNDET} from the detection of {Name_DET}, Filters: {output_suffix}  Run: 18-30")
print("The number of single detections in the input runs:", counter_particles_single)
print("The number of particles that are of the wanted type:", counter_particles_reaction)
print("The number of particles that are of the wanted type that have made it through the filters:", counter_particles_filtered)
if STRIP_HITS:
    print("Strip hits:", strip_hits)
print("The total time of excecution", total_time)
print("END")
