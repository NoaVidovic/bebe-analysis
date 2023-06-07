# -*- coding: utf-8 -*-
"""
Created on Fri Jun 17 00:26:56 2022

@author: pmf1
This is a script that calculates excitation spectra of various nuclei from coincident particle detections. 
No need to separate "particle=entry" and "event=entry" since we are only using events which have two detected particles (cnc=2).
"""

""" import predefined functions """
from datetime import timedelta
from sys import argv
from time import time
from ROOT import TFile, TH1F, TH2S  # type: ignore

import numpy as np

from el import el_name

""" 
Define program-wide info
    !ADJUST!: S_RUN_LIST, S_RUN
        (True,False): DEBUG, STRIP_HITS, REMOVE_INTERSTRIP, REMOVE_CONFLICTS , SAME_PARTICLES, ASSYMETRY_CHECK, ptype1_DET2
        choose projectile-target composition and the reaction channel
        choose particle filters and their options
            (True,False):DETECTORS,STRIP_RANGE,STRIP_PAIRS,STRIP_PAIRS_2, ROMANO_CUT; QthetaC_CUT, Q_CUT, Ex_compound_CUT
            DET1, DET2, detectors_used_1, detectors_used_2
            front_low, front_high, dE_low, dE_high
            f_romano, cut_number, cut_number_used, cut_number_used_S
            f_QthetaS, cut_number, cut_number_used, cut_number_used_S
            Q_value_MIN, Q_value_MAX
            Ex_compound_MIN,Ex_compound_MAX
"""
DEBUG = False

CALIB = 'RUTH90'
S_RUN_LIST = [ f'run{n_run}_{CALIB}_particles(E=E,F1-4)_ptype(4He,6He,6Li,7Li,8Li,9Be,10Be)_CUT_cnrg' for n_run in range(18, 30) ]

ptype1 = 204 if len(argv) < 2 else int(argv[1])
ptype2 = 204 if len(argv) < 3 else int(argv[2])

det1_raw = '2' if len(argv) < 5 else argv[3]
det2_raw = '2' if len(argv) < 5 else argv[4]

TREE_NAME = 'tree' #usually we used names "tree" or "T"

STRIP_HITS = False #do we count the number of hits in different strips
if STRIP_HITS:
    strip_hits = {}

#are we removing potential interstrip events (REMOVE_INTERSTRIP=True) and/or conflicting events (REMOVE_CONFLICTS=True)
REMOVE_INTERSTRIP = True
REMOVE_CONFLICTS = True

#are the two detected particles of the same type so we have to be carefull about not counting the same events twice
SAME_PARTICLES = ptype1 == ptype2

#are we checking for the assimetry in detection; each particle comes from a different detector
# ASSYMETRY_CHECK = False
# if ASSYMETRY_CHECK: 
#     # which of the two possible combinations are we considering; for example 4He comes from A detector and 6He from B or vice a versa
#     ptype1_DET2 = True
    
""" REACTION DETAILS """
#choosing the projectile-target composition
Be9_Be9 = True


#choosing the reaction channel by the combination of detected particles
He4He4 = ptype1 == ptype2 == 204
He6He4 = ptype1 == 204 and ptype2 == 206
    

""" 
9Be+9Be projectile-target composition
SCENARIOUS:
    SCENARIO 1; detected particles coming from a compound particle decay
    SCENARIO 2 and 3; one detected particle coming from a compound particle decay, and the other comes from the reaction mechanism
REACTION CHANNELS:
    9Be(9Be,4He4He); 8Be,10Be,14C Excitation
        SCENARIO 1; 4He+4He from 8Be, the remaining nucleus is 10Be --> 8Be and 10Be Excitation spectra
        SCENARIO 2 and 3; 4He+14C exit channel, the other 4He comes from 14C decay --> 14C Excitation spectra
    9Be(9Be,4He6He); 8Be,10Be,12C,14C Excitation
        SCENARIO1; 4He+6He from 10Be, the remaining nucleus is 8Be --> 8Be and 10Be Excitation spectra
        SCENARIO2; 6He + 12C exit channel, 4He comes from 12C decay --> 12C Excitation spectra
        SCENARIO3; 4He+14C exit channel, 6He comes from 14C decay --> 14C Excitation spectra
""" 
if Be9_Be9:
    Energy_PROJECTILE_HALFtarget = 53.778 #projectile energy at half target where the reaction approximately occurs
    ptype_PROJECTILE, ptype_TARGET = 409, 409  # projectile and target particle type

    Mass_PROJECTILE, Mass_TARGET = ptype_PROJECTILE % 100, ptype_TARGET % 100  # projectile and target mass
    Name_PROJECTILE, Name_TARGET = el_name(ptype_PROJECTILE), el_name(ptype_TARGET)  # projectile and target name

    Energy_BEAM = '54' #projectile energy before entering the target

    M_DET1, M_DET2 = ptype1 % 100, ptype2 % 100
    Name_DET1, Name_DET2 = el_name(ptype1), el_name(ptype2)

    # for scenario 1
    M_COMPOUND = M_DET1 + M_DET2
    Name_UNDET_compound = el_name((ptype1 // 100 + ptype2 // 100) * 100 + int(M_COMPOUND))

    M_UNDET_THIRD = Mass_PROJECTILE + Mass_TARGET - M_COMPOUND
    
    if He4He4:
        #SCENARIO 1; 
        E_THRESHOLD = -0.09184  # threshold energy for 8Be decay into 4He+4He
        Q0 = 5.240  # Qvalue for the 3particle reaction when all the particles are in their ground states; Q = E_det1 + E_det2 + E_undet - E_projectile
        
        #SCENARIO 2 and 3; for same particles both scenarios are the same, but the values are duplicated to be in accordance with other reaction channels
        Name_UNDET1 = '14C' #name of the undetected particle
        M_UNDET1 = 14.0   #mass of the undetected particle
        Q0_single1 = 17.252   # Q value for 4He+14C reaction when all the particles are in their ground states; Q = E_det + E_undet - E_projectile
        
        Name_UNDET2 = '14C'  #name of the undetected particle
        M_UNDET2 = 14.0 #mass of the undetected particle
        Q0_single2 = 17.252   # Q value for 4He+14C reaction when all the particles are in their ground states; Q = E_det + E_undet - E_projectile 
    elif He6He4:
        E_THRESHOLD= 7.4095 # threshold energy for 10Be decay into 4He+6He
        Q0= -2.262
        
        Name_UNDET1 = '12C' 
        M_UNDET1 = 12.0
        Q0_single1 = 5.104
        
        Name_UNDET2 = '14C'  
        M_UNDET2 = 14.0
        Q0_single2 = 17.252      
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
ROMANO_CUT = False  #do we filter by cuts on Romano plot
QthetaC_CUT = False #do we filter by cuts on Q vs. angle of compount particle plot
Q_CUT = False #do we filter by ranges on Q spectrum
Ex_compound_CUT = False #do we filter by ranges on the excitation energy of the compount particle spectrum

filters_used = [] #field of string markings of all filters used, to add to output file name

#choose which detectors to use for each particle, detector hits always recorded in the ascending order i.e. DET1<DET2
if DETECTORS: 
    det_names = 'ABCD'

    DET1 = [ int(c) for c in det1_raw ] #for the first detected particle
    DET2 = [ int(c) for c in det2_raw ] #for the second detected particle
    detectors_used_1 = ''.join([ det_names[i-1] for i in sorted(DET1) ]) 
    detectors_used_2 = ''.join([ det_names[i-1] for i in sorted(DET2) ]) 
    
    filters_used.append(f'{Name_DET1}[{detectors_used_1}]{Name_DET2}[{detectors_used_2}]_')
else:
    filters_used.append(f'{Name_DET1}{Name_DET2}_')


#bottom and upper front and dE strip number that we let into further calculation                
if STRIP_RANGE:
    front_low, front_high = 33,48
    dE_low, dE_high = 145,160
    
    strip_range_s = 'E=[%i,%i],dE=[%i,%i]_'%(front_low,front_high, dE_low,dE_high)
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

#Romano plot CUT
if ROMANO_CUT:
    f_romano = TFile('cut_romano_4He6He_detB_3cuts.root')
    cut_romano = []
    cut_number = 3 #number of cuts made
    for i in range(cut_number):
        cut_romano.append(f_romano.Get('cut%i;1'%(i+1)))
    cut_number_used = [0] #list of cuts that we are using of all that we made 
    cut_number_used_S = "1" #marking of used cuts
    cut_romano_s= 'romanoCUT[cuts%s]_'%(cut_number_used_S)
    filters_used.append(cut_romano_s)

#Q value vs angle of compound particle CUT
if QthetaC_CUT:
    f_QthetaS = TFile('cut_QthetaS_4He[B]4He[C]_4cut.root')
    cut_QthetaS = []
    cut_number = 4
    for i in range(cut_number):
        cut_QthetaS.append(f_QthetaS.Get('cut%i;1'%(i+1)))
    cut_number_used = [1] #list of cuts that we are using of all that we made 
    cut_number_used_S = "2" #marking of used cuts
    cut_QthetaC_s= 'QthetaCUT[cuts%s]_'%(cut_number_used_S)
    filters_used.append(cut_QthetaC_s)
#Q-Q0 value limits, only one interval for now, if more intervals are needed make a field of min and max values and adapt the further code    
if Q_CUT: 
    Q_value_MIN = 1.0
    Q_value_MAX = 2.7
    cut_Q_s= 'QCUT[%f,%f]_'%(Q_value_MIN,Q_value_MAX)
    filters_used.append(cut_Q_s)
#Ex of compound nucleus value limits, only one interval for now, if more intervals are needed make a field of min and max values and adapt the further code  
if Ex_compound_CUT:
    Ex_compound_MIN = 1.0
    Ex_compound_MAX = 50.0
    cut_ExCompound_s= f'ExCompoundCUT[{Ex_compound_MIN},{Ex_compound_MAX}]_'
    filters_used.append(cut_ExCompound_s)

if REMOVE_INTERSTRIP:
    filters_used.append("NOinter_")
if REMOVE_CONFLICTS:
    filters_used.append("NOconf_")
    
filters_used.append("cnc=2_")

output_suffix = ''.join(filters_used)
out_name = f"cnc/coinc_{CALIB}_{output_suffix}run18-30.root"


    
""" Defining histograms"""
#defining limits and bin size
res_energy_detected, low_energy_detected, high_energy_detected = 2200, -1, 54 #25keV resolution
res_energy_dE, low_energy_dE, high_energy_dE  = 1200,0,30 #25keV resolution
res_energy_Ex, low_energy_Ex, high_energy_Ex  = 2400,-10,50 #25keV resolution

res_romano_E, low_romano_E, high_romano_E  = 8000, -20, 60 #10keV resolution
res_romano_P, low_romano_P, high_romano_P  = 12000, -100, 500 #50keV resolution

res_Q, low_Q, high_Q  = 500, -25, 25  # 100 keV resolution
if ptype1 == ptype2:
    res_Q = 2000  # 25 keV

res_theta, low_theta, high_theta  = 225, 10, 65 #0.2 degrees resolution


h1 = TH2S('dEdet_vs_Edet', "deltaE vs E of chosen detected particles", res_energy_detected, low_energy_detected, high_energy_detected , res_energy_dE, low_energy_dE, high_energy_dE ) # dE-E of chosen detected particles
h2 = TH2S('EDET1_vs_EDET2', "total energy vs total energy of chosen detected particles ", res_energy_detected, low_energy_detected, high_energy_detected , res_energy_detected, low_energy_detected, high_energy_detected ) # E-E of chosen detected particles
h3 = TH2S('Romano', "E_til vs P_til of chosen detected particles", res_romano_P, low_romano_P, high_romano_P, res_romano_E, low_romano_E, high_romano_E) # romano plot of chosen detected particles
h3colour = TH2S('Romano_colour', "E_til vs P_til of chosen detected particles", res_romano_P, low_romano_P, high_romano_P, res_romano_E, low_romano_E, high_romano_E) # romano plot of chosen detected particles
h3colour.SetOption("COLZ")  
h4 = TH1F('Q0-Q', 'Q0 - Q spectra i.e. Excitation of third particle ',res_Q, low_Q, high_Q) # Q0-Q  spectrum of chosen detected particles
h5 = TH2S('thetaCompound_vs_Q0-Q', " angle %s vs Q0-Q"%(Name_UNDET_compound), res_Q, low_Q, high_Q, res_theta, low_theta, high_theta) # Q0-Q vs theta of compound nucleus spectrum of chosen detected particles
h5colour = TH2S('thetaCompound_vs_Q0-Q_colour', " angle %s vs Q0-Q"%(Name_UNDET_compound), res_Q, low_Q, high_Q, res_theta, low_theta, high_theta) # Q0-Q  vs theta of compound nucleus spectrum of chosen detected particles
h5colour.SetOption("COLZ")
h6 = TH1F('ExCompound','Ex of compound nucleus %s'%(Name_UNDET_compound),res_energy_Ex, low_energy_Ex, high_energy_Ex) # Excitation of compound nucleus
h7 = TH2S('ExCompound_vs_ExNotDET1', "Ex(%s) vs Ex(%s)"%(Name_UNDET_compound, Name_UNDET1), res_energy_Ex, low_energy_Ex, high_energy_Ex, res_energy_Ex, low_energy_Ex, high_energy_Ex) # Excitation of compound nucleus vs Excitation of 1st possible single scenario
h7colour = TH2S('ExCompound_vs_ExNotDET1_colour', "Ex(%s) vs Ex(%s)"%(Name_UNDET_compound, Name_UNDET1), res_energy_Ex, low_energy_Ex, high_energy_Ex, res_energy_Ex, low_energy_Ex, high_energy_Ex)# Excitation of compound nucleus vs Excitation of 1st possible single scenario
h7colour.SetOption("COLZ")
h8 = TH2S('ExCompound_vs_ExNotDET2', "Ex(%s) vs Ex(%s)"%(Name_UNDET_compound, Name_UNDET2), res_energy_Ex, low_energy_Ex, high_energy_Ex, res_energy_Ex, low_energy_Ex, high_energy_Ex) # Excitation of compound nucleus vs Excitation of 2nd possible single scenario
h8colour = TH2S('ExCompound_vs_ExNotDET2_colour', "Ex(%s) vs Ex(%s)"%(Name_UNDET_compound, Name_UNDET2), res_energy_Ex, low_energy_Ex, high_energy_Ex, res_energy_Ex, low_energy_Ex, high_energy_Ex) # Excitation of compound nucleus vs Excitation of 2nd possible single scenario
h8colour.SetOption("COLZ")
if SAME_PARTICLES: #combining particles from two ingle scnearios becase the scenarios are the same for same particle detections
    h_same = TH2S('ExCompound_vs_ExNotDET2_same', "Ex(%s) vs Ex(%s) addition of previous two graphs for same particles"%(Name_UNDET_compound, Name_UNDET2), res_energy_Ex, low_energy_Ex, high_energy_Ex, res_energy_Ex, low_energy_Ex, high_energy_Ex) # Excitation of compound nucleus vs Excitation of 2nd possible single scenario
    h_same_colour = TH2S('ExCompound_vs_ExNotDET2_colour', "Ex(%s) vs Ex(%s)"%(Name_UNDET_compound, Name_UNDET2), res_energy_Ex, low_energy_Ex, high_energy_Ex, res_energy_Ex, low_energy_Ex, high_energy_Ex) # Excitation of compound nucleus vs Excitation of 2nd possible single scenario
    h_same_colour.SetOption("COLZ")
h9 = TH2S('ExNotDET1_vs_ExNotDET2', "Ex(%s) vs Ex(%s)"%(Name_UNDET1, Name_UNDET2), res_energy_Ex, low_energy_Ex, high_energy_Ex, res_energy_Ex, low_energy_Ex, high_energy_Ex) # Excitation of 1st possible single scenario vs Excitation of 2nd possible single scenario
h9colour = TH2S('ExNotDET1_vs_ExNotDET2_colour', "Ex(%s) vs Ex(%s)"%(Name_UNDET1, Name_UNDET2), res_energy_Ex, low_energy_Ex, high_energy_Ex, res_energy_Ex, low_energy_Ex, high_energy_Ex) # Excitation of 1st possible single scenario vs Excitation of 2nd possible single scenario
h9colour.SetOption("COLZ")
h10 = TH1F('ExNotDET1','Ex %s'%(Name_UNDET1),res_energy_Ex, low_energy_Ex, high_energy_Ex) # Excitation of 1st possible single scenario 
h11 = TH1F('ExNotDET2','Ex %s'%(Name_UNDET2),res_energy_Ex, low_energy_Ex, high_energy_Ex) # Excitation of 2nd possible single scenario 
if SAME_PARTICLES: 
    h_same_projection = TH1F('ExNotDET','Ex %s addition of previous two graphs for same particles '%(Name_UNDET1),res_energy_Ex, low_energy_Ex, high_energy_Ex) # Excitation of 1st possible single scenario



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
Romano plot values E_til and P_til which are conected with the third undetected particle energy and impuls
"""
 
def romano(Energy_PROJECTILE_HALFtarget, Mass_PROJECTILE, m1, m2, E_DET1, E_DET2, theta1, theta2, phi1, phi2):
    e_til = Energy_PROJECTILE_HALFtarget  - E_DET1 - E_DET2
    
    pp_par = np.sqrt(2*Mass_PROJECTILE*Energy_PROJECTILE_HALFtarget )
    p1_par = np.sqrt(2*m1*E_DET1)*np.cos(theta1)
    p2_par = np.sqrt(2*m2*E_DET2)*np.cos(theta2)
    
    p3_par = pp_par - p1_par - p2_par
    
    p1_ver = np.sqrt(2*m1*E_DET1)*np.sin(theta1)
    p2_ver = np.sqrt(2*m2*E_DET2)*np.sin(theta2)
    
    p3_ver_squared = p1_ver**2 + p2_ver**2 + 2*p1_ver*p2_ver*np.cos(phi1-phi2)
    
    p3_squared = p3_par**2 + p3_ver_squared
    
    p_til = p3_squared / 2
    
    return (p_til, e_til)
    

"""
Excitation energy calculations
    - Excitation of the compound nucleus
    - Excitation of possible single detection scenarios
"""    
def calculate_Ex_compound(E_THRESHOLD, m1, m2, E_DET1, E_DET2, theta1, theta2, phi1, phi2):
    m_red=m1*m2/(m1+m2)
    costheta_rel=np.cos(theta1)*np.cos(theta2)+np.sin(theta1)*np.sin(theta2)*np.cos(phi1-phi2)
    temp1=E_DET1/m1
    temp2=E_DET2/m2
    temp3=2*np.sqrt(temp1*temp2)*costheta_rel
    
    E_rel=m_red*(temp1+temp2-temp3)
    
    Ex=E_rel+E_THRESHOLD
    
    return Ex

def calculate_Ex_single(Q0, Mass_PROJECTILE, m_det, m_undet, Energy_PROJECTILE_HALFtarget , E_det, theta_det):
    temp1 = ((m_undet - Mass_PROJECTILE)/m_undet) * Energy_PROJECTILE_HALFtarget 
    temp2 = ((m_undet + m_det)/m_undet) * E_det
    temp3 = Mass_PROJECTILE * m_det * Energy_PROJECTILE_HALFtarget  * E_det
    temp4 = (2.0/m_undet) * np.sqrt(temp3) * np.cos(theta_det)

    Ex = Q0 + temp1 - temp2 + temp4
    
    return Ex
    
    
    
"""
Angle of the compound nucleus
"""   

def calculate_theta_compound(E_x, E_THRESHOLD, m1, m2, E_DET1, E_DET2, theta1, theta2, phi1, phi2, m_compound):
    E_compound = E_DET1 + E_DET2 + E_THRESHOLD - E_x
    
    temp1 = np.sqrt(2*m1*E_DET1)*np.cos(theta1) + np.sqrt(2*m2*E_DET2)*np.cos(theta2)
    temp2 = np.sqrt(2*m1*E_DET1)*np.sin(theta1)*np.sin(phi1)
    temp3 = np.sqrt(2*m2*E_DET2)*np.sin(theta2)*np.sin(phi2)
    
    theta_compound = np.arccos(temp1/np.sqrt(2*m_compound*E_compound))
    phi_compound = np.arcsin((temp2+temp3)/(np.sqrt(2*m_compound*E_compound)*np.sin(theta_compound)))
    
    theta_compound_deg = theta_compound*180/np.pi
    phi_compound_deg = phi_compound*180/np.pi
        
    
    return (theta_compound_deg, phi_compound_deg)
    




""" EXCECUTION"""
print(f"Detections: {Name_DET1}+{Name_DET2}, Filters used: {output_suffix}, run18-30")
    
#counter of total programme excecution time
total_time = 0
start_time = time() #start the excecution timer 


#counters
counter_particles = 0 #total number of particles in the used runs
counter_particles_coincident = 0 #Number of particles in pairs of the chosen type
counter_particles_filter = 0 #Number of particles that have gone through the chosen filters


    
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
        counter_particles +=1
    
        if i % 50000 == 0: # print the timer info every 50000 entries
            print_time(start_time_run,i,n_entries)
        
        t.GetEntry(i) # get an entry
        
        
        #do we remove conflictung events
        if REMOVE_CONFLICTS:
            if t.conflicts > 0:
                continue
        #do we remove potential interstrip events   
        if REMOVE_INTERSTRIP:
            if t.interstrip > 0:
                continue
        #use only events with 2 particle detections    
        if t.cnc != 2:
            continue
        
        #use only particles from chosen detectors, #particles can come from different detectors but  DET1<DET2  is always true so there is no need for reverse option
        if DETECTORS:
            if t.detector[0] not in DET1 or t.detector[1] not in DET2:
                continue
            
        for i in range(t.cnc):
            front_number = i*3
            dE_number = i*3+2
            #don't include particles that are not in the strip range that we want to use
            if  STRIP_RANGE:
                if dE_low <= t.adc[dE_number] <= dE_high and front_low <= t.adc[front_number] <= front_high:
                    pass
                else:
                    continue
                
            #don't include particles that are not in the strip pairs
            if t.adc[dE_number] not in list_pairs[t.adc[front_number]]:
                continue

        "Calculations"
        #choose particle types
        if t.ptype[0] == ptype1 and t.ptype[1] == ptype2: 
            counter_particles_coincident += 2
            
            mass = [] #field for masses of detected particles
            mass.append(M_DET1)
            mass.append(M_DET2)
            
            #Excitation energy and angle of the compound nucleus
            Ex_compound = calculate_Ex_compound(E_THRESHOLD, mass[0], mass[1], t.cnrg[0], t.cnrg[1], np.deg2rad(t.theta[0]),
                                        np.deg2rad(t.theta[1]),np.deg2rad(t.phi[0]),np.deg2rad(t.phi[1]))
            theta_compound,phi_compound = calculate_theta_compound(Ex_compound, E_THRESHOLD, mass[0], mass[1], t.cnrg[0], t.cnrg[1], np.deg2rad(t.theta[0]),
                                        np.deg2rad(t.theta[1]),np.deg2rad(t.phi[0]),np.deg2rad(t.phi[1]), M_COMPOUND)
                                        
            #compound excitation cut                            
            if Ex_compound_CUT:
                if not Ex_compound_MIN < Ex_compound < Ex_compound_MAX:
                    continue                          
            
            #romano plot values
            P_til, E_til = romano(Energy_PROJECTILE_HALFtarget , Mass_PROJECTILE, mass[0], mass[1], t.cnrg[0], t.cnrg[1], np.deg2rad(t.theta[0]),
                                        np.deg2rad(t.theta[1]),np.deg2rad(t.phi[0]),np.deg2rad(t.phi[1]))
            #Q value
            Q_value = 0. - E_til + P_til/M_UNDET_THIRD
            
            #Excitation energy of the third undetected particle
            Ex_third = Q0 - Q_value #energija pobuđenja treće nedetektirane čestice

            #Q value cut
            if Q_CUT:
                 if not Q_value_MIN < (Q0 - Q_value) < Q_value_MAX:
                    continue                           
            

            #romano cut
            if ROMANO_CUT:
                temp_counter = 0
                for i in range(len(cut_number_used)):
                    if cut_romano[cut_number_used[i]].IsInside(P_til,E_til):
                        temp_counter+=1
                if temp_counter == 0:
                    continue
            #Q vs theta compound cut
            if QthetaC_CUT:
                temp_counter = 0
                for i in range(len(cut_number_used)):
                    if cut_QthetaS[cut_number_used[i]].IsInside(theta_compound,Q0 - Q_value):
                        temp_counter+=1
                if temp_counter == 0:
                    continue            
            
            #excitation energies for possible single scenarios
            Ex1 = calculate_Ex_single(Q0_single1, Mass_PROJECTILE, mass[0], M_UNDET1, Energy_PROJECTILE_HALFtarget , t.cnrg[0], np.deg2rad(t.theta[0]))
            Ex2 = calculate_Ex_single(Q0_single2, Mass_PROJECTILE, mass[1], M_UNDET2, Energy_PROJECTILE_HALFtarget , t.cnrg[1], np.deg2rad(t.theta[1]))
            
            #do we write strip hits   
            if STRIP_HITS:
                for i in range(t.mult):
                    try:
                        strip_hits[t.adc[i]] += 1
                    except:#setting up the counter (dictionary member) for the strip at its first mention
                        strip_hits[t.adc[i]] = 1                 
                        
            counter_particles_filter +=2
            
            h1.Fill(t.nrg[0],t.nrg[2]), h1.Fill(t.nrg[3],t.nrg[5])             
            h2.Fill(t.cnrg[0],t.cnrg[1])
            h3.Fill(P_til,E_til) 
            h3colour.Fill(P_til,E_til) 
            h4.Fill(Ex_third)
            h5.Fill(Ex_third, theta_compound)
            h5colour.Fill(Ex_third, theta_compound)
            h6.Fill(Ex_compound)
            h7.Fill(Ex1, Ex_compound)
            h7colour.Fill(Ex1, Ex_compound)
            h8.Fill(Ex2, Ex_compound)
            h8colour.Fill(Ex2, Ex_compound)
            if SAME_PARTICLES:
                h_same.Add(h7, h8, 1.0, 1.0)
                h_same_colour.Add(h7colour, h8colour, 1.0, 1.0)
            h9.Fill(Ex2, Ex1)
            h9colour.Fill(Ex2, Ex1)
            h10.Fill(Ex1)
            h11.Fill(Ex2)
            if SAME_PARTICLES:
                h_same_projection.Add(h10, h11, 1.0, 1.0)         
        
        #ensuring that we do not count the same particle pairs twice or that we check the assimetry by skipping the reverse pairing of particle type and detector 
        if SAME_PARTICLES: 
            continue

        # if ASSYMETRY_CHECK:
        #     continue
        
        #reverse order of particle types
        if t.ptype[0] == ptype2 and t.ptype[1] == ptype1:
            counter_particles_coincident += 2
            
            mass = [] #field for masses of detected particles
            mass.append(M_DET2) #inverting the ordering of masses
            mass.append(M_DET1)
  
            #Excitation energy and angle of the compound nucleus
            Ex_compound = calculate_Ex_compound(E_THRESHOLD, mass[0], mass[1], t.cnrg[0], t.cnrg[1], np.deg2rad(t.theta[0]),
                                        np.deg2rad(t.theta[1]),np.deg2rad(t.phi[0]),np.deg2rad(t.phi[1]))
            theta_compound,phi_compound = calculate_theta_compound(Ex_compound, E_THRESHOLD, mass[0], mass[1], t.cnrg[0], t.cnrg[1], np.deg2rad(t.theta[0]),
                                        np.deg2rad(t.theta[1]),np.deg2rad(t.phi[0]),np.deg2rad(t.phi[1]), M_COMPOUND)
                                        

            #compound excitation cut                            
            if Ex_compound_CUT:
                if not Ex_compound_MIN < Ex_compound < Ex_compound_MAX:
                    continue  

            
            #romano plot values
            P_til, E_til = romano(Energy_PROJECTILE_HALFtarget , Mass_PROJECTILE, mass[0], mass[1], t.cnrg[0], t.cnrg[1], np.deg2rad(t.theta[0]),
                                        np.deg2rad(t.theta[1]),np.deg2rad(t.phi[0]),np.deg2rad(t.phi[1]))
            #Q value
            Q_value = 0. - E_til + P_til/M_UNDET_THIRD
            
            #Excitation energy of the third undetected particle
            Ex_third = Q0 - Q_value #energija pobuđenja treće nedetektirane čestice

            #Q value cut
            if Q_CUT:
                 if not Q_value_MIN < (Q0 - Q_value) < Q_value_MAX:
                    continue                             
            

            #romano cut
            if ROMANO_CUT:
                temp_counter = 0
                for i in range(len(cut_number_used)):
                    if cut_romano[cut_number_used[i]].IsInside(P_til,E_til):
                        temp_counter+=1
                if temp_counter == 0:
                    continue
            
            if QthetaC_CUT:
                temp_counter = 0
                for i in range(len(cut_number_used)):
                    if cut_QthetaS[cut_number_used[i]].IsInside(theta_compound,Q0 - Q_value):
                        temp_counter+=1
                if temp_counter == 0:
                    continue 
                
                
            #excitation energies for possible single scenarios
            Ex1 = calculate_Ex_single(Q0_single1, Mass_PROJECTILE, mass[0], M_UNDET1, Energy_PROJECTILE_HALFtarget , t.cnrg[0], np.deg2rad(t.theta[0]))
            Ex2 = calculate_Ex_single(Q0_single2, Mass_PROJECTILE, mass[1], M_UNDET2, Energy_PROJECTILE_HALFtarget , t.cnrg[1], np.deg2rad(t.theta[1]))

            #do we write strip hits   
            if STRIP_HITS:
                for i in range(t.mult):
                    try:
                        strip_hits[t.adc[i]] += 1
                    except:#setting up the counter (dictionary member) for the strip at its first mention
                        strip_hits[t.adc[i]] = 1                               
            
            counter_particles_filter +=2
            
            h1.Fill(t.nrg[0],t.nrg[2]), h1.Fill(t.nrg[3],t.nrg[5])             
            h2.Fill(t.cnrg[0],t.cnrg[1])
            h3.Fill(P_til,E_til) 
            h3colour.Fill(P_til,E_til) 
            h4.Fill(Ex_third)
            h5.Fill(Ex_third, theta_compound)
            h5colour.Fill(Ex_third, theta_compound)
            h6.Fill(Ex_compound)
            h7.Fill(Ex1, Ex_compound)
            h7colour.Fill(Ex1, Ex_compound)
            h8.Fill(Ex2, Ex_compound)
            h8colour.Fill(Ex2, Ex_compound)
            h9.Fill(Ex2, Ex1)
            h9colour.Fill(Ex2, Ex1)
            h10.Fill(Ex1)
            h11.Fill(Ex2)
                
#when there are multiple runs, we calculate the total elapsed time
total_time = round(time()-start_time)
total_time = timedelta(0,total_time)  
  
    
g = TFile(out_name,"recreate")  
h1.Write()
h2.Write()
h3.Write()
h3colour.Write()
h4.Write()
h5.Write()
h5colour.Write()
h6.Write()
h7.Write()
h7colour.Write()
h8.Write()
h8colour.Write()
if SAME_PARTICLES:
    h_same.Write()
    h_same_colour.Write()
h9.Write()
h9colour.Write()
h10.Write()
h11.Write()
if SAME_PARTICLES:
    h_same_projection.Write()
g.Close()      
            
print("Detections: %s+%s, Filters used: %s, run18-30"%(Name_DET1, Name_DET2, output_suffix) )
print("Number of particles in the used runs:", counter_particles)
print("   Number of particles in pairs of the chosen type:", counter_particles_coincident)
print("      In chosen filters %s there are:"%output_suffix, counter_particles_filter  )
if STRIP_HITS:
    print("Number of strip hits:", strip_hits)
print("The total time of execution", total_time)
print("END")
