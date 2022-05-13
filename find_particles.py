# -*- coding: utf-8 -*-
"""
Created on Mon Apr 29 04:49:04 2019

@author: Deni

This is a script that takes the information about the number and the strip position of the detector hits in an event and deduces which combination of 
those hits belong to certain particles. A necessery step to be able to operate with those particles further on, adding the appropriate angles,
correcting the energies and calculating the final energy spectra. 

"""

""" import predefined functions """
from ROOT import TFile, TTree, TH2S, gPad, TCanvas
from numpy import zeros, matrix, shape
from itertools import product
from time import time
from datetime import timedelta
from sys import argv


""" 
Define program-wide info
    !ADJUST!: N_RUN, CALIB,NRG_TOLERANCE, AMPL_TOLERANCE, 
        HAS_dE (True,False), HAS_ENERGY (True,False), DEBUG(True,False), PARTICLE_ENTRY(True,False), EVENT_ENTRY(True,False)
        All filters (True,False), MULTIPLICITY_RULE(True,False), S_SUFFIX2
        INTERSTRIP_MARKING (True,False), CONFLICT_MARKING (True,False)        
"""



S_RUN = 'run26' if len(argv) < 2 else argv[1]  # S_RUN is the first argument of the scripts, 'run26' if not given
CALIB = 'TMIN2Be' # the name of the calibration
S_RUN = f'{S_RUN}_NRG({CALIB})'
TREE_NAME = 'tree' #usually we used names "tree" or "T"


DEBUG = False # do we want an output check of certain parts of the code

HAS_ENERGY = True # do we have already the energies of the hits, or do we only have amplitudes

if HAS_ENERGY: NRG_TOLERANCE = 0.05 # the difference between the FRONT and BACK energies that we allow while constituting a particle
else: AMPL_TOLERANCE = 0.3 # the difference between the FRONT and BACK amplitudes that we allow while constituting a particle

#Define the output tree entry structure
PARTICLE_ENTRY = False  # will an entry in the outgoing file be a single particle or an entire event?
#maximum number of particles that we write into a single root entry
if PARTICLE_ENTRY: 
    MAX_PARTICLES_ENTRY = 1 
    S_SUFFIX = 'P=E' # suffix to add to the output files and/or histograms
else:
    MAX_PARTICLES_ENTRY = 30 # any number that is higher than the expected maximum number of particles in an experimental event, 10 is enough when the conflict filter is on, more if it is not
    S_SUFFIX = 'E=E' # suffix to add to the output files and/or histograms

HAS_dE = True # is there a deltaE detector in the run

N_ADC = 192 # total number of detector strips/channels
MAX_mult = 150 # maximum number of hits in the event that we count as a physical event, designed to discard the PULSAR events
if HAS_dE:
    MIN_mult = 3 # minimum number of hits in the event needed to define a physical event, dE, FRONT and BACK info needed to constitue a single particle
else: 
    MIN_mult = 2 # minimum number of hits in the event needed to define a physical event, FRONT and BACK info needed to constitue a single particle

    
#Turning filters ON and OFF
FILTER_1 = True    #removes PULSAR and NOISE, meaning events with too many or too few hits
FILTER_Interstrip_1 = False #removes possible interstrip events by removing neighbouring dE hits prior to combination formation
FILTER_2 = True #do the hits in a combination belong to the same detector, meaning is the combination geometrically possible
FILTER_3 = True #do the front and back hits of a combination match in amplitude/energy
FILTER_4 = True #are the hits in a combination geometrically close, meaning do they belong to nearby strips (ADCs)
FILTER_Interstrip_2 = False #removes possible interstrip events by removing combinations with neighbouring dE hits
FILTER_5 = False #removing combinations in an event which share a front/back/dE strip(ADC)
MULTIPLICITY_RULE = True # do we compare multiplicity od the event and the number of combinations prior to excecuting FILTER5

if not FILTER_Interstrip_2: # if we turn OF this filter, do we mark the particles that this filter would remove
    INTERSTRIP_MARKING = True
if not FILTER_5:  # if we turn OF this filter, do we mark the particles that this filter would remove
    CONFLICT_MARKING = True

S_SUFFIX2 = 'F1-4' #suffix to denote the filters we have used


#a matric containing the detector number and the associated front/back/dE strips (ADCs) 
#MATRIX STRUCTURE: 1st front strip, last front strip, 1st back, last back, 1st delta E, last dE, detector number
DETECTORS_STRIPS_MATRIX = matrix([
                                  [1,16,17,32,129,144,1], #detA 
                                  [33,48,49,64,145,160,2], #detB
                                  [65,80,81,96,161,176,3], #detC
                                  [97,112,113,128,177,192,4], #detD                  
                                 ])
                    
                    
OUT_RUN = f'{S_RUN}_NRG({CALIB})_particles({S_SUFFIX}, {S_SUFFIX2})'

                    
                    
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
Functions for hit categorization:
    get_detector_adcType: determines the detector number and type (front/back/dE)
    categorize_hits: separating hits into lists according to strip(adc) type (front/back/dE), making a dictionary with all the information (key=strip)
    get_all_potential_combinations: making all possible (front,back,dE) combinations as a product of strip(adc) hits in an event
"""      
              
# getting detector number and type from the matrix, marking adc type: 0 = front, 1 = back, 2 = dE
def get_detector_adcType(adc):
    for i in range(shape(DETECTORS_STRIPS_MATRIX)[0]):
        if adc >= DETECTORS_STRIPS_MATRIX[i,0] and adc <= DETECTORS_STRIPS_MATRIX[i,1]:
            return [0, DETECTORS_STRIPS_MATRIX[i,6]]
        elif adc >= DETECTORS_STRIPS_MATRIX[i,2] and adc <= DETECTORS_STRIPS_MATRIX[i,3]:
            return [1, DETECTORS_STRIPS_MATRIX[i,6]]
        elif adc >= DETECTORS_STRIPS_MATRIX[i,4] and adc <= DETECTORS_STRIPS_MATRIX[i,5]:
            return [2, DETECTORS_STRIPS_MATRIX[i,6]]
    return [-1, -1]
    

# categorizing hits in an event acording to type, making a dictionary with the information, strip(adc) is the dict key because every strip can only be hit once during an event
def categorize_hits(t):           
    hits_dict = {}
    lst_front = []
    lst_back = []
    lst_dE = []
    
    for j in range(t.mult):
        tmp = {}
        tmp['ampl'] = t.ampl[j]
        if HAS_ENERGY:
            tmp['nrg']=t.nrg[j]
            
        adc_type, adc_detector = get_detector_adcType(t.adc[j])
        tmp['detector'] = adc_detector   
        
        hits_dict[t.adc[j]] = tmp 
        
        if adc_type == 0:
            lst_front.append(t.adc[j])
        elif adc_type == 1:
            lst_back.append(t.adc[j])
        else:
            lst_dE.append(t.adc[j])
            
    return lst_front, lst_back, lst_dE, hits_dict


# making all possible (front,back,dE) combinations as a product of strip(adc) hits in an event
def get_all_potential_combinations(lst_front, lst_back, lst_dE):
    if HAS_dE:
        combinations = product(lst_front,lst_back,lst_dE)
    else:
        combinations = product(lst_front,lst_back)

    return combinations
    
    
"""
Functions for filters
"""

#FILTER1  
#removes event with too many and too few hits
def is_mult_ok(mult):    
    return MIN_mult <= mult <= MAX_mult
    
    
#FILTER_Interstrip_1  
#removes possible interstrip events by removing neighbouring dE hits prior to combination formation
def is_neighbouring_dE_1(list_dE, hits_dictionary):
    
    good_dE = list_dE[:] #initially all hits are seen as good hits
    bad_dE = []
    
    #removing neighbouring dE hits from the starting list
    for i in range(len(list_dE)-1):
        a = list_dE[i]
        a_det = hits_dictionary[a]['detector']
        for j in range(i+1,len(list_dE)):
            b = list_dE[j]
            b_det = hits_dictionary[b]['detector']
            
            #if they are neighbouring and on the same detector, we remove both of them
            if abs(a-b) == 1 and a_det == b_det:
                
                if list_dE[i] in good_dE:
                    bad_dE.append(list_dE[i])
                    good_dE.remove(list_dE[i])
                if list_dE[j] in good_dE:
                    bad_dE.append(list_dE[j])
                    good_dE.remove(list_dE[j])
                    
    return good_dE, bad_dE
    
    
#FILTER2
# do the hits in a combination belong to the same detector, meaning is the combination geometrically possible (same row of the above defined matrix)
def is_combination_possible_1(combination):
    matrix = DETECTORS_STRIPS_MATRIX
    front = combination[0]
    back = combination[1]
    if HAS_dE:
        dE = combination[2]
    for i in range(shape(matrix)[0]):
        if HAS_dE:
            if front >= matrix[i,0] and front <= matrix[i,1] and back >= matrix[i,2] and back <= matrix[i,3] and dE >= matrix[i,4] and dE <= matrix[i,5] :
                return True
        elif front >= matrix[i,0] and front <= matrix[i,1] and back >= matrix[i,2] and back <= matrix[i,3]:
            return True
    return False
    
#FILTER3
#do the front and back hits of a combination match in amplitude/energy
    
def front_and_back_match(combination, hits_dict): #taking either amplitude or energy as a measure (key), and the appropriate limit for matching front and back hits
    if HAS_ENERGY: 
        key = 'nrg'
        limit = NRG_TOLERANCE
    else:
        key = 'ampl'
        limit = AMPL_TOLERANCE
        
    return is_difference_ok(combination, hits_dict, key, limit)
    
# calculating the difference in front and back hits                   
def is_difference_ok(combination, dictionary, key, limit):
    front = dictionary[combination[0]][key]
    back = dictionary[combination[1]][key]
    diff = 2.*abs(float(front-back))/(front+back)
    #diff = abs(float(front-back)/front)
    if diff <= limit:
        return True
    if DEBUG:
        print("FILTER3, difference of f and b hit:", dictionary[combination[0]], dictionary[combination[1]], diff)
    return False 


    
#FILTER4  
#are the hits in a combination geometrically justified, meaning do they belong to nearby strips (ADCs)
# geometric neighbourhood determined by dE-E matching, currently wider area around that information used for simpler code and more leniency 
def is_combination_possible_2(combination):
    matrix = DETECTORS_STRIPS_MATRIX
    front = combination[0]
    dE = combination[2]
    
    if front >= matrix[0,0] and front <= matrix[0,1]  and dE >= front + 127 and dE <= front + 129 : 
        return True
        
    elif front >= matrix[1,0] and front <= matrix[1,1]  and dE >= front + 111 and dE <= front + 113 :
        return True
        
    elif front >= matrix[2,0] and front <= matrix[2,1]  and dE >= front + 95 and dE <= front + 98 :
        return True
    
    elif front >= matrix[3,0] and front <= matrix[3,1]  and dE >= front + 78 and dE <= front + 81 :
        return True
    
    return False


#FILTER_Interstrip_2  
#removes possible interstrip events by removing combinations with neighbouring dE hits
def is_neighbouring_dE_2(candidates, hits_dictionary):

    no_interstrip = candidates[:] #initially all combinations are seen as good combinations
    interstrip = []
    
    #removing combinations with neighbouring dE hits
    for i in range(len(candidates)-1):
        a = candidates[i]
        a_dE_det = hits_dictionary[a[2]]['detector']
        for j in range(i+1,len(candidates)):
            b = candidates[j]
            b_dE_det = hits_dictionary[b[2]]['detector']
            
            #if dE are neighbouring and on the same detector, we remove both of them
            if abs(a[2]-b[2]) == 1 and a_dE_det == b_dE_det:  
                
                if candidates[i] in no_interstrip:
                    interstrip.append(candidates[i])
                    no_interstrip.remove(candidates[i])
                if candidates[j] in no_interstrip:
                    interstrip.append(candidates[j])
                    no_interstrip.remove(candidates[j])
                    
    return no_interstrip, interstrip


#FILTER5    
#removing combinations in an event which share a front/back/dE strip(ADC)
def filter_out_conflicts(candidates):
      
    good_combinations = candidates[:] #initially all combinations are seen as good combinations
    bad_combinations = []
    
                        
    # removing those that have the same front/back/dE 
    for i in range(len(candidates)-1):
        a = set(candidates[i])
        for j in range(i+1,len(candidates)):
            b = set(candidates[j])
            
            # if two combinations share any one of the tri strips, we remove both of them
            if a & b:
                if candidates[i] in good_combinations:
                    bad_combinations.append(candidates[i])
                    good_combinations.remove(candidates[i])
                if candidates[j] in good_combinations:
                    bad_combinations.append(candidates[j])
                    good_combinations.remove(candidates[j])
        
    return good_combinations, bad_combinations
    
    
    
""" EXCECUTION """ 

"""open the file and get the tree"""
myfile = TFile(f'./nrg/{S_RUN}.root') #open
if DEBUG:
    myfile.ls() #check the file contents
t = myfile.Get(TREE_NAME)  #get the tree
if DEBUG:
    t.Print() #check the tree contents
    
    
"""creating the output root file; defining name, output tree, output tree variables, output tree structure """
out_file = TFile(f'./particles/{OUT_RUN}.root', "recreate") # output root file
out_t = TTree('tree', f'{OUT_RUN}') # output tree


# output tree variables
# variables that exist in the starting tree
mult = zeros(1, dtype='short') # multiplicity of an event - number of hits, value of type Short
adc = zeros(MIN_mult * MAX_PARTICLES_ENTRY, dtype='short') # channel number, field of values of type Short
ampl = zeros(MIN_mult * MAX_PARTICLES_ENTRY, dtype='short') # amplitude, field of values of type Short
if HAS_ENERGY:
    nrg = zeros(MIN_mult * MAX_PARTICLES_ENTRY, dtype='float') # energy, field of values of type Float
#new variables
event = zeros(1, dtype= 'int')# the number of events in the starting run, value of type Int
detector = zeros(MAX_PARTICLES_ENTRY, dtype='short')# which detector the particle went into; 1,2,3 or 4, value of type Short
cnc = zeros(1, dtype='short') # the number of particles in an event (coincidences), value of type Short
if CONFLICT_MARKING:
    conflicts = zeros(1, dtype='short') # the number of conflicting particles in an event, value of type Short
if INTERSTRIP_MARKING:
    interstrip = zeros(1, dtype='short') # the number of potential interstrip particles in an event, value of type Short

# output tree structure; adding new branches
out_t.Branch('event', event, 'event/I') # event number branch
out_t.Branch('cnc', cnc, 'cnc/S') # coincidences branch    
out_t.Branch('mult',mult,'mult/S') # multiplicity branch
if PARTICLE_ENTRY:  
    out_t.Branch('adc',adc,'adc[3]/S') # channel number branch, field of lenght 3
    out_t.Branch('ampl',ampl,'ampl[3]/S') # amplitude branch, field of lenght 3
    if  HAS_ENERGY:
        out_t.Branch('nrg', nrg, 'nrg[3]/D') # energy branch, field of lenght 3
    out_t.Branch('detector', detector, 'detector/S') # detector branch
else:
    out_t.Branch('adc',adc,'adc[mult]/S') # channel number branch, field of lenght "mult"
    out_t.Branch('ampl',ampl,'ampl[mult]/S') # amplitude branch, field of lenght "mult"
    if  HAS_ENERGY:
        out_t.Branch('nrg', nrg, 'nrg[mult]/D') # energy branch, field of lenght "mult"
    out_t.Branch('detector', detector, 'detector[cnc]/S') # detector branch, field of lenght "cnc"

if CONFLICT_MARKING:
    out_t.Branch('conflicts', conflicts, 'conflicts/S') #conflict branch
if INTERSTRIP_MARKING:
    out_t.Branch('interstrip', interstrip, 'interstrip/S') #interstrip branch    

"""  
definition of histograms
"""

#dE-E histograms that record the filtered out combinations
if HAS_dE: #amplitude histograms
    histogram2 = TH2S('filter2', "deltaE-E", 3000, 0, 3000, 3000, 0, 3000) # filtered out combinations by filter2
    histogram2.SetOption("COLZ") # contrast by number 
    histogram3 = TH2S('filter3', "deltaE-E", 3000, 0, 3000, 3000, 0, 3000) # filtered out combinations by filter3
    histogram3.SetOption("COLZ") # contrast by number 
    histogram4 = TH2S('filter4', "deltaE-E", 3000, 0, 3000, 3000, 0, 3000) # filtered out combinations by filter4
    histogram4.SetOption("COLZ") # contrast by number 
    histogram8 = TH2S('filter Interstrip2', "deltaE-E", 3000, 0, 3000, 3000, 0, 3000) # filtered out combinations by Interstrip2 filter
    histogram8.SetOption("COLZ") # contrast by number        
    histogram5 = TH2S('filter5', "deltaE-E", 3000, 0, 3000, 3000, 0, 3000) # filtered out combinations by filter5
    histogram5.SetOption("COLZ") # contrast by number 
    histogram6 = TH2S('marked conflicts', "deltaE-E", 3000, 0, 3000, 3000, 0, 3000) # for conflicting combinations which are marked
    histogram6.SetOption("COLZ") # contrast by number 
    histogram7 = TH2S('marked interstrip', "deltaE-E", 3000, 0, 3000, 3000, 0, 3000) # for potential interstrip combinations wich are marked
    histogram7.SetOption("COLZ") # contrast by number 
    if HAS_ENERGY: #energy histograms
        histogram2_nrg = TH2S('filter2_NRG', "deltaE-E", 6000, 0, 60, 6000, 0, 60) # filtered out combinations by filter2
        histogram2_nrg.SetOption("COLZ") # contrast by number 
        histogram3_nrg = TH2S('filter3_NRG', "deltaE-E", 6000, 0, 60, 6000, 0, 60) # filtered out combinations by filter3
        histogram3_nrg.SetOption("COLZ") # contrast by number 
        histogram4_nrg = TH2S('filter4_NRG', "deltaE-E", 6000, 0, 60, 6000, 0, 60) # filtered out combinations by filter4
        histogram4_nrg.SetOption("COLZ") # contrast by number 
        histogram8_nrg = TH2S('filter Interstrip2_NRG', "deltaE-E", 6000, 0, 60, 6000, 0, 60) # filtered out combinations by Interstrip2 filter
        histogram8_nrg.SetOption("COLZ") # contrast by number 
        histogram5_nrg = TH2S('filter5_NRG', "deltaE-E", 6000, 0, 60, 6000, 0, 60) # filtered out combinations by filter5
        histogram5_nrg.SetOption("COLZ") # contrast by number 
        histogram6_nrg = TH2S('marked conflicts_NRG', "deltaE-E", 6000, 0, 60, 6000, 0, 60) # for conflicting combinations which are marked
        histogram6_nrg.SetOption("COLZ") # contrast by number 
        histogram7_nrg = TH2S('marked interstrip_NRG', "deltaE-E", 6000, 0, 60, 6000, 0, 60) # for potential interstrip combinations wich are marked
        histogram7_nrg.SetOption("COLZ") # contrast by number                 

# a histogram for vizualizing the front and back match with borders that define which combinations are filtered out by the "front and back match" filter3.        
if not HAS_ENERGY:
    histogram = TH2S('signal match check', "Eb-Ef", 3000, 0, 3000, 3000, 0, 3000) 
    histogram_border = TH2S('hist6', "borders", 3000, 0, 3000, 3000, 0, 3000) # borders of front and back mathing
else:
    histogram = TH2S('signal match check', "Eb-Ef", 3000, 0, 60, 3000, 0, 60)
    histogram_border = TH2S('hist6', "borders", 3000, 0, 60, 3000, 0, 60) # borders of front and back mathing
histogram.SetOption("COLZ") # contrast by number

for i in range(3000): #border calculation
    x = x_E = y1 = y2 = 0.0     #.0 to avoid integer calculation

    x = i   #drawing a bunch of points, x is Ef, range of amplitudes 0-3000, range of energies 0-60
    
    if not HAS_ENERGY: #Borders are defined by inverting the expressions in "is_difference_ok" to a form E_b = a*E_f  
        a = AMPL_TOLERANCE/2
        y1 =  (1+a)/(1-a)* x  #uper border
        y2 =  (1-a)/(1+a)* x #lower border
        histogram_border.Fill(x, y1)
        histogram_border.Fill(x, y2)
    if  HAS_ENERGY:
        a = NRG_TOLERANCE/2
        x_E = x / 50.0 #scaled for energy histogram
        y1 =  (1+a)/(1-a) * x_E  #uper border
        y2 =  (1-a)/(1+a) * x_E #lower border
        histogram_border.Fill(x_E, y1)
        histogram_border.Fill(x_E, y2)


""" go trough the run data """
n_entries = t.GetEntriesFast() # get the  number of entries in the tree
if DEBUG:
    print("number of entries: ", n_entries)
    # n_entries = 30

start_time = time() #start the timer  

#initializing various counters
event_counter = 0 #counting the number of remaining events in the out run  
num_particles = 0 # number of detected particles
num_pulsar_noise = 0 # number of events with too few (noise) and too many (pulsar) hits
num_bad_signals = 0 # number of combinations whose front and back signals do not match
num_impossible_combinations = 0 # number of combinations whose strips(ADCs) belong to different detectors
num_geometricallyBad_combinations = 0 #number of combinations whose strips(ADCs) are on the same detector but geometrically far apart
num_conflicted_adcs = 0 # number of combinations that share a front/back/dE strip(ADC)
num_conflicted_adcs2 = 0 #additional check for conflict counting
num_interstrip_1 = 0 # number of neighbouring dE hits
num_interstrip_2 = 0 # number of combinations with neighbouring dE hits after going through the previous filters
num_coincidences = {} # dictionary for writing the number of coincidences

for i in range(n_entries): 
    if i % 50000 == 0: # print the timer info every 50000 entries
        print_time(start_time,i,n_entries)
            
    t.GetEntry(i) # get an entry(event)
    
    if DEBUG:
        print("NEXT EVENT:", event_counter, "multiplicity:", t.mult)

     # 1. FILTER        
    # immidiately remove events with too many and too few particles
    if FILTER_1: 
        if not is_mult_ok(t.mult):
            num_pulsar_noise += 1
            if DEBUG:
                print("FILTER1, number of hits out of range:", mult)
            continue
        
    # categorizing hits in an event
    lst_front, lst_back, lst_dE, hits_dict = categorize_hits(t)
    
    if DEBUG:
        print("Hits", hits_dict, "f:", lst_front, "b:", lst_back, "dE:", lst_dE)

     # FILTER Interstrip 1
    #removes possible interstrip events by removing neighbouring dE hits prior to combination formation
    if FILTER_Interstrip_1:
        bad_dE = []
        good_dE = []
        
        good_dE, bad_dE = is_neighbouring_dE_1(lst_dE, hits_dict)
        
        lst_dE = good_dE
        num_interstrip_1 += len(bad_dE)



    # making all possible (front,back,dE) combinations as a product of strip(adc) hits in an event
    combinations = get_all_potential_combinations(lst_front, lst_back, lst_dE)
    candidates = []  #a variable to save intermediate good combinations
    

    # going through all combinations and eliminating the "bad ones
    for combination in combinations:
        if DEBUG:
            print("Reviewing combination:", combination)
        
        
         # 2. FILTER
        # eliminating combinations with strips(ADCs) from different detectors
        if FILTER_2: 
            if not is_combination_possible_1(combination):                
                num_impossible_combinations += 1
                if HAS_dE:
                    histogram2.Fill(hits_dict[combination[0]]['ampl'],hits_dict[combination[2]]['ampl']) #put filtered out combinations in an dE-E histogram 
                    if HAS_ENERGY:
                       histogram2_nrg.Fill(hits_dict[combination[0]]['nrg'],hits_dict[combination[2]]['nrg']) #put filtered out combinations in an dE-E histogram  
                if DEBUG:
                    print("FILTER2, geometrically impossible combination:", combination)
                continue


         # 3. FILTER
        # front and back match
        if FILTER_3:           
            if not front_and_back_match(combination, hits_dict):
                num_bad_signals += 1 
                if not HAS_ENERGY:
                    histogram.Fill(hits_dict[combination[0]]['ampl'],hits_dict[combination[1]]['ampl']) #put filtered out combinations in an Ef-Eb histogram
                if  HAS_ENERGY:
                    histogram.Fill(hits_dict[combination[0]]['nrg'],hits_dict[combination[1]]['nrg']) #put filtered out combinations in an Ef-Eb histogram
                if HAS_dE:
                    histogram3.Fill(hits_dict[combination[0]]['ampl'],hits_dict[combination[2]]['ampl']) #put filtered out combinations in an dE-E histogram
                    if HAS_ENERGY:
                        histogram3_nrg.Fill(hits_dict[combination[0]]['nrg'],hits_dict[combination[2]]['nrg']) #put filtered out combinations in an dE-E histogram                
                continue  
        
        
        # putting all remaining combinations after the first 3 filters in the same Ef-Eb histogram
        if not HAS_ENERGY:
            histogram.Fill(hits_dict[combination[0]]['ampl'],hits_dict[combination[1]]['ampl']) 
        if HAS_ENERGY:
            histogram.Fill(hits_dict[combination[0]]['nrg'],hits_dict[combination[1]]['nrg']) 
        
         # 4. FILTER
        # geometrical neighbourhood
        if FILTER_4:
            if HAS_dE:
                if not is_combination_possible_2(combination):
                    num_geometricallyBad_combinations += 1                     
                    histogram4.Fill(hits_dict[combination[0]]['ampl'],hits_dict[combination[2]]['ampl']) #put filtered out combinations in an dE-E histogram
                    if HAS_ENERGY:
                        histogram4_nrg.Fill(hits_dict[combination[0]]['nrg'],hits_dict[combination[2]]['nrg']) #put filtered out combinations in an dE-E histogram
                
                    if DEBUG:
                        print("FILTER4, geometrically improbable combination", combination)
                    
                    continue
            
            

        #at this point the combination has passed all the filters except for the Interstrip2 and Filter5(conflict)
        candidates.append(combination)
        if DEBUG:
            print("A particle candidate (combination)", combination, "passed filters 1-4")
       

    if DEBUG:
        print("All particle candidates (combinations) in an event", candidates)


    num_candidates = len(candidates) #number of particle candidates at this point
    
    # if there are no particle candidates, go to the next event       
    if num_candidates < 1:
        if DEBUG:
            print("No good combinations after the first 4 filters")
        continue
    
    
     # FILTER Interstrip 2
    if FILTER_Interstrip_2:
        if num_candidates >= 2:
            #if there are combinations with neighbouring dE hits, we remove them
            no_interstrip_combinations, interstrip_combinations = is_neighbouring_dE_2(candidates, hits_dict)
            
            candidates = []
            candidates = no_interstrip_combinations[:]
            num_candidates = len(candidates)
            num_interstrip_2 += len(interstrip_combinations)
            
            for combination in interstrip_combinations:
                if HAS_dE:
                    histogram8.Fill(hits_dict[combination[0]]['ampl'],hits_dict[combination[2]]['ampl'])
                    if HAS_ENERGY:
                        histogram8_nrg.Fill(hits_dict[combination[0]]['nrg'],hits_dict[combination[2]]['nrg'])
                if DEBUG:
                    print("FILTER Interstrip 2, possible interstrip combinations", interstrip_combinations)
     
            if num_candidates < 1:
                if DEBUG:
                    print("No good combinations after Interstrip_2 filter")
                continue
    # if we are not filtering out interstrip particles but are marking them before letting them through    
    if (not FILTER_Interstrip_2) and INTERSTRIP_MARKING: 
        num_interstrip_combinations = 0               
        if num_candidates >= 2: 
            no_interstrip_combinations, interstrip_combinations = is_neighbouring_dE_2(candidates, hits_dict)
            num_interstrip_combinations = len(interstrip_combinations)
            num_interstrip_2 += len(interstrip_combinations)
            
            
            for combination in interstrip_combinations:
                if HAS_dE:
                    histogram7.Fill(hits_dict[combination[0]]['ampl'],hits_dict[combination[2]]['ampl'])
                    if HAS_ENERGY:
                        histogram7_nrg.Fill(hits_dict[combination[0]]['nrg'],hits_dict[combination[2]]['nrg'])
            if DEBUG:
                print("Marked potential interstrip combinations ", interstrip_combinations)
            
    
    
     # 5. FILTER
    # removing conflicts if there is more than one combination
    if FILTER_5:         
        if num_candidates == 1:
                good_combinations = candidates[:]  
                bad_combinations = []
        else:
            if MULTIPLICITY_RULE: #using the multiplicity rule
                if t.mult < num_candidates*3: 
                    good_combinations, bad_combinations = filter_out_conflicts(candidates)
                    
                
                elif t.mult >= num_candidates*3: 
                        good_combinations = candidates[:]  
                        bad_combinations = []
                        
            if not MULTIPLICITY_RULE: 
                if num_candidates>=2:
                    good_combinations, bad_combinations = filter_out_conflicts(candidates)
                    
            for combination in bad_combinations:
                if HAS_dE:
                    histogram5.Fill(hits_dict[combination[0]]['ampl'],hits_dict[combination[2]]['ampl'])
                    if HAS_ENERGY:
                        histogram5_nrg.Fill(hits_dict[combination[0]]['nrg'],hits_dict[combination[2]]['nrg'])
                if DEBUG:
                    print("FILTER5, conflict of combination's strips(ADCs)", bad_combinations)
        
            
    if not FILTER_5: # if we are not filtering out interstrip particles  
        if not CONFLICT_MARKING:
            good_combinations = candidates[:]  
            bad_combinations = []
            
        if CONFLICT_MARKING: #but we are marking them before letting them through
            if num_candidates == 1:
                    good_combinations = candidates[:]  
                    bad_combinations = []

            else: 
                good_combinations, bad_combinations = filter_out_conflicts(candidates)
                
                for combination in bad_combinations:
                    if HAS_dE:
                        histogram6.Fill(hits_dict[combination[0]]['ampl'],hits_dict[combination[2]]['ampl'])
                        if HAS_ENERGY:
                            histogram6_nrg.Fill(hits_dict[combination[0]]['nrg'],hits_dict[combination[2]]['nrg'])
                if DEBUG:
                    print("Marked conflicts", bad_combinations)
            


    num_good_combinations = len(good_combinations)
    num_conflicted_adcs += num_candidates - num_good_combinations
    num_bad_combinations = len(bad_combinations) #sanity check for number of conflicts, num_conflicted_adcs and num_conflicted_adcs2 must be equal
    num_conflicted_adcs2 += num_bad_combinations
        
    
    
    """ Filters done, writing down the remaining combinations """
    # if there are no particle candidates, go to the next event 
    if not CONFLICT_MARKING:
        if num_good_combinations < 1:
            if DEBUG:
                print("No good combinations (particles)")
            continue

    if DEBUG:
        print("FILTERS DONE")
        print("Good combinations (particles)", good_combinations)
        print(" Conflicting combinations", bad_combinations)
    
    #if we are marking and writing down conflicts, then the candidates that were present prior to the Filter5 are denoted as good combinations
    if CONFLICT_MARKING:
        num_good_combinations = num_candidates
        good_combinations = candidates [:]

    # counting the number of events that have certains numbers of coincidences (particles)
    try:
        num_coincidences[num_good_combinations] += 1
    except:
        num_coincidences[num_good_combinations] = 1 #setting up the dictionary member (counter) for a given coincidence number
        

    #raising the number of "valuable events" by 1, meaning events with defined particles 
    event_counter += 1
    
    #particle and hits counters for writing output data in the EVENT=ENTRY mode 
    particles_counter = 0
    hits_counter = 0
    
    
    # variables that are attributed to all of the particles in an event: event, mult, cnc
    event[0] = event_counter #the cardinal number of the event 
    mult[0] =  num_good_combinations* MIN_mult #the number of hits in an event
    cnc[0] = num_good_combinations #the number of coincidences in an event
    if CONFLICT_MARKING:
        conflicts[0] = num_bad_combinations #the number of conflicts in an event
    if INTERSTRIP_MARKING:
        interstrip[0] = num_interstrip_combinations #the number of possible interstrip particles in an event
        
    # going through combinations(particles) one by one
    for combination in good_combinations:          
        # variables that are different for different particles in an event: detector, adc, ampl, nrg                                  
        
        if PARTICLE_ENTRY:
            detector[0] = hits_dict[combination[0]]['detector'] #detector number for the particles in the event
        else:
            detector[0 + particles_counter] = hits_dict[combination[0]]['detector'] #detector number for the particles in the event
                                
        for j in range(3):
            adc[j + hits_counter] = combination[j]
            ampl[j + hits_counter] = hits_dict[combination[j]]['ampl']
            if HAS_ENERGY:
                nrg[j + hits_counter] = hits_dict[combination[j]]['nrg']
                            
        num_particles += 1
        if PARTICLE_ENTRY:
            out_t.Fill() # filling the out tree inside the particle loop; particle is the output entry
        else: #adding extra information in the adc,ampl,nrg and detector fields for subsequent particles
            particles_counter += 1 #detector[0] for particle 1, detector[1] for particle 2, and so on
            hits_counter += 3 #adc,ampl,nrg[0,1,2] for particle1, adc,ampl,nrg[3,4,5] for particle2, and so on
    
    if not PARTICLE_ENTRY:
        out_t.Fill() # filling the out tree outside the particle loop; event is the output entry as in the raw root files
    
    
    
out_file.Write()


# ispis rezultata
print("Starting number of events(entries)", n_entries)
if PARTICLE_ENTRY:
    print("Final number of particles (entries)", num_particles)
    print("Final number of events", event_counter)
else:
    print("Final number of particles", num_particles)
    print("Final number of events(entries)", event_counter)
print("Coincidences", num_coincidences)
print("Filtered out dE because of interstrip:", num_interstrip_1)
print("Filtered out pulsar and noise events", num_pulsar_noise)
print("Filtered out combinations: ")
print("  impossible combinations", num_impossible_combinations)
print("  front and back match", num_bad_signals)
print("  geometrically improbable", num_geometricallyBad_combinations)
print("  potential dE interstrip", num_interstrip_2)
print("  strip(ADC) conflicts", num_conflicted_adcs, num_conflicted_adcs2)
if INTERSTRIP_MARKING:
    print("Potential interstrip particles marked but not filtered out.")
if CONFLICT_MARKING:
    print("Conflicts marked but not filtered out.")
    if INTERSTRIP_MARKING:
        print("The conflict number contains the potential interstrip number.")





""" writing histograms inside the output root file and drawing an Ef-Eb histogram onto a canvas"""
print("Writing and drawing histograms")
if HAS_dE: #dE-E histograms of filtered out particles
    histogram2.Write()
    histogram3.Write(), histogram4.Write(), histogram8.Write(), histogram5.Write(), histogram6.Write(), histogram7.Write()
    if HAS_ENERGY:
        histogram2_nrg.Write(), histogram3_nrg.Write(), histogram4_nrg.Write(), histogram8_nrg.Write(), histogram5_nrg.Write(), histogram6_nrg.Write(), histogram7_nrg.Write()

#writing Ef-Eb histogram into the output root file, no visual differentiation between borders and combinations
histogram.Add(histogram, histogram_border, 1.0, 1.0) 
histogram.Write()            
    
# drawing the Ef-Eb histogram onto a canvas with visual differentiation
can = TCanvas('hist','Histogram',1980,1080)
histogram.Draw()
histogram_border.Draw("SAME PMC")
gPad.Update()
can.SaveAs(f"./particles/signal_match_{OUT_RUN}.png") # saving to a file

print("DONE")

""" closing input and output root files """
myfile.Close()
out_file.Close()




