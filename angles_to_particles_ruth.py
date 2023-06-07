# -*- coding: utf-8 -*-
"""
Created on Tue Oct 29 06:21:21 2019

@author: Deni

This is a script that looks at the positions of detector centers and from it calculates the positions of all the other detector parts (pixels).
This enables us to assign the geometrical information to the particles detected at these pixels. 
    
"""


""" import predefined functions """
import math as m

from datetime import timedelta
from numpy import zeros, deg2rad, rad2deg, cos, sin
from sys import argv
from time import time
from ROOT import TFile, TTree


""" 
Define program-wide info
    !ADJUST!: N_RUN, CALIB, PARTICLES, DETECTORS_INFO
        HAS_dE (True,False), HAS_ENERGY (True,False), PARTICLE_ENTRY(True,False), EVENT_ENTRY(True,False), INTERSTRIP_MARKING (True,False), CONFLICT_MARKING (True,False) 
        DEBUG (True,False)
"""

CALIB = 'RUTH' # the name of the calibration
PARTICLES = 'particles(E=E,F1-4)'#info about the particle finding process
S_RUN = 'run26' if len(argv) < 2 else argv[1]  # S_RUN is the first argument of the scripts, 'run26' if not given
S_RUN = f'{S_RUN}_NRG({CALIB})_{PARTICLES}'
TREE_NAME = 'tree' #usually we used names "tree" or "T"

DEBUG = False # do we want an output check of certain parts of the code

HAS_ENERGY = True # do we have already the energies of the hits, or do we only have amplitudes

    
#Define the output tree entry structure
PARTICLE_ENTRY = False  # will an entry in the outgoing file be a single particle or an entire event?

#maximum number of particles that we write into a single root entry
if PARTICLE_ENTRY: 
    MAX_PARTICLES_ENTRY = 1 
else:
    MAX_PARTICLES_ENTRY = 30 # any number that is higher than the expected maximum number of particles in an experimental event, 10 is enough when the conflict filter is on, more if it is not

CONFLICT_MARKING = True # do we have conflicted particles who share the same adc strip marked
INTERSTRIP_MARKING = True  # do we have potentital interstrip particles marked

HAS_dE = True # is there a deltaE detector in the run

if HAS_dE:
    MIN_mult = 3 # minimum number of hits in the event needed to define a physical event, dE, FRONT and BACK info needed to constitue a single particle
else: 
    MIN_mult = 2 # minimum number of hits in the event needed to define a physical event, FRONT and BACK info needed to constitue a single particle

#Information about the detectors: 
# angles and radial distance of the detector center from the target
# number of stripc (adc channels) per detector and the detector width
DETECTORS_INFO={'theta_center':[52.125,20,24.9,52.975], 'phi_center': [0,0,180,180], 'distance_center':[0.147,0.312,0.299,0.147], 'num_strips':16, 'detector_width':0.05 , }    

# the name of the output file 
OUT_RUN = S_RUN.replace('NRG(', '').replace(')_', '_')

 
    
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
Functions for transformation of spherical coordinates r, theta, phi into cartesian coordinates x,y,z and vice a versa
"""
def spher_2_cart(r, theta, phi):
    x = r * sin(theta) * cos(phi)
    y = r * sin(theta) * sin(phi)
    z = r * cos(theta)
  
    return (x, y, z)    

def cart_2_spher(x,y,z):
    r = m.sqrt(x**2 + y**2 + z**2)
    theta = m.atan(m.sqrt(x**2 + y**2) / z)
    phi = m.atan2(y, x)    

    return (r, theta, phi)
    
"""
A function to calculate the cartesian coordinates of the detector centers, transfer angles of the detector centers from degrees to radians and append that info to "detectors_info"
"""

def calculate_detector_centers(detectors_info):
    if detectors_info is not None:
        num_detectors = len(detectors_info['theta_center'])
        detector_center_cart = []
        theta_center_rad = zeros(num_detectors)
        phi_center_rad = zeros(num_detectors) 
        r_center = zeros(num_detectors)            
        
        for i in range(num_detectors):
            theta_center_rad[i] = deg2rad(detectors_info['theta_center'][i])
            phi_center_rad[i] = deg2rad(detectors_info['phi_center'][i])
            r_center[i] = detectors_info['distance_center'][i]
            
            detector_center_cart.append(spher_2_cart(r_center[i], theta_center_rad[i], phi_center_rad[i]))
            
        detectors_info['detector_center_cart'] = detector_center_cart
        detectors_info['theta_center_rad'] = theta_center_rad
        detectors_info['phi_center_rad'] = phi_center_rad
        
        if DEBUG:
            print(detectors_info)


"""
Functions to calculate angles of all pixels:
    pixel_offset: calculates pixel offset (a,b) from the center of the detector
    get_one_pixel_coordinates: calculates spherical coordinates for a pixel
    get_all_pixel_coordinates: calculates and records spherical coordinates of all pixels

Details of the calculation:    
    pixels are enumerated from 0 to n-1
    bottom right pixel is (0,0)
    front strips(adcs) constitute columns
    back strips(adcs) constitute rows
    one dimension of the detector coincides with the "y" axis, while the other is in the "x-z" plane
    "a" and "b" change signs depending on the part of the detector where the pixel is situated
    cos and sin of "thet_center" are always positive, while cos(phi_center) changes sign while moving from one side of the beam (z axis) to the other
"""
       
def pixel_offset(row, column, num_strips, pixel_width):
    center = num_strips / 2.    
    a = (center - column - 0.5) * pixel_width #in "x-z" plane
    b = (center - row - 0.5) * pixel_width #along the "y" axis
    
    return (a, b)
       
def get_one_pixel_coordinates(row, column, num_strips, pixel_width, detector_center_cart, theta_center_rad, phi_center_rad):
    a, b = pixel_offset(row, column, num_strips, pixel_width)
    x = detector_center_cart[0] + a * m.cos(theta_center_rad) 
    y = detector_center_cart[1] + b
    z = detector_center_cart[2] - a * m.sin(theta_center_rad) * m.cos(phi_center_rad)
    r, theta, phi = cart_2_spher(x, y, z)

    cartesian = [x, y, z]
    spherical = [r, theta, phi]
    return (cartesian, spherical)
        
def get_all_pixel_coordinates(detectors_info):
    num_detectors = len(detectors_info['theta_center'])
    num_strips = detectors_info['num_strips']
    pixel_width = detectors_info['detector_width']/num_strips
    
    pixel_coordinates = []
    for i in range(num_detectors):
        pixel_coordinates.append([])
        theta_center_rad = detectors_info['theta_center_rad'][i]
        phi_center_rad = detectors_info['phi_center_rad'][i]
        detector_center_cart = detectors_info['detector_center_cart'][i]
        for row in range(num_strips):
            pixel_coordinates[i].append([])
            for column in range(num_strips):
                cartesian, spherical = get_one_pixel_coordinates(row, column, num_strips, pixel_width, detector_center_cart, theta_center_rad, phi_center_rad)
                pixel_coordinates[i][row].append([spherical[0], rad2deg(spherical[1]), rad2deg(spherical[2])])
    
    return pixel_coordinates
    


""" EXCECUTION """ 
       
"""open the file and get the tree"""
myfile = TFile(f'./particles/{S_RUN}.root') #open
if DEBUG:
    myfile.ls() #check the file contents

t = myfile.Get(TREE_NAME)  #get the tree
if DEBUG:
    t.Print() #check the tree contents
            

"""creating the output root file; defining name, output tree, output tree variables, output tree structure """
out_file = TFile(f'./angles/{OUT_RUN}.root', 'recreate') # output root file
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

#new variables           
r = zeros(MAX_PARTICLES_ENTRY, dtype='float') # distance from the target to the pixel where the particle is detected, field of values of type Float
theta = zeros(MAX_PARTICLES_ENTRY, dtype='float') # angle from the target to the pixel where the particle is detected, field of values of type Float
phi = zeros(MAX_PARTICLES_ENTRY, dtype='float') # angle from the target to the pixel where the particle is detected, field of values of type Float

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
    out_t.Branch('r', r, 'r/D') # particle distance branch
    out_t.Branch('theta', theta, 'theta/D') # particle angle branch
    out_t.Branch('phi', phi, 'phi/D') # particle angle branch
else:
    out_t.Branch('adc',adc,'adc[mult]/S') # channel number branch, field of lenght "mult"
    out_t.Branch('ampl',ampl,'ampl[mult]/S') # amplitude branch, field of lenght "mult"
    if  HAS_ENERGY:
        out_t.Branch('nrg', nrg, 'nrg[mult]/D') # energy branch, field of lenght "mult"
    out_t.Branch('detector', detector, 'detector[cnc]/S') # detector branch, field of lenght "cnc"  
    out_t.Branch('r', r, 'r[cnc]/D') # particle distance branch
    out_t.Branch('theta', theta, 'theta[cnc]/D') # particle angle branch
    out_t.Branch('phi', phi, 'phi[cnc]/D') # particle angle branch

if CONFLICT_MARKING:
    out_t.Branch('conflicts', conflicts, 'conflicts/S') #conflict branch
if INTERSTRIP_MARKING:
    out_t.Branch('interstrip', interstrip, 'interstrip/S') #interstrip branch    



""" go trough the run data """
n_entries = t.GetEntriesFast() # get the  number of entries in the tree
if DEBUG:
    print('number of entries: ', n_entries)
    n_entries = 30
    
start_time = time() #start the timer 


#calculating and defining various needed variables
calculate_detector_centers(DETECTORS_INFO) #cartesian coordinates of detector centres
pixel_coordinates = get_all_pixel_coordinates(DETECTORS_INFO) #all pixel coordinates
num_strips = DETECTORS_INFO['num_strips'] #number of strips(adcs) on a detector
row = zeros(MAX_PARTICLES_ENTRY, dtype='int') # variable for rows which correspond to horizontal (back) strips
column = zeros(MAX_PARTICLES_ENTRY, dtype='int') # variable for columns which correspond to vertical (front) strips


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
    
    # variables that are different for different particles in an event: detector, adc, ampl, nrg, r, theta, phi    
    if PARTICLE_ENTRY:
        detector[0]=t.detector
        
        #column(front strip) and row(back strip) of the detected particle, set to have values [0,15] by taking the remainder of the division by 16
        column[0] = (t.adc[0] - 1) % num_strips
        row[0] = (t.adc[1] - 1) % num_strips # (16 - t.adc[1]) % num_strips,  if we want to invert the ordering of back strips

        #distance and angles of the pixel where the particle was detected
        r[0] = pixel_coordinates[detector-1][row[0]][column[0]][0]
        theta[0] = pixel_coordinates[detector-1][row[0]][column[0]][1]
        phi[0] = pixel_coordinates[detector-1][row[0]][column[0]][2]
        #copying adc, ampl and energy info from the previous run          
        for i in range(3):
            adc[i] = t.adc[i]
            ampl[i] = t.ampl[i]
            if HAS_ENERGY:
                nrg[i] = t.nrg[i]
                
        out_t.Fill() # filling the out tree; particle is the output entry
    else:
        if cnc[0] > MAX_PARTICLES_ENTRY:
            continue

        #going through all the particles in an event
        for i in range(cnc[0]):
            front_counter = 3*i #counter that goes though field positions that contain information about the front strip
            back_counter = 3*i + 1 #counter that goes through field positions that contain information about the back strip
            
            detector[i] = t.detector[i]
            
            #column(front strip) and row(back strip) of the detected particle, set to have values [0,15] by taking the remainder of the division by 16
            column[i] = (t.adc[front_counter] - 1) % num_strips
            row[i] = (t.adc[back_counter] - 1) % num_strips # (16 - t.adc[back_counter] % num_strips),  if we want to invert the ordering of back strips
            
            #distance and angles of the pixel where the particle was detected
            r[i] = pixel_coordinates[detector[i]-1][row[i]][column[i]][0]
            theta[i] = pixel_coordinates[detector[i]-1][row[i]][column[i]][1]
            phi[i] = pixel_coordinates[detector[i]-1][row[i]][column[i]][2]
        
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
