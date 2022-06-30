# -*- coding: utf-8 -*-
"""
@author: Deni

This is a script for making histograms of channels, amplitudes and energies of hits. Primarily used for checking the quality of data in several steps of the analysis.
For example checking the noise thresholds for the channels. 
"""


""" import predefined functions """
from ROOT import TCanvas, gPad, TFile, TH1F, TH2S
from time import time
from datetime import timedelta
from sys import argv

""" 
Define program-wide info
    !ADJUST!: S_RUN, RUN_INFO, TREE_NAME, N_ADC, DRAW_HIST (True,False), WRITE_HIST (True,False), NRG(True,False)
"""
S_RUN = argv[1] if len(argv) > 1 else 'run26' # name of the run, string form
TREE_NAME = 'tree' #usually we used names "tree" or "T"
RUN_INFO = '9Be on 9Be' # projectile and target
N_ADC = 192 # number of ADC channels (strips) and matching histograms
DRAW_HIST = False # do we require drawing of histograms and saving in .png form
WRITE_HIST = True # do we require drawing of histograms and saving in .root form
NRG = False #do we want to use (do we have) info on energy (True) or amplitude (False)

S_PREFIX = 'HIST(adc,ampl)' # prefix to add to the output files 
if NRG:
    S_PREFIX = 'HIST(adc,nrg)' # prefix to add to the output files 

"""
A function to follow the status of program execution, gives: elapsed time, estimate of remaining time, number of events evaluated up to that point 
"""
def print_time(start_time, step, total):
    elapsed_time = round(time()-start_time)
    percentage = (step*1.)/total
    remaining_time = 0
    if percentage != 0:
        remaining_time = round(elapsed_time * (1-percentage)/percentage)
        remaining_time = timedelta(0,remaining_time)
    elapsed_time = timedelta(0,elapsed_time)
    print("{:.1f}%, elapsed {}, remaining {}, event {} out of {}".format(percentage*100, elapsed_time, remaining_time, step, total))


"""  
definition of histograms: "histos" (ampl/nrg,N), "histo_0" (adc,N), "histo_1" (adc,ampl/nrg)  
    !ADJUST!: put histogram ranges of interest and choose binning
"""
histos = [] #list for N_ADC channel histograms, each records the number of hits with a certain amplitude/energy for that specific channel: N-ampl/energy  (y-x)
for i in range(N_ADC): # for; 0 to N_ADC-1
    if not NRG: #for amplitudes
        name = "adc_{}_(ampl,N)".format(i+1) # define names for histograms
        title = "N-ampl histogram: {}, ADC-{}".format(RUN_INFO,i+1) # define titles for histograms
        histos.append(TH1F(name,title,3500,0,3500)) # create channel histogram and append to the list
    if NRG: #same as above with different range for energies
        name = "adc_{}_(nrg,N)".format(i+1) 
        title = "N-nrg histogram: {}, ADC-{}".format(RUN_INFO,i+1) 
        histos.append(TH1F(name,title,6000,0,60)) 
    
histo_0 = TH1F('(adc,N)', "N-adc histogram: {}".format(RUN_INFO), 192, 0, 192) # records the number of hits in all channels: N-adc  (y-x)
histo_1 = TH2S('(adc,ampl)', "ampl-adc histogram: {}".format(RUN_INFO), 192, 0, 192, 3500, 0, 3500) # 2D histogram records amplitude of all hits in all channels: ampl-adc  (y-x) 
if NRG:  #different range for energies
    histo_2 = TH2S('(adc,nrg)', "nrg-adc histogram: {}".format(RUN_INFO), 192, 0, 192, 6000, 0, 60) # 2D histogram records energy of all hits in all channels: nrg-adc  (y-x)
    
""" 
functions for histogram drawing and saving
    draw_hists: draw histograms and save in .png form: NOT FINISHED!!!
    write_hists: make histograms and save in .root form 
"""
def draw_hists(hists): #NOT FINISHED
    can = TCanvas('hists','Histograms',800,600) #defining canvas
    for i in range(N_ADC):
        hists[i].Draw() # drawing a histogram
        gPad.Update()
        can.SaveAs("{}_adc{}_{}.png".format(S_PREFIX,i+1,S_RUN)) # saving a histogram in a png file
    can1 = TCanvas('hists1','Histograms1',800,600)
    histo_0.Draw()
    histo_1.Draw()
    if NRG:
        histo_2.Draw()
    can1.SaveAs("{}_{}.png".format(S_PREFIX, S_RUN)) # sinmanje u file

def write_hists(hists):
    f = TFile("{}_{}.root".format(S_PREFIX, S_RUN),"recreate")
    
    histo_0.Write()
    histo_1.Write()
    if NRG:
        histo_2.Write()    
    for i in range(N_ADC):
        hists[i].Write()
        
    f.Close()


"""open the file and get the tree"""
myfile = TFile('{}.root'.format(S_RUN)) #open
#myfile.ls() #check the file contents
t = myfile.Get('{}'.format(TREE_NAME))  #get the tree
#t.Print() #check the tree contents
    

""" go trough the run data """
n_entries = t.GetEntriesFast() # get the  number of entries in the tree

start_time = time() #start the timer

for i in range(n_entries): # for: 0 to n_entries-1

    if i % 10000 == 0: # print the timer info every 10000 entries
        print_time(start_time,i,n_entries)
    
    t.GetEntry(i) # get an entry(event)

    for j in range(t.wm): # go trough all the hits in an entry(event)
        
        histo_0.Fill(t.wadc[j])  # fill the histogram with the channel which recorded the hit
        #if t.wampl[j]>2000: 
        histo_1.Fill(t.wadc[j],t.wampl[j])  # fill the 2D histogram with the channel and amplitude of the hit
        
        if not NRG:
            histos[t.wadc[j]-1].Fill(t.wampl[j]) # fill the appropriate channel histogram with amplitude of the hit
            
        if NRG:
            histos[t.wadc[j]-1].Fill(t.wnrg[j]) # fill the appropriate channel histogram with energy of the hit
            #if t.wampl[j]>2000: 
            histo_2.Fill(t.wadc[j],t.wnrg[j])  # fill the 2D histogram with the channel and energy of the hit

#writing the histograms
if WRITE_HIST:
    write_hists(histos)

#drawing the histograms
if DRAW_HIST: 
    draw_hists(histos)


