import numpy as np

from scipy.optimize import curve_fit
from sys import argv
from ROOT import TFile

STRIP_USED = 47 if len(argv) < 3 else int(argv[2])
f = TFile(f'./singles/single_Ex9Be_from9Be_detB_E=[{STRIP_USED},{STRIP_USED}],dE=[145,160]_bestPairs_run18-30.root')
h = f.Get('Ex;1')


n_bins = h.GetNbinsX()

# dodajemo -1 na početak liste jer je bin numbering 1-based a ne 0-based
h_data = np.array([-1] + [ h.GetBinContent(i) for i in range(n_bins) ])

zero_point = h.GetXaxis().FindBin(0)  # redni broj bina na energiji od 0 MeV
step_size = 0.05  # MeV
state1_nrg = 2.43  # MeV


""" FUNCTIONS """

""" Cuts out all elements of a that are not around p """
def cut(a, p):
    if p not in a:
        return None
    
    index = np.where(a == p)[0][0]
    i_low, i_high = index, index
    
    for i in range(index-1, -1, -1):
        if a[i] == a[i_low] - 1:
            i_low = i
        else:
            break
            
    for i in range(index+1, len(a)):
        if a[i] == a[i_high] + 1:
            i_high = i
        else:
            break
            
    return a[i_low:i_high+1]


""" Gaussian function for fitting """
def gauss(x, *pars):
    A, mu, sigma = pars
    return A * np.exp(-(x-mu)**2/(2.*sigma**2))


""" Turns bin index to energy """
def nrg(x):
    return (x - zero_point) * step_size


""" Gets mean of fit around a peak """
def get_fit(peak_height, peak_index, peak_cutoff):
    peak_nbh = np.where(h_data/peak_height >= peak_cutoff)[0]
    peak_nbh = cut(peak_nbh, peak_index)
    
    peak_nrg = nrg(peak_nbh)
    peak_count = h_data[peak_nbh]
    pars0 = [ peak_height, nrg(peak_index), step_size ]
    
    return curve_fit(gauss, peak_nrg, peak_count, p0=pars0)[0][1]  # mean distribucije


peak0_height = max(h_data)
peak0_index = np.where(h_data == peak0_height)[0][0]
x0 = get_fit(peak0_height, peak0_index, 0.1)

m = h.GetXaxis().FindBin(state1_nrg)
peak1_height = max(h_data[m-10:m+10])
peak1_index = np.where(h_data == peak1_height)[0][0]
x1 = get_fit(peak1_height, peak1_index, 0.1)

a = state1_nrg / (x1 - x0)
b = -a*x0

print(f'{a}\t{b}')

with open('strip_9be_correction', 'a') as f:
    f.write(f'{STRIP_USED}\t{a}\t{b}\n')
