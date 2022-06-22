import numpy as np

from scipy.optimize import curve_fit
from sys import argv
from ROOT import TFile

STRIP_USED = 47 if len(argv) < 2 else int(argv[1])
STEP_SIZE = 0.05  # MeV
SEARCH_WIDTH = 30  # * step_size = search_width in MeV

f = TFile(f'./singles/calib/6he/single_Ex12C_from6He_E=[{STRIP_USED},{STRIP_USED}],dE=[112,192]_bestPairs_run18-30.root')
h = f.Get('Ex;1')


n_bins = h.GetNbinsX()

# dodajemo -1 na početak liste jer je bin numbering 1-based a ne 0-based
h_data = np.array([-1] + [ h.GetBinContent(i) for i in range(n_bins) ])

zero_point = h.GetXaxis().FindBin(0)  # redni broj bina na energiji od 0 MeV
state0_nrg = 0  # MeV
state1_nrg = 4.44  # MeV


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
    return (x - zero_point) * STEP_SIZE


""" Gets mean of fit around a peak """
def get_fit(search_nrg, peak_cutoff=0.1):
    m = h.GetXaxis().FindBin(search_nrg)
    peak_height = max(h_data[m-SEARCH_WIDTH:m+SEARCH_WIDTH])

    tmp = np.where(h_data == peak_height)[0]
    i = 0
    peak_index = tmp[i]
    while peak_index not in range(m-SEARCH_WIDTH, m+SEARCH_WIDTH):
        i += 1
        peak_index = tmp[i]

    print(f'peak around {search_nrg} MeV: E={nrg(peak_index)}, N={peak_height}')

    peak_nbh = np.where(h_data/peak_height >= peak_cutoff)[0]
    peak_nbh = cut(peak_nbh, peak_index)
    
    peak_nrg = nrg(peak_nbh)
    peak_count = h_data[peak_nbh]
    pars0 = [ peak_height, nrg(peak_index), STEP_SIZE ]
    
    mean = curve_fit(gauss, peak_nrg, peak_count, p0=pars0)[0][1]  # mean distribucije
    print(f'\tfit mean at {mean} MeV')

    return mean

E0 = get_fit(state0_nrg, peak_cutoff=0.1)
E1 = get_fit(state1_nrg, peak_cutoff=0.1)

a = (state1_nrg - state0_nrg) / (E1 - E0)
b = state0_nrg - a*E0

print(f'{a}\t{b}')

with open('strip_6he_correction', 'a') as f:
    f.write(f'{STRIP_USED}\t{a}\t{b}\n')
