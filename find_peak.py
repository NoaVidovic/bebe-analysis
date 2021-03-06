import numpy as np

from scipy.optimize import curve_fit
from sys import argv
from ROOT import TFile

try:
    IN_NAME = argv[1]
    SEARCH_NRG = float(argv[2])  # MeV

    STEP_SIZE = 0.05  # MeV
    PEAK_CUTOFF = -1 if len(argv) < 4 else float(argv[3])
    SEARCH_WIDTH = 0.2 if len(argv) < 5 else float(argv[4])  # in MeV
    SEARCH_WIDTH = int(SEARCH_WIDTH / STEP_SIZE)  # in steps
except:
    print('arg1: file name')
    print('arg2: search energy')
    print('arg4: search width (optional, default=0.5 [MeV])')
    print('arg3: peak cutoff ratio (optional, default is automatic search)')

f = TFile(IN_NAME)
h = f.Get("Ex';1")  # if using _corr_ files


n_bins = h.GetNbinsX()

# dodajemo -1 na početak liste jer je bin numbering 1-based a ne 0-based
h_data = np.array([-1] + [ h.GetBinContent(i) for i in range(n_bins) ])
zero_point = h.GetXaxis().FindBin(0)  # redni broj bina na energiji od 0 MeV


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

    print(f'peak around {search_nrg:.3f} MeV: E={nrg(peak_index):.3f}, N={peak_height:.0f}')

    peak_nbh = np.where(h_data/peak_height >= peak_cutoff)[0]
    peak_nbh = cut(peak_nbh, peak_index)
    
    peak_nrg = nrg(peak_nbh)
    peak_count = h_data[peak_nbh]
    pars0 = [ peak_height, nrg(peak_index), STEP_SIZE ]
    
    A, mean, sigma = curve_fit(gauss, peak_nrg, peak_count, p0=pars0)[0]  # mean distribucije
    print(f'\tfit ampl {A:.3f} MeV')
    print(f'\tfit mean at {mean:.3f} MeV')
    print(f'\tfit stdev {sigma:.3f} MeV')

    return (A, mean, sigma)

if PEAK_CUTOFF == -1:
    pc = 0.1
    
    while True:
        try:
            print(f'peak cutoff = {pc:.2f}: ', end='')
            _, x, s = get_fit(SEARCH_NRG, pc)

            if not (SEARCH_NRG - SEARCH_WIDTH*STEP_SIZE <= x <= SEARCH_NRG + SEARCH_WIDTH*STEP_SIZE):
                if abs(pc - 0.9) < 0.01:
                    break

                pc += 0.1
            else:
                if pc > 0.1:
                    pc -= 0.1
                break
        except:
            if abs(pc - 0.9) < 0.01:
                break

            pc += 0.1

    print('Fine tuning…')

    while True:
        try:
            print(f'peak cutoff = {pc:.2f}: ', end='')
            _, x, s = get_fit(SEARCH_NRG, pc)

            if SEARCH_NRG - SEARCH_WIDTH*STEP_SIZE <= x <= SEARCH_NRG + SEARCH_WIDTH*STEP_SIZE:
                break
            else:
                pc += 0.01
                if pc >= 1.:
                    print('Unable to fit')
                    break
        except:
            pc += 0.01
            if pc >= 1.:
                print('Unable to fit')
                break
else:
    _, x, s = get_fit(SEARCH_NRG, PEAK_CUTOFF)

print(f'MEAN: {x}')
print(f'STDEV: {s}')
