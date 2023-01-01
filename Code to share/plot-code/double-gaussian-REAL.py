# program to fit a double gaussian to some 2D-data
# now instead of using mock-data, using binned data from root-file.

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import uproot

# data filename
baseFilePath = r"C:\Users\benja\OneDrive\Documents\15c project\USTC\all_ana_roots"
# filename = "SigmaSigmabar_out_99k"    # This is 99k exclusive MC-events.
# filename = "inclusiveMC_raw"          # This is 225 M inclusive MC-events, no filtering
filename = "REAL-DATA_ana_filtered"   # This is 225 M REAL DATA filtered for Sigma-masses

# make your best guess:
initialGuess = [500, 1.19, 0.001, 100, 1.19, 0.005]
nBins = 80
mSigma0  = 1.192642 # Sigma0 mass, PDG 2022 (GeV)
mSigma0_err = 0.024 # standard deviation, PDG 2022 (GeV)
n_err = 3   # width of sigma-window for mass-filtering

# define the double gaussian function
def double_gaussian(x, A1, mu1, sigma1, A2, mu2, sigma2):
    return A1*np.exp(-(x-mu1)**2/(2*sigma1**2)) + A2*np.exp(-(x-mu2)**2/(2*sigma2**2))

# read data from root-file
file = uproot.open(baseFilePath + "\\" + filename + ".root")
# navigate to tree named "fit4c"
tree = file["fit4c"]
# navigate to branch named "m_mSigma0"
branch = tree["m_m_Sigma0"]
# get data
data = branch.array()
data = np.array(data)
print("Number of events: ", len(data))
data = data[data < 1.28]
data = data[data > 1.15]
# bin the data and store in x and y
y, x = np.histogram(data, bins=nBins, density=False)
# shift x to the middle of the bins
x = (x[1:] + x[:-1])/2
print("Maximum counts: ", max(y))
print(sum(y[x > 1.22])) 

# linspace for fitted functions
xspace = np.linspace(min(x), max(x), 1000)

# suppress runtime warnings
np.seterr(all='ignore')

# fit the data, but constrain sigma1, sigma2, A1 and A2 to be positive
popt, pcov = curve_fit(double_gaussian, x, y, bounds=([0, -np.inf, 0, 0, -np.inf, 0], np.inf), p0=initialGuess)
# store fitted parameters and their uncertainties
fit_A1, fit_mu1, fit_sigma1, fit_A2, fit_mu2, fit_sigma2 = popt
sig_A1, sig_mu1, sig_sigma1, sig_A2, sig_mu2, sig_sigma2 = np.sqrt(np.diag(pcov))

# plot the data and the fit
plt.hist(data, bins=nBins, density=False, label='data (hist)', alpha=0.3, edgecolor='black', linewidth=1)
plt.plot(x, y, 'b+', label='data')
#plt.show()
plt.plot(xspace, double_gaussian(xspace, *popt), 'r-', label='double-gaussian fit', alpha=0.5)


# print the fit parameters and their uncertainties with 4 decimals. Also print their deviation from the test parameters
print(f'A1 = {fit_A1:.4f} +/- {sig_A1:.4f} ')
print(f'mu1 = {fit_mu1:.4f} +/- {sig_mu1:.4f} ')
print(f'sigma1 = {fit_sigma1:.4f} +/- {sig_sigma1:.4f} ')
print(f'A2 = {fit_A2:.4f} +/- {sig_A2:.4f} ')
print(f'mu2 = {fit_mu2:.4f} +/- {sig_mu2:.4f} ')
print(f'sigma2 = {fit_sigma2:.4f} +/- {sig_sigma2:.4f} ')

# also write this info in the plot in a textbox
plt.figtext(0.65, 0.6, f'A1 = {fit_A1:.4f} +/- {sig_A1:.4f} ', fontsize=8)
plt.figtext(0.65, 0.55, f'mu1 = {fit_mu1:.4f} +/- {sig_mu1:.4f} ', fontsize=8)
plt.figtext(0.65, 0.5, f'sigma1 = {fit_sigma1:.4f} +/- {sig_sigma1:.4f} ', fontsize=8)
plt.figtext(0.65, 0.45, f'A2 = {fit_A2:.4f} +/- {sig_A2:.4f} ', fontsize=8)
plt.figtext(0.65, 0.4, f'mu2 = {fit_mu2:.4f} +/- {sig_mu2:.4f} ', fontsize=8)
plt.figtext(0.65, 0.35, f'sigma2 = {fit_sigma2:.4f} +/- {sig_sigma2:.4f} ', fontsize=8)


# now also plot the two gaussians separately, with their area colored
plt.plot(xspace, double_gaussian(xspace, popt[0], popt[1], popt[2], 0, 0, 0), 'g--', label='gaussian 1', alpha=0.5)
plt.plot(xspace, double_gaussian(xspace, 0, 0, 0, popt[3], popt[4], popt[5]), 'y--', label='gaussian 2', alpha=0.5)
plt.fill_between(xspace, double_gaussian(xspace, popt[0], popt[1], popt[2], 0, 0, 0), color='g', alpha=0.3)
plt.fill_between(xspace, double_gaussian(xspace, 0, 0, 0, popt[3], popt[4], popt[5]), color='y', alpha=0.3)

plt.legend()
plt.xlabel('m(mSigma0) [GeV/c$^2$]')
plt.ylabel('counts')
plt.title('Double Gaussian fit to m(mSigma0) real data')
#plt.show()
plt.savefig('C:/Users/benja/OneDrive/Documents/15c project/USTC/plot-code/imgs/double-gaussian-REAL-sigma.png')
print(f"Saved plot to C:/Users/benja/OneDrive/Documents/15c project/USTC/plot-code/imgs/double-gaussian-REAL-sigma.png")
