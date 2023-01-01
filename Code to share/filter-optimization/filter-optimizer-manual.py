# This program aims to optimize the filtering of the inclusive MC data.
#
# It will return the number of signal and background events, as well as
# our figure of merit: S / sqrt(S + B)
# 
# It first performs a filtering according to some criteria, and then 
# stores the figure of merit and number of events to a list.
#
# BASELINE: this filtering gave sub 1% BG with approx 15% efficiency.
# m_mSigma0>=1.180642 && m_mSigma0<=1.214642 && m_mSigmabar0<=1.214642 && m_mSigmabar0>=1.180642 && mLambda<=1.125683 && mLambda>=1.105683 && mLambdabar<=1.125683 && mLambdabar>=1.105683 && chi2<=60
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

# NOTE: Doesn't work. Root-output is bad format. Using .card instead.
import numpy as np
import matplotlib.pyplot as plt
import uproot

SIGNAL_ID = 2   # The ID of the signal reaction, retrieved from topoana-execution.
sigma_window = 0.01 # +/- this value from estimated value
lambda_window = 0.0082
chi2_cut = 60 # cut on chi2, less than this


mSigma_est = 1.1918 # estimate of sigma0 mass from double gaussian fit
mLambda_est= 1.1159 # estimate of lambda mass from double gaussian fit

# We want to maximize the FOM.

# # REFERENCE VALUES
# mSigma0_pdg  = 1.192642 
# mSigma0_err_pdg = 0.024/1000 # standard deviation, PDG 2022 (GeV)
# mLambda_pdg = 1.115683 
# mLambda_err_pdg = 0.006/1000 # standard deviation, PDG 2022 (GeV)

# n_err = 3   # width of sigma-window for mass-filtering

# sigma_est_err = 0.0033 # standard deviation of narrow fit, approx. (GeV)
# lambda_est_err = 0.001 # standard deviation of narrow fit, approx. (GeV)
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

print("START.")
# data filename
baseFilePath = r"C:\Users\benja\OneDrive\Documents\15c project\USTC\all_ana_roots"
filename = "inclusiveMC_raw"
# read data from root-file
file = uproot.open(baseFilePath + "\\" + filename + ".root")

# navigate to tree named "fit4c"
tree = file["fit4c"]
# navigate to branch named "m_mSigma0"
branch_mSig = tree["m_mSigma0"]
branch_mSigbar = tree["m_mSigmabar0"]
branch_mLambda = tree["mLambda"]
branch_mLambdabar = tree["mLambdabar"]
branch_chi2 = tree["chi2"]
branch_iDcyTr = tree["iDcyTr"]  # Stores indices of decay-type.

# get data for filtering
mSig = branch_mSig.array()
mSigbar = branch_mSigbar.array()
mLambda = branch_mLambda.array()
mLambdabar = branch_mLambdabar.array()
chi2 = branch_chi2.array()
iDcyTr = branch_iDcyTr.array()

# Possibly add requirements for the other items too, e.g. mSigma0-bar.
filterList_sigma = (mSig < mSigma_est + sigma_window) & (mSig > mSigma_est - sigma_window) & (mSigbar < mSigma_est + sigma_window) & (mSigbar > mSigma_est - sigma_window)
filterList_lambda= (mLambda < mLambda_est + lambda_window) & (mLambda > mLambda_est - lambda_window) & (mLambdabar < mLambda_est + lambda_window) & (mLambdabar > mLambda_est - lambda_window)
filterList_chi2  = (chi2 < chi2_cut)
filterList = (filterList_sigma & filterList_lambda & filterList_chi2)

# ^ list of True/False for each event, True if event passes filter.

# Now we want to filter the data according to the above criteria. We then store B and S.
file.close()    # not needed anymore.

iDcyTr_filtered = iDcyTr[filterList]
B = sum(iDcyTr_filtered != SIGNAL_ID)
S = sum(iDcyTr_filtered == SIGNAL_ID)

FOM = S/np.sqrt(S+B)

print(f"Total events: {len(iDcyTr)}")

print(f"""{"-"*30}
Background: {B} 
    Signal: {S} 
       FOM: {FOM:.6f}
Relative BG: {B/(B+S)*100:.6f} %
Using: sigma_window = {sigma_window}, lambda_window = {lambda_window}, chi2_cut = {chi2_cut}  
{"-"*30}""")
print(f"Efficiency: {S/23032*100:.6f} % \t (relative to 23032 signal events, no filtering after ES)")



print("DONE.")


# Todo: sweep over the parameters, and store the FOM to a list.
# Then plot the swept parameter vs FOM. Do one parameter at a time for now (though I could use scipy?)