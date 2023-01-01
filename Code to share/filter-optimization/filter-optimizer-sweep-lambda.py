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
SIGNAL_ID = 2   # The ID of the signal reaction, retrieved from topoana-execution.
sigma_window = 0.01 # +/- this value from estimated value
lambda_window = 0.0083
chi2_cut = 60 # cut on chi2, less than this
weight = 10 # weight of background events, 1 = same weight as signal events

cut_max = 0.04 # (already cut 200 in ES.cxx)
cut_min = 0.001 
# get a linspace of chi2 values to sweep over
sweep_range = np.linspace(cut_min, cut_max, 50)

mSigma_est = 1.1918 # estimate of sigma0 mass from double gaussian fit
mLambda_est= 1.1159 # estimate of lambda mass from double gaussian fit



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


FOMs = []   # list of FOMs
for lambda_window in sweep_range:
    # Possibly add requirements for the other items too, e.g. mSigma0-bar.
    filterList_sigma = (mSig < mSigma_est + sigma_window) & (mSig > mSigma_est - sigma_window) & (mSigbar < mSigma_est + sigma_window) & (mSigbar > mSigma_est - sigma_window)
    filterList_lambda= (mLambda < mLambda_est + lambda_window) & (mLambda > mLambda_est - lambda_window) & (mLambdabar < mLambda_est + lambda_window) & (mLambdabar > mLambda_est - lambda_window)
    filterList_chi2  = (chi2 < chi2_cut)
    filterList = (filterList_sigma & filterList_lambda & filterList_chi2)

    # ^ list of True/False for each event, True if event passes filter.

    # Now we want to filter the data according to the above criteria. We then store B and S.

    iDcyTr_filtered = iDcyTr[filterList]
    B = sum(iDcyTr_filtered != SIGNAL_ID)
    S = sum(iDcyTr_filtered == SIGNAL_ID)

    FOM = S/np.sqrt(S+weight*B) # check +/- 1 sigma of FOM from max, pick lowest - gives high suppression while maintaining good statistical uncertainty.
    # FOM = S/(S+B)
    # relative_BG = B/(S+B)*100 # in percent
    # efficiency = S/23032*100 # in percent, relative to 23032 signal events, no filtering after ES
    # FOM = 100-relative_BG if efficiency > 70 else 0 # if BG is too high, set FOM to 0.


    dFdS = (2*B+S)/(2*(S+B)**(3/2))

    FOMs.append(FOM)

    #print(f"Total events: {len(iDcyTr)}")

#     print(f"""{"-"*30}
#     Background: {B} 
#         Signal: {S} 
#         FOM: {FOM:.6f} +/- {dFdS:.6f}
# Window size: {lambda_window:.6f} GeV
#     {"-"*30}""")
#     print(f"Efficiency: {S/23032*100:.6f} % \t (relative to 23032 signal events, no filtering after ES)")


file.close()    # not needed anymore.
print("DONE.")

# Now plot the chi2_sweep vs FOMs
plt.plot(sweep_range, FOMs, label="FOM")
plt.xlabel("lambda window [GeV]", fontsize = 12)
plt.ylabel("FOM", fontsize = 12)
plt.grid()
# find intersection of this line with FOMs
# find the first index of FOMs that is greater than 0.98*max_FOMFOM
index98 = next((index for index, value in enumerate(FOMs) if value > 0.98*max(FOMs)), None)
# Get the FOM value at this index
FOM98 = FOMs[index98]
# mark this point on the plot
plt.scatter(sweep_range[index98], FOM98, marker="+", color="green", label=f"FOM>0.98 max : {FOM98:.6f} at {sweep_range[index98]:.6f} GeV")

# plot a horizontal red dashed line at max FOM - 5% max
plt.axhline(0.98*max(FOMs), color="red", linestyle="--", label=f"98% of max FOM: {max(FOMs)*0.98:.6f}", linewidth=0.8)

# and mark the highest FOM
max_FOM = max(FOMs)
max_FOM_index = FOMs.index(max_FOM)
plt.scatter(sweep_range[max_FOM_index], max_FOM, marker="x", color="red", label=f"Max FOM: {max_FOM:.6f} at {sweep_range[max_FOM_index]:.6f} GeV")
plt.legend()
print(f"Max FOM: {max_FOM:.6f} at {sweep_range[max_FOM_index]:.6f} GeV")
print(f"First FOM > 0.98 max : {FOM98:.6f} at {sweep_range[index98]:.6f} GeV")

plt.title("FOM vs $w_{\Lambda}$", fontsize=12)
plt.figtext(0.5, 0.35, f'$w_\Sigma^0$ = {sigma_window:.4f} GeV, $c_\chi$ = {chi2_cut:.4f}', fontsize=10)
plt.savefig(r"C:\Users\benja\OneDrive\Documents\15c project\USTC\all_ana_roots"+f"/filter-optimization/plots/FOM_vs_wLambda", dpi=300)
plt.show()
print("Saved figure.")


# Todo: sweep over the parameters, and store the FOM to a list.
# Then plot the swept parameter vs FOM. Do one parameter at a time for now (though I could use scipy?)