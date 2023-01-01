# program to fit a double gaussian to some 2D-data
# now instead of using mock-data, using binned data from root-file.


# NOTE: Doesn't work. Root-output is bad format. Using .card instead.
import numpy as np
import matplotlib.pyplot as plt
import uproot
import awkward as ak

sigma_window = 0.01 # +/- this value from estimated value
lambda_window = 0.0082
chi2_cut = 60 # cut on chi2, less than this

mSigma_est = 1.19175 # estimate of sigma0 mass from double gaussian fit
mLambda_est= 1.115883 # estimate of lambda mass from double gaussian fit

# data filename
baseFilePath = r"C:\Users\benja\OneDrive\Documents\15c project\USTC\all_ana_roots"
filename = "inclusiveMC_raw"
# filename = "REAL-DATA_ana.root"
# filename = "filtered"
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

# get data for filtering
mSig = branch_mSig.array()
mSigbar = branch_mSigbar.array()
mLambda = branch_mLambda.array()
mLambdabar = branch_mLambdabar.array()
chi2 = branch_chi2.array()

# Possibly add requirements for the other items too, e.g. mSigma0-bar.
filterList_sigma = (mSig < mSigma_est + sigma_window) & (mSig > mSigma_est - sigma_window) & (mSigbar < mSigma_est + sigma_window) & (mSigbar > mSigma_est - sigma_window)
filterList_lambda= (mLambda < mLambda_est + lambda_window) & (mLambda > mLambda_est - lambda_window) & (mLambdabar < mLambda_est + lambda_window) & (mLambdabar > mLambda_est - lambda_window)
filterList_chi2  = (chi2 < chi2_cut)
filterList = (filterList_sigma & filterList_lambda & filterList_chi2)
print(f"Number of events passing filter: {sum(filterList)}")
# ^ list of True/False for each event, True if event passes filter.

# Print the parameters used for cut
print("sigma_window: " + str(sigma_window))
print("lambda_window: " + str(lambda_window))
print("chi2_cut: " + str(chi2_cut))

dataDict = {}
# Now just store the filtered data to a new .root-file and we're done.
for branch in tree:
    # print(branch.name)
    d = branch.array()
    dataDict[branch.name] = d[filterList]
    # print("ok")
file.close()

# now store this to a new root-file
newFile = uproot.recreate(baseFilePath + "\\" + filename + "_filtered.root")
newFile["fit4c"] = dataDict

# print(dataDict)

print("DONE.")

newFile.close() 
print("Saved: " + baseFilePath + "\\" + filename + "_filtered.root")
