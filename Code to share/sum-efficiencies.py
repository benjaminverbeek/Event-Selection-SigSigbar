# Quick script to sum efficiencies from outfiles.

# Find all files with end .bosslog

import os

dirPath = "/ustcfs/BES3User/2021/benjamin/runarea/test/cmtfile/"

files = os.listdir(dirPath)

files = [f for f in files if f.split(".")[-1] == "bosslog"]

print("Going through all output log-files...")

sumList = [0]*8 # there are 8 filtering steps.

for file in files[:]:

    with open(dirPath + file, 'r') as f:
        lines = f.readlines()[-20:]    # the relevant info is in the last 100 lines for sure.
        lines = [line.strip('\n') for line in lines]

        for i, line in enumerate(lines):
            if "ApplicationMgr       INFO Application Manager Stopped successfully" in line:
                startIndex = i + 1 # skip current line
            elif "ApplicationMgr       INFO Application Manager Finalized successfully" in line:
                endIndex = i - 1 # skip current line
        
        for i, line in enumerate(lines[startIndex:endIndex]):
            value = int(line.split()[-1])
            sumList[i] += value


#print(sumList)
tot, nGood, nGam, pid, lambd, lambdbar, pass4c, passAll = sumList

print("*** TOTAL EFFICIENCIES OF ALL FOUND .bosslog-FILES ***")
print("total number:         ", tot         , f"  {tot/tot*100:.3f} % "     , f"  {tot/tot*100:.3f} % "     )
print("nGood==4, nCharge==0: ", nGood       , f"  {nGood/tot*100:.3f} % "   , f"  {nGood/tot*100:.3f} % "   )
print("nGam>=2:              ", nGam        , f"  {nGam/tot*100:.3f} % "    , f"  {nGam/nGood*100:.3f} % "    )
print("Pass Pid:             ", pid         , f"  {pid/tot*100:.3f} % "     , f"  {pid/nGam*100:.3f} % "     )
print("Pass lambda tag:      ", lambd       , f"  {lambd/tot*100:.3f} % "   , f"  {lambd/pid*100:.3f} % "   )
print("Pass lambdabar tag:   ", lambdbar    , f"  {lambdbar/tot*100:.3f} % ", f"  {lambdbar/lambd*100:.3f} % ")
print("Pass 4C:              ", pass4c      , f"  {pass4c/tot*100:.3f} % "  , f"  {pass4c/lambdbar*100:.3f} % "  )
print("J/psi->Sig0 Sigbar0:  ", passAll     , f"  {passAll/tot*100:.3f} % " , f"  {passAll/pass4c*100:.3f} % " )
