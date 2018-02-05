# Plot the (FE) efficiencies

import numpy as np
from matplotlib import pyplot as plt 
import pandas as pd

datafile = "Peaks.dat"
#Eg[MeV?] 	 nEvents 	 nCounts 	 cntFE 	 cntSE 	 cntDE 	 cnt511 	 cntCompt=Rest
data = np.loadtxt(datafile)
data =  data[np.lexsort(np.fliplr(data).T)] # sort by gamma energy
# print data[:,0]

# df = pd.read_csv(datafile, index_col=None,header=0, sep= "\t")
# df=df.rename(columns={'# Eg[MeV?]':'Eg', 'cntCompt=Rest':'rest'})
# df=df.sort_values(by="Eg")
# # print list(df)
# data = df.as_matrix(columns=None)
# print data


FE_eff = data[:,3] / data[:,1]
SE = (data[:,4]) / data[:,1]
DE = (data[:,5]) / data[:,1]
c511 = (data[:,6]) / data[:,1]
Rest = (data[:,7]) / data[:,1]
Eff_tot = FE_eff + SE + DE + Rest
# nCounts

# print "cntFE peak: ", data[:,3]

fig, ax = plt.subplots()
plt.errorbar(data[:,0], Eff_tot, fmt="o", label="Eff_tot, sim.")
plt.errorbar(data[:,0], FE_eff, fmt="o", label="Eff-FE, sim.")
plt.errorbar(data[:,0], SE, fmt="o", label="Eff_SE, sim.")
plt.errorbar(data[:,0], DE, fmt="o", label="Eff_DE, sim.")
plt.errorbar(data[:,0], Rest, fmt="o", label="Eff_Rest, sim.")


plt.xlabel('Eg [MeV]')
plt.ylabel('eff_FE')

handles, labels = ax.get_legend_handles_labels()
ax.legend(handles, labels, numpoints=1)

plt.show()

