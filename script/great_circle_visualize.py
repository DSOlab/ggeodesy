#! /usr/bin/python

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

csv_file = "great_circle.csv"

##  read in great-circle reults, in csv format
df = pd.read_csv(csv_file)

##  plot differences for each column with ref to Vincenty
"""
plt.ylim([-10,10])
plt.ylabel('Algorithm Diff. from Vincenty (m)')
plt.xlabel('Average Distance Between Two Points (km)')
for col in range(1, df.shape[1]):
    cSeriesObj = df.iloc[:, col]
    absDiffs = np.absolute(cSeriesObj.to_numpy() - df['Vincenty'].to_numpy())
    plt.scatter(df['Av.Distance(km)'], absDiffs, alpha=0.5)
    plt.plot(df['Av.Distance(km)'], absDiffs)
plt.legend(tuple(list(df.columns[1:])), loc='upper left')
plt.show()
"""

fig, axis = plt.subplots(df.shape[1]-1, 1, sharex=True)
##  Remove horizontal space between axes
#fig.subplots_adjust(hspace=0)
plt.ylabel('Algorithm Diff. from Vincenty (m)')
# plt.xlabel('Average Distance Between Two Points (km)')
axisn = 0
for col in range(1, df.shape[1]):
    cSeriesObj = df.iloc[:, col]
    absDiffs = np.absolute(cSeriesObj.to_numpy() - df['Vincenty'].to_numpy())
    if max(absDiffs)>1e4: axis[axisn].set_yscale('log')
    axis[axisn].scatter(df['Av.Distance(km)'], absDiffs, alpha=0.5)
    axis[axisn].plot(df['Av.Distance(km)'], absDiffs)
    for lnat in [ 10**x for x in range(0,5) if 10**x<max(absDiffs)]: axis[axisn].axhline(lnat)
    axis[axisn].set_ylabel(df.columns[axisn+1])
    axisn += 1
plt.show()
