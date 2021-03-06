#!/usr/bin/env python

import os
import re
import sys
import time

import numpy as np
import matplotlib.pyplot as plt

if len(sys.argv) < 2:
    exit("Usage: ./plotModels [type]")

typ = str(sys.argv[1])

fig = plt.figure(1)
ax = fig.add_subplot(111)
ax.set_xlim([0,243])
    
true = typ + "_true.txt"
    
plt.ion()
plt.show()
ax.set_autoscale_on(False)
for file in os.listdir('./'):

    if not file.startswith(typ + "_0"):
        continue

    if typ == "mu":
        ax.set_ylim([10, 20])
    elif typ == "theta":
        ax.set_ylim([0.15, 0.45])
    else:
        exit("UNKNOWN TYPE")
    
    itr = re.findall(r'\d+', file)
    plt.title('Iteration: ' + itr[0])
    values = np.loadtxt(file)
    trueVal = np.loadtxt(true)
    plt.plot(values)
    plt.plot(trueVal)
    plt.draw()
    plt.cla()    