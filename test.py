from matplotlib import use
use("Agg")
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt

#CONSTANTS
DEG=np.pi/180
RAD=1/DEG

#LOAD DATA
elements=np.loadtxt("elements.dat")

#FIGURE REGION
fig,axs=plt.subplots(5,sharex=True,figsize=(8,12))

#ELEMENTS
ts=elements[:,0]
ees=elements[:,2]
aes=elements[:,1]/(1-ees)
ies=elements[:,3]*RAD
Wes=elements[:,4]*RAD
wes=elements[:,5]*RAD

#OPTIONS
args=dict(marker='None',ls='-',lw=1,color='b');

#PLOTS
i=0
axs[i].plot(ts,aes,**args);i+=1
axs[i].plot(ts,ees,**args);i+=1
axs[i].plot(ts,ies,**args);i+=1
axs[i].plot(ts,Wes,**args);i+=1
axs[i].plot(ts,wes,**args);i+=1

#DECORATION
i=0
axs[i].set_ylabel(r"$a$ ($R_\oplus$)");
axs[i].set_yticks(axs[i].get_yticks()[1:-1])
i+=1
axs[i].set_ylabel(r"$e$");
axs[i].set_yticks(axs[i].get_yticks()[1:-1])
i+=1
axs[i].set_ylabel(r"$i$ (deg)");
axs[i].set_yticks(axs[i].get_yticks()[1:-1])
i+=1
axs[i].set_ylabel(r"$\Omega$ (deg)");
axs[i].set_yticks(axs[i].get_yticks()[1:-1])
i+=1
axs[i].set_ylabel(r"$\omega$ (deg)");
axs[i].set_yticks(axs[i].get_yticks()[1:-1])
i+=1

#GLOBAL
axs[-1].set_xlim((ts[0],ts[-1]))
axs[-1].set_xlabel(r"$t$ (hours)")

#SAVE FIGURE
plt.tight_layout()
fig.subplots_adjust(hspace=0)
fig.savefig("elements.png")
