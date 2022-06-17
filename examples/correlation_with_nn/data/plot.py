#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np
import csv
import matplotlib.pyplot as plt


# In[4]:


with open("data08.csv","r") as file:
    data = [[float(i) for i in line] for line in csv.reader(file)]
    data = np.asarray(data)
    temps,data = data[:,0],data[:,1:]


# In[5]:


plt.figure(figsize=(9,9))
for temp,y in zip(temps,data):
    plt.plot([temp]*len(y), y,"x",label="{:.2}".format(temp))
plt.legend()
plt.title("Nearest Neighbour Corrleation",fontsize=24)
plt.ylabel("Monte-Carlo estimated correlation",fontsize=18)
plt.xlabel("Temperature in units of KbT (Boltzman's constant)", fontsize=18)
plt.show()


# The data collected was sampled 10 times per ising simulation. If we like we can take the mean of the samples in each simulation by averaging over those 10 samples.

# In[6]:


# Average the data across each simulation
assert data.shape[1] % 10 == 0
n_sims_per_temp = data.shape[1] // 10
data_avg = data.reshape((-1,n_sims_per_temp,10)).mean(axis=-1)


# The error in the mean is $\sigma / \sqrt{n-1}$, where $n$ is our `n_sims_per_temp` variable defined above. 

# In[7]:


mean = data_avg.mean(axis=-1)
err  = data_avg.std(axis=-1)/np.sqrt(n_sims_per_temp-1)


# In[8]:


plt.figure(figsize=(9,9))
for temp,y in zip(temps,data_avg):
    plt.plot([temp]*len(y), y,"x",label="{:.2}".format(temp))
plt.errorbar(temps,mean,yerr=err,xerr=np.zeros(err.shape),fmt="k.",
            barsabove=False,elinewidth=2)
plt.legend()
plt.title("Nearest Neighbour Corrleation",fontsize=24)
plt.ylabel("Monte-Carlo estimated correlation, avg over run",fontsize=18)
plt.xlabel("Temperature in units of KbT (Boltzman's constant)", fontsize=18)
plt.show()


# The error bars are quite small. 
