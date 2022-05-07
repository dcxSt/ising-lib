import csv
import numpy as np
import matplotlib.pyplot as plt

with open("data.csv","r") as file:
    data = [[float(i) for i in line] for line in csv.reader(file)]
data = np.asarray(data)
mu,var,temp = data[:,0],data[:,1],data[:,2]

plt.figure(figsize=(8,5))
plt.plot(temp,mu,"r.")
# plt.plot(temp,mu,"-",color="black",alpha=0.3,linewidth=0.5)
plt.errorbar(x=temp,y=mu,xerr=np.zeros(len(temp)),yerr=np.sqrt(var),ecolor="black",elinewidth=0.5)

# x,y = np.arange(10),np.arange(10)
# plt.errorbar(x=x,y=x,xerr=np.zeros(10),yerr=np.ones(10),fmt="k-",linewidth=0.5,elinewidth=0.5)
# plt.plot(x,y,"r.")
plt.show(block=True)

