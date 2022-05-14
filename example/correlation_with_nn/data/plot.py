# TODO: implement this in rust... eventually
import csv
import numpy as np
import matplotlib.pyplot as plt

with open("data.csv","r") as file:
    data = [[float(i) for i in line] for line in csv.reader(file)]
data = np.asarray(data)
temp,ydata = data[:,0],data[:,1:]
mu = ydata.mean(axis=1)
# TODO: compute variance
# var = 

# TODO:make this plot pretty, put it in the readme of this example
# TODO: first step is to get better time data, better resolution
plt.figure(figsize=(8,5))
for idx in range(ydata.shape[1]):
    plt.plot(temp,ydata[:,idx],"g.")
plt.plot(temp,mu,"ro")
# plt.plot(temp,mu,"-",color="black",alpha=0.3,linewidth=0.5)
# plt.errorbar(x=temp,y=mu,xerr=np.zeros(len(temp)),yerr=np.sqrt(var),ecolor="black",elinewidth=0.5)


# x,y = np.arange(10),np.arange(10)
# plt.errorbar(x=x,y=x,xerr=np.zeros(10),yerr=np.ones(10),fmt="k-",linewidth=0.5,elinewidth=0.5)
# plt.plot(x,y,"r.")
plt.show(block=True)

