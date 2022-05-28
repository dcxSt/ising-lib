# TODO: implement this in rust... eventually
import csv
import numpy as np
import matplotlib.pyplot as plt

with open("data06_56.csv","r") as file:
    data = [[float(i) for i in line] for line in csv.reader(file)]
data = np.asarray(data)
temp,ydata = data[:,0],data[:,1:]
data_dict = {temp[idx]:{"samples":ydata[idx,:]} for idx in range(len(temp))} # key = temp, value = dict {samples, mean, var}
for idx,t in enumerate(temp):
    samples = data_dict[t]["samples"]
    mu = samples.mean() 
    var = np.power(samples - mu,2).mean()
    data_dict[t]["mu"] = mu
    data_dict[t]["var"] = var

# TODO:make this plot pretty, put it in the readme of this example
# TODO: first step is to get better time data, better resolution
plt.figure(figsize=(8,5))
plt.title("Nearest Neighbour Correlations\nMonte Carlo Sample")
for t,data in data_dict.items():
    plt.plot([t]*len(data["samples"]),data["samples"],"g.",alpha=0.2)
    plt.plot([t],data["mu"],"bx")
    plt.errorbar(x=[t],y=data["mu"],xerr=0,yerr=np.sqrt(data["var"]),color="red")

plt.xlabel("Temperature")
plt.ylabel("Nearest Neighbour Correlation")
plt.grid()

plt.tight_layout()

# plt.plot(temp,mu,"-",color="black",alpha=0.3,linewidth=0.5)
# plt.errorbar(x=temp,y=mu,xerr=np.zeros(len(temp)),yerr=np.sqrt(var),ecolor="black",elinewidth=0.5)


# x,y = np.arange(10),np.arange(10)
# plt.errorbar(x=x,y=x,xerr=np.zeros(10),yerr=np.ones(10),fmt="k-",linewidth=0.5,elinewidth=0.5)
# plt.plot(x,y,"r.")
plt.show(block=True)




