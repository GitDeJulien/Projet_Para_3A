import numpy as np
import matplotlib.pyplot as plt
import sys

nb_tot_proc = sys.argv[1]
max_time_metis = []
max_time_scotch = []

xxproc = [1,2,3,4]

with open('metis/nbproc.1.plt', 'r') as file:
    lines = file.readlines()
    charge_tot_metis = lines[0]
    
with open('scotch/nbproc.1.plt', 'r') as file:
    lines = file.readlines()
    charge_tot_scotch = lines[0]
    
max_time_metis.append(float(charge_tot_metis))
max_time_scotch.append(float(charge_tot_scotch))


for i in range(2,int(nb_tot_proc)+1):
    with open(f'metis/nbproc.{i}.plt', 'r') as file:
        data = [float(line.strip()) for line in file]
        
    # Calculate the mean
    max_time_metis.append(max(data))
    
    with open(f'scotch/nbproc.{i}.plt', 'r') as file:
        data = [float(line.strip()) for line in file]
        
    # Calculate the mean
    max_time_scotch.append(max(data))
    
charge_tot_metis = float(charge_tot_metis)
charge_tot_scotch = float(charge_tot_scotch)


max_time_metis = np.array(max_time_metis)
max_time_scotch = np.array(max_time_scotch)
xxproc = np.array(xxproc)



plt.figure()
plt.plot(xxproc,1.*xxproc, ':k', label="Ideal")
plt.plot(xxproc,charge_tot_metis/max_time_metis, '-ob', label="Metis")
plt.plot(xxproc,charge_tot_scotch/max_time_scotch, '-^r', label="Scotch")
plt.xlabel("Number of proc")
plt.ylabel("Speed Up")

plt.title("Speed Up plot for personal computeur")
plt.legend()
plt.savefig("SpeedUp_pc_perso")

