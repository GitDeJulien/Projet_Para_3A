import numpy as np
import matplotlib.pyplot as plt

time_value_D = []
time_value_R = []


#### NX #####

Nx_value = [10,20,40,60,80,100,150,200,250,300]
        
for i in Nx_value:
    with open(f'BC.Dirichlet.Nx.{i}.dat', 'r') as file:
        line = file.readlines()
        time_value_D.append(float(line[0]))
        
for i in Nx_value:
    with open(f'BC.Robin.Nx.{i}.dat', 'r') as file:
        line = file.readlines()
        time_value_R.append(float(line[0]))


Nx_value = np.array(Nx_value)
time_value_D = np.array(time_value_D)
time_value_R = np.array(time_value_R)

# # Affichage Nx
plt.figure()
plt.plot(Nx_value**2, time_value_D, '-ob', label='Dirichlet')
plt.plot(Nx_value**2, time_value_R, '-^r', label='Robin')
plt.title("Speed Up taille maillage avec 4 proc")
plt.xlabel("N (nombre de points)")
plt.ylabel("Time")
#plt.xscale('log')
#plt.xlim([-100, 100000])
plt.legend()
plt.savefig("Speed_Up_Nx_4proc_pc_local")

#### OVERLAP #####

# overlap_value = [2,3,4,5,6,7,8,9,10,11,12,13,14]


# for i in overlap_value:
#     with open(f'BC.Dirichlet.overlap.{i}.dat', 'r') as file:
#         line = file.readlines()
#         time_value_D.append(float(line[0]))
        
# for i in overlap_value:
#     with open(f'BC.Robin.overlap.{i}.dat', 'r') as file:
#         line = file.readlines()
#         time_value_R.append(float(line[0]))


# overlap = np.array(overlap_value)
# time_value_D = np.array(time_value_D)
# time_value_R = np.array(time_value_R)


# ## Affichage overlap
# plt.figure()
# plt.plot(overlap, time_value_D, '-ob', label='Dirichlet')
# plt.plot(overlap, time_value_R, '-^r', label='Robin')

# plt.title("Speed Up pour le recouvrement avec 4 proc")
# plt.xlabel("Recouvrement (nombre de lignes)")
# plt.ylabel("Time")
# #plt.xlim([overlap[0], overlap[-1]])
# plt.legend()
# plt.savefig("Speed_Up_recouvrement_4proc_pc_local")