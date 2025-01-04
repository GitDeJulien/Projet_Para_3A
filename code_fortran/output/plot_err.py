import numpy as np
import matplotlib.pyplot as plt

# Assuming the file is named "data.txt"

# Load the file and extract the data
data = np.loadtxt("err_D.dat")

# Split the data into three separate arrays
temps_value = data[:, 1]  # Second column
err_value_D = data[:, 2]  # Third column

# Load the file and extract the data
data = np.loadtxt("err_R.dat")

# Split the data into three separate arrays
err_value_R = data[:, 2]  # Third column


# # Affichage error
plt.figure()
plt.plot(temps_value, err_value_D, '-b', label='Dirichlet')
plt.plot(temps_value, err_value_R, '-r', label='Robin')
plt.title("Evolution de l'erreur (problème instationnaire)")
plt.xlabel("Temps")
plt.ylabel("Erreur")
plt.yscale('log')
#plt.xlim([-100, 100000])
plt.legend()
plt.savefig("error_pb_3")