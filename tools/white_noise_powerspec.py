# power spectrum calculator by Francesco Sinigaglia

import numpy as np
import scipy as sp
from nbodykit.lab import *
from nbodykit import setup_logging, style
import nbodykit
import matplotlib.pyplot as plt


# **********************************************
# INPUT PARAMETERS
Nft1 = 128   # Number of cells per side of the box 1
Nft2 = 128   # Number of cells per side of the box 2                          
file_name_dm = 'Name of file 1'   # Name of tracers file
file_name_tr = 'Name of file 2'   # Name of tracers file          
Lbox = 500                    # Lenght of the side of the box, given in Mpc/h

outfolder = ''

correct_dm_for_shotnoise = 'no'   # Shotnoise correction DM
correct_tr_for_shotnoise = 'no'  # Shotnoise correction TR


convert_to_delta_dm = 'yes'
convert_to_delta_tr = 'yes'

dm_prec = 'double'
tr_prec = 'double'

write_power_spectrum_dm_to_file = 'no'
write_power_spectrum_tr_to_file = 'no'

dim = '1d'

# **********************************************
# **********************************************
# Function write_to_file
def write_ps_to_file(k, Pk, modes, field_property):
    filename =  outfolder + 'Pk_%s_Nft%d.txt' %(field_property, Nft)
    f = open(filename, 'w')
    f.write('#  k    Pk    modes \n')
    for i in range(len(k)):
        f.write(str(k[i]) + '   ' + str(Pk[i]) + '    ' + str(int(modes[i]/2)) + '\n')
    f.close()
        

# *************************************************************
# MAIN BODY
# Import the fields 
 
if dm_prec == 'single':
    file_bin_dm = open(file_name_dm, "rb")
    dm_array = np.fromfile(file_bin_dm, dtype=np.float32)

else:
    file_bin_dm = open(file_name_dm, "r")
    dm_array = np.fromfile(file_bin_dm, dtype=np.float64)

if tr_prec == 'single':
    file_bin_tr = open(file_name_tr, "rb")
    tr_array = np.fromfile(file_bin_tr, dtype=np.float32)
else:
    file_bin_tr = open(file_name_tr, "r")
    tr_array = np.fromfile(file_bin_tr, dtype=np.float64)


# *************************************************************
# START P(k) COMPUTATION

mean_dm = np.mean(dm_array)
mean_tr = np.mean(tr_array)

# Go to overdensities, if requested in the input params
if convert_to_delta_dm == 'yes':
    dm_array = dm_array/mean_dm - 1.
    
if convert_to_delta_tr == 'yes':
    tr_array = tr_array/mean_tr - 1.


# NOW COMPUTE THE POWER SPECTRA
# Define the fundamental mode
kf = 2*np.pi/Lbox   


# Compute the shot-noise 
shotnoise_dm = Lbox**3/np.sum(dm_array)
shotnoise_tr = Lbox**3/np.sum(tr_array)

# Convert the 1-d fields in 3x3 grids
mesh_dm = np.reshape(dm_array, (Nft1,Nft1,Nft1))
mesh_tr = np.reshape(tr_array, (Nft2,Nft2,Nft2))

# Convert the grids in meshes readable by the nbodykit methods
mesh_dm = nbodykit.source.mesh.array.ArrayMesh(mesh_dm, BoxSize=Lbox)
mesh_tr = nbodykit.source.mesh.array.ArrayMesh(mesh_tr, BoxSize=Lbox)

# Compute power spectrum of DM
r_dm = FFTPower(mesh_dm, mode=dim, dk=kf)
Pk_dm = r_dm.power
power_spectrum_dm = Pk_dm['power'].real 

# Compute power spectrum of TR
r_tr = FFTPower(mesh_tr, mode=dim, dk=kf)
Pk_tr = r_tr.power
power_spectrum_tr = Pk_tr['power'].real

#Rescale the P(k) to the mean densities
power_spectrum_dm = power_spectrum_dm
power_spectrum_tr = power_spectrum_tr

# Eventually correct for shotnoise 
if correct_dm_for_shotnoise == 'yes':  
    power_spectrum_dm = power_spectrum_dm - shotnoise_dm
power_spectrum_dm[0] = 1. # to avoid having negative values of P(k)

if correct_tr_for_shotnoise == 'yes':  
    power_spectrum_tr = power_spectrum_tr - shotnoise_tr
power_spectrum_tr[0] = 1. # to avoid having negative values of P(k)

# Define wavenumbers and modes
k_dm = Pk_dm['k']
modes_dm = Pk_dm['modes']

k_tr = Pk_tr['k']
modes_tr = Pk_tr['modes']

# Write the Pks to file
if write_power_spectrum_dm_to_file == 'yes':
    write_ps_to_file(k_dm, power_spectrum_dm, modes_dm, 'dm')
if write_power_spectrum_tr_to_file == 'yes':
    write_ps_to_file(k_tr, power_spectrum_tr, modes_tr, 'tracers')

# Plot the power spectra
plt.plot(k_dm, power_spectrum_dm, 'red', linestyle='solid', label='DM')
plt.plot(k_tr, power_spectrum_tr, 'green', linestyle='dashed', label='Tracers')
plt.xscale('log')
plt.yscale('log')
plt.xlabel('k')
plt.ylabel('P(k)')
plt.title('Power spectra')
plt.xlim([0.02, 0.85]) # Depend on mesh size and physical volume
plt.legend()
plt.show()


# Plot ratios
#plt.plot(k_dm, power_spectrum_tr/power_spectrum_dm, 'red', linestyle='solid')#, label='FGPA/Ref')
#plt.plot(k_dm, np.ones(len(k_dm)))
#plt.fill_between(k_tr, 0.99*np.ones(len(k_tr)), 1.01*np.ones(len(k_tr)),'gray', alpha=0.7)
#plt.fill_between(k_tr, 0.98*np.ones(len(k_tr)), 1.02*np.ones(len(k_tr)),'gray', alpha=0.5)
#plt.fill_between(k_tr, 0.95*np.ones(len(k_tr)), 1.05*np.ones(len(k_tr)),'gray', alpha=0.3)
#plt.xscale('log')
#plt.yscale('linear')
#plt.xlabel('k')
#plt.ylabel('P(k)')
#plt.title(r'Ratio P(k)$_{\rm{FGPA}}$/P(k)$_{\rm{ref}}$ %d³' %Nft)
#plt.xlim([0.1,16.08])
#plt.legend()
#plt.show()



print('   ')
print('Done!')

