# divide ic file created by MUSIC2 into several files. 
# assume that ic file contains only PartType 0 and 1 and the number of particles of each component is identical.

import numpy as np
import h5py
import sys

FILENAME = sys.argv[1]
NFILES = int(sys.argv[2])

print("we will divide {}.hdf5 into {} files.".format(FILENAME, NFILES))
print("loading...")

infile = h5py.File("{}.hdf5".format(FILENAME), 'r')
NumPart_Total = infile["Header"].attrs["NumPart_Total"]

print("NumPart_Total: {}".format(NumPart_Total))

reminder = NumPart_Total[0] % NFILES
quotient = (NumPart_Total[0] - reminder) / NFILES
begin = [int(quotient * j) for j in range(NFILES+1)]
begin[-1] += reminder

print("Divided files will contain 2 x {} particles".format(quotient))

sys.stdout.flush()

# %% write particle data

for file_i in range(NFILES):
  ic = h5py.File("{}.{}.hdf5".format(FILENAME, file_i), 'w')

  print("writing {}.{}.hdf5".format(FILENAME, file_i), flush=True)

  PartType = [0, 1]
  for type_i in PartType:
    print("copy data of PartType {} from {} to {}".format(type_i, begin[file_i], begin[file_i+1]))
    # create group
    ic.create_group('PartType{}'.format(type_i))

    # copy data
    tmp_data = np.array(infile['PartType{}/ParticleIDs'.format(type_i)])[begin[file_i]:begin[file_i+1]]
    ic.create_dataset("PartType{}/ParticleIDs".format(type_i), data = tmp_data)
    print("IDs done", flush=True)

    tmp_data = np.array(infile['PartType{}/Masses'.format(type_i)])[begin[file_i]:begin[file_i+1]]
    ic.create_dataset("PartType{}/Masses".format(type_i), data = tmp_data)
    print("Masses done", flush=True)
    
    tmp_data = np.array(infile['PartType{}/Coordinates'.format(type_i)])[begin[file_i]:begin[file_i+1]]
    ic.create_dataset("PartType{}/Coordinates".format(type_i), data = tmp_data)
    print("Coordinates done", flush=True)
    
    tmp_data = np.array(infile['PartType{}/Velocities'.format(type_i)])[begin[file_i]:begin[file_i+1]]
    ic.create_dataset("PartType{}/Velocities".format(type_i), data = tmp_data)
    print("Velocities done", flush=True)

    ic.flush()
  
  ic.create_group('Header')

  print("copy Header")

  ic["Header"].attrs["BoxSize"] = infile["Header"].attrs["BoxSize"]
  ic["Header"].attrs["Flag_Cooling"] = infile["Header"].attrs["Flag_Cooling"]
  ic["Header"].attrs["Flag_Entropy_ICs"] = infile["Header"].attrs["Flag_Entropy_ICs"]
  ic["Header"].attrs["Flag_Feedback"] = infile["Header"].attrs["Flag_Feedback"]
  ic["Header"].attrs["Flag_Metals"] = infile["Header"].attrs["Flag_Metals"]
  ic["Header"].attrs["Flag_Sfr"] = infile["Header"].attrs["Flag_Sfr"]
  ic["Header"].attrs["Flag_StellarAge"] = infile["Header"].attrs["Flag_StellarAge"]
  ic["Header"].attrs["HubbleParam"] = infile["Header"].attrs["HubbleParam"]
  ic["Header"].attrs["MassTable"] = infile["Header"].attrs["MassTable"]
  ic["Header"].attrs["Omega0"] = infile["Header"].attrs["Omega0"]
  ic["Header"].attrs["OmegaLambda"] = infile["Header"].attrs["OmegaLambda"]
  ic["Header"].attrs["Redshift"] = infile["Header"].attrs["Redshift"]
  ic["Header"].attrs["Time"] = infile["Header"].attrs["Time"]
  ic["Header"].attrs["NumPart_Total"] = infile["Header"].attrs["NumPart_Total"]
  ic["Header"].attrs["NumPart_Total_HighWord"] = infile["Header"].attrs["NumPart_Total_HighWord"]
  
  ic["Header"].attrs["NumFilesPerSnapshot"] = NFILES
  ic["Header"].attrs["NumPart_ThisFile"] = [len(tmp_data), len(tmp_data), 0, 0, 0, 0]

  ic.flush()
  ic.close()

  print("{}.{}.hdf5 closed".format(FILENAME, file_i))

  sys.stdout.flush()

infile.close()

print("finish.")
