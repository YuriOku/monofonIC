# %%
import numpy as np
import camb
from camb import model, initialpower

# Set parameters
class param:
  z = 49
  kmax = 1e3
  kmin = 1e-4
  npoints = 1000
  cosmology = "Planck18"
  nonlinear=True

# Set up cosmologal parameters
class Cosmo:
  def __init__(self, cosmology):
    if cosmology == "WMAP7":
      self.H0 = 70.2
      self.ombh2 = 0.02255
      self.omch2 = 0.1126
      self.As = 2.430e-9
      self.ns = 0.968
    elif cosmology == "WMAP9":
      self.H0 = 69.33
      self.ombh2 = 0.02266
      self.omch2 = 0.1157
      self.As = 2.427e-9
      self.ns = 0.971
    elif cosmology == "Planck15":
      self.H0 = 67.74
      self.ombh2 = 0.02230
      self.omch2 = 0.1188
      self.As = 2.142e-9
      self.ns = 0.9667
    elif cosmology == "Planck18":
      self.H0 = 67.66
      self.ombh2 = 0.02242
      self.omch2 = 0.11933
      self.As = 2.105e-9
      self.ns = 0.9665
    else:
      pass

cosmo = Cosmo(param.cosmology)

# Set up parameters for CAMB
pars = camb.CAMBparams()
pars.set_cosmology(H0=cosmo.H0, ombh2=cosmo.ombh2, omch2=cosmo.omch2)
pars.InitPower.set_params(As=cosmo.As, ns=cosmo.ns)
pars.set_matter_power(redshifts=[param.z], nonlinear=param.nonlinear, kmax=param.kmax)

# Linear spectra
results = camb.get_results(pars)
kh, z, pk = results.get_matter_power_spectrum(minkh=param.kmin, maxkh=param.kmax, npoints = param.npoints)

# # # Plot 
# import matplotlib.pyplot as plt
# for i, (redshift, line) in enumerate(zip(z,['-','--'])):
#     plt.loglog(kh, pk[i,:], color='k', ls = line)
# plt.xlabel('k/h Mpc')
# plt.ylabel('P(k)')
# plt.show()
# plt.close()

# save to file
powerspec = np.vstack((kh, pk))
np.savetxt("powerspec_z{}_{}.txt".format(param.z, param.cosmology), powerspec.T)
pkm = 8 * np.pi**3 * pk
powerspec_music = np.vstack((kh, pkm, pkm, pkm, pkm, pkm, pkm, pkm, pkm, pkm, pkm, pkm, pkm))
np.savetxt("powerspec_music_z{}_{}.txt".format(param.z, param.cosmology), powerspec_music.T)
delta2 = 4 * np.pi * kh**3 * pk / (2 * np.pi)**3
powerspec_ngenic = np.log10(np.vstack((kh, delta2)))
np.savetxt("powerspec_ngenic_z{}_{}.txt".format(param.z, param.cosmology), powerspec_ngenic.T)
# %%
