from libnegf import NEGF
import numpy as np
import matplotlib.pyplot as plt
import scipy

negf = NEGF()

######################################################################
# Assign matrices: you can either load from file or pass a csr matrix
# Uncomment the correspoinding section for debug
###################################
# Load data
# negf.read_hs("HR.dat", "HI.dat", 0)
###################################
# Explicit assignment
mat = np.zeros(shape=(100,100), dtype='complex128')
for ii in range(80):
    mat[ii,ii-1] = 1.0
    mat[ii-1,ii] = 1.0
for ii in range(81,100):
    mat[ii,ii-1] = 1.0
    mat[ii-1,ii] = 1.0
mat[0,80] = 1.0
mat[80,0] = 1.0

# I'm going to add a double barrier, just for fun!
mat[25,25] = 1.0
mat[35,35] = 1.0

mat_csr = scipy.sparse.csr_matrix(mat)
mat_csr.sort_indices()
negf.set_h(mat_csr)
###################################
negf.set_s_id(100)
negf.init_structure(2, np.array([80,100]), np.array([60,80]), 
        4, np.array([15, 30, 45, 60]), np.array([4,1]))
negf.params.emin = -3.0
negf.params.emax = 3.0
negf.params.estep = 0.01
negf.params.mu[0] = 0.1
negf.set_params()
#For debug, uncomment to print TNegf container
#negf.print_tnegf()

negf.solve_density()
dm = negf.get_dm()
plt.imshow(np.real(dm.todense()))
plt.show()

# Set also some local DOS
negf.set_ldos_intervals(np.array([1, 31, 1]), np.array([60, 60, 30]))
negf.solve_landauer()

#Get transmission, dos and energies as numpy object
energies, foo = np.real(negf.get_energies())
trans = negf.get_transmission()
ldos = negf.get_ldos()
currents = negf.get_currents()
print('Currents',currents)
plt.plot(energies, trans[0,:])
plt.show()
plt.plot(energies, ldos[0,:])
plt.plot(energies, ldos[1,:])
plt.plot(energies, ldos[2,:])
plt.show()

#Also write to file for debug
negf.write_tun_and_dos()
