from libnegf import NEGF
import numpy as np
import matplotlib.pyplot as plt

negf = NEGF()
negf.read_hs("HR.dat", "HI.dat", 0)
negf.set_s_id(100)
negf.init_structure(2, np.array([80,100]), np.array([60,80]), 
        2, np.array([30, 60]), np.array([2,1]))
negf.params.emin = -3.0
negf.params.emax = 3.0
negf.params.estep = 0.01
negf.set_params()
#For debug, uncomment to print TNegf container
#negf.print_tnegf()

# Set also some local DOS
negf.set_ldos_intervals(np.array([1, 31, 1]), np.array([60, 60, 30]))
negf.solve_landauer()

#Get transmission as numoy object
trans = negf.get_transmission()
ldos = negf.get_ldos()
plt.plot(trans[0,:])
plt.show()
plt.plot(ldos[0,:])
plt.plot(ldos[1,:])
plt.plot(ldos[2,:])
plt.show()

#Also write to file for debug
negf.write_tun_and_dos()
