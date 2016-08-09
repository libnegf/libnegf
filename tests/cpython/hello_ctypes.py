from libnegf import NEGF
import numpy as np

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

negf.solve_landauer()
negf.write_tun_and_dos()
