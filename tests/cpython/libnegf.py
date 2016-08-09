from ctypes import *
from numpy.ctypeslib import ndpointer
from numpy import dtype

MAXCONT = 10
INTTYPE = 'int32'

class NEGF:


    class LNParams(Structure):
        _fields_ = [
            ("verbose", c_int),
            ("readoldsgf", c_int),
            ("spin", c_int),
            ("kpoint", c_int),
            ("iteration", c_int),
            ("g_spin", c_double),
            ("delta", c_double),
            ("dos_delta", c_double),
            ("eneconv", c_double),
            ("wght", c_double),
            ("ec", c_double),
            ("ev", c_double),
            ("emin", c_double),
            ("emax", c_double),
            ("estep", c_double),
            ("mu_n", c_double * MAXCONT),
            ("mu_p", c_double * MAXCONT),
            ("mu", c_double * MAXCONT),
            ("contact_dos", c_double * MAXCONT),
            ("fictcont", c_int * MAXCONT),
            ("kbt", c_double * MAXCONT),
            ("np_n", c_int * 2),
            ("np_p", c_int * 2),
            ("np_real", c_int * 11),
            ("n_kt", c_int),
            ("n_poles", c_int),
            ("ni", c_int * MAXCONT),
            ("nf", c_int * MAXCONT),
            ("dore", c_char)]


    def __init__(self, path = "./libnegf_x86_64.so"):
    
        #Initialize and store handler reference in self._href  
        self._lib = cdll.LoadLibrary("./libnegf_x86_64.so") 
        self._handler_size = c_int()
        self._lib.negf_gethandlersize(byref(self._handler_size))
        self._handler = (c_int * self._handler_size.value)()
        self._href = pointer(self._handler)
        self._href_type = POINTER(c_int * self._handler_size.value)
        self._lib.negf_init_session.argtypes = [self._href_type]
        self._lib.negf_init_session(self._href)
        self._lib.negf_init.argtypes = [self._href_type]
        self._lib.negf_init(self._href)

        #Init parameters to default
        self.params = NEGF.LNParams()
        self._lib.negf_get_params(self._href, byref(self.params))

    def __del__(self):
  
        #We only need to clean up the library
        self._lib.negf_destruct_libnegf(self._href)  

    def get_params(self):
        """
        Get parameters from libnegf instance and update
        the class member. For debug
        or to get default values. 
        """
        self._lib.negf_get_params.argtypes = [
                self._href_type,
                POINTER(NEGF.LNParams)
                ]
        self._lib.negf_get_params(self._href, pointer(self.params))

    def set_params(self):
        """
        Set the parameters from class member to libnegf.
        This is always called before a "solve" function
        """
        self._lib.negf_set_params.argtypes = [
                self._href_type,
                POINTER(NEGF.LNParams)
                ]
        self._lib.negf_set_params(self._href, byref(self.params))
    
    def read_input(self):
        """
        Parse negf.in
        """  
        self._lib.negf_read_input(self._href)

    def solve_landauer(self):
        """
        Solve the Landauer problem: calculate tunnelling and 
        (eventually) LDOS
        """
        self.set_params()
        self._lib.negf_solve_landauer(self._href)

    def read_hs(self, re_fname, im_fname, target):
        """
        Read H and S from file.

        Args:
            re_fname (string): real part path
            im_fname (string): string with imaginary part path
            target (int): 0 for hamiltonian, 1 for overlap
        """
        re_f = c_char_p(re_fname)
        im_f = c_char_p(im_fname)
        tt = c_int(target)
        self._lib.negf_read_hs.argtypes = [self._href_type, 
            c_char_p, c_char_p, c_int]
        self._lib.negf_read_hs(self._href, re_f, im_f, tt)

    def set_s_id(self, nrow):
        """
        Set the overlap matrix as identity matrix.

        Args:
            nrow (int): number of rows
        """
        self._lib.negf_set_s_id.argtypes = [self._href_type, c_int]
        self._lib.negf_set_s_id(self._href, c_int(nrow))

    def init_structure(self, ncont, contend, surfend, npl, plend, cblk):
        """
        Initialize the geometrical structure.

        Args:
            ncont (int): number of contacts
            contend (numpy.ndarray): end of contact indexes (fortran indexing)
            surfend (numpy.ndarray): end of surface indexes (fortran indexing)
            npl (int): number of PLs
            plend (numpy.ndarray): end of PL indexes (fortran indexing)
            cblk (numpy.ndarray): indexes of blocks interacting with contacts
                (fortran indexing)
        """
        self._lib.negf_init_structure.argtypes = [
                self._href_type,
                c_int,
                ndpointer(c_int),
                ndpointer(c_int),
                c_int,
                ndpointer(c_int),
                ndpointer(c_int)
                ]
        self._lib.negf_init_structure(self._href,
                c_int(ncont), 
                contend.astype(dtype=INTTYPE), 
                surfend.astype(dtype=INTTYPE),
                c_int(npl), 
                plend.astype(dtype=INTTYPE), 
                cblk.astype(dtype=INTTYPE))

    def write_tun_and_dos(self):
        """
        Write tunnelling and LDOS to file (for debugging)
        """
        self._lib.negf_write_tunneling_and_dos(self._href)

    def print_tnegf(self):
        """
        Write all infos on TNegf container, for debug
        """
        self._lib.negf_print_tnegf.argtypes = [self._href_type]
        self._lib.negf_print_tnegf(self._href)
      
