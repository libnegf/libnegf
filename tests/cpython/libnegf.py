from ctypes import *
import numpy as np
from numpy.ctypeslib import ndpointer
from numpy import dtype
from scipy.sparse import csr_matrix

MAXCONT = 10
MAXLAYERS = 10000
INTTYPE = 'int32'
REALTYPE = 'float64'

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
            ("ni", c_int * MAXLAYERS),
            ("nf", c_int * MAXLAYERS),
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

    def solve_density(self):
        """
        Solve the density problem for an all electron problem (dft like)
        """
        self.set_params()
        self._lib.negf_solve_density_dft(self._href)

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
                contend.astype(dtype=INTTYPE, copy=False), 
                surfend.astype(dtype=INTTYPE, copy=False),
                c_int(npl), 
                plend.astype(dtype=INTTYPE, copy=False), 
                cblk.astype(dtype=INTTYPE, copy=False))


    def set_h(self, mat):
        """
        Set H from a scipy.sparse.csr_matrix
        NOTE: libnegf is picky about the order of the PL blocks
        in the sparse matrix as well. There is no automatic reordering,
        you should have a well-sorted matrix.
        
        Args:
            mat (complex csr_matrix): input Hamiltonian
        """
        self._lib.negf_set_h_cp.argtypes = [
            self._href_type,
            c_int,
            ndpointer(c_double),
            ndpointer(c_double),
            ndpointer(c_int),
            ndpointer(c_int)]
        mat_re = np.array(np.real(mat.data))
        mat_im = np.array(np.imag(mat.data))
        self._lib.negf_set_h_cp(self._href,
            c_int(mat.shape[0]),
            mat_re,
            mat_im,
            mat.indices + 1,
            mat.indptr + 1
        )
        
    def set_s(self, mat):
        """
        Set S from a scipy.sparse.csr_matrix
        NOTE: libnegf is picky about the order of the PL blocks
        in the sparse matrix as well. There is no automatic reordering,
        you should have a well-sorted matrix.
        
        Args:
            mat (complex csr_matrix): input Overlap
        """
        self._lib.negf_set_s_cp.argtypes = [
            self._href_type,
            c_int,
            ndpointer(c_double),
            ndpointer(c_double),
            ndpointer(c_int),
            ndpointer(c_int)]
        mat_re = np.array(np.real(mat.data))
        mat_im = np.array(np.imag(mat.data))
        self._lib.negf_set_s_cp(self._href,
            c_int(mat.shape[0]),
            mat_re,
            mat_im,
            mat.indices + 1,
            mat.indptr + 1
        )

    def get_energies(self):
        """
        Get a local copy of energies array

        Returns:
            real_en (array): real part of points on 
            energy axis. This quantity may change in libnegf in 
            runtime when performing different integrals (contour, 
            real axis)
            im_en (array): imaginary part (same as above)
        """
        self._lib.negf_get_energies.argtypes = [
                self._href_type,
                POINTER(c_int),
                ndpointer(c_double),
                ndpointer(c_double),
                c_int
                ]
        npoints = c_int()
        self._lib.negf_get_energies(self._href,
                byref(npoints), 
                np.zeros(1,dtype=REALTYPE), 
                np.zeros(1,dtype=REALTYPE), 0)
        re_en = np.zeros(npoints.value, dtype=REALTYPE)
        im_en = np.zeros(npoints.value, dtype=REALTYPE)
        self._lib.negf_get_energies(self._href,
                byref(npoints), re_en, im_en, 1)
        return re_en, im_en

    
    def get_currents(self):
        """
        Get a local copy of currents array

        Returns:
            currents (array): array of currents for
            each possible lead pair defined in input
        """
        self._lib.negf_get_currents.argtypes = [
                self._href_type,
                POINTER(c_int),
                ndpointer(c_double),
                c_int
                ]
        npoints = c_int()
        self._lib.negf_get_currents(self._href,
                byref(npoints), 
                np.zeros(1,dtype=REALTYPE), 0)
        currents = np.zeros(npoints.value, dtype=REALTYPE)
        self._lib.negf_get_currents(self._href,
                byref(npoints), currents, 1)
        return currents
    
    
    def get_dm(self):
        """
        Get a local copy of CSR sparse density matrix

        Returns:
            dm (scipy sparse): density matrix
        """
        self._lib.negf_get_dm.argtypes = [
                self._href_type,
                POINTER(c_int),
                POINTER(c_int),
                ndpointer(c_int),
                ndpointer(c_int),
                ndpointer(c_double),
                ndpointer(c_double),
                c_int
                ]
        nnz = c_int()
        nrow = c_int()
        self._lib.negf_get_dm(self._href,
                byref(nnz),
                byref(nrow),
                np.zeros(1,dtype=INTTYPE), 
                np.zeros(1,dtype=INTTYPE), 
                np.zeros(1,dtype=REALTYPE), 
                np.zeros(1,dtype=REALTYPE), 0)
        rowpnt =np.zeros(nrow.value + 1, dtype=INTTYPE)
        colind =np.zeros(nnz.value, dtype=INTTYPE)
        re_dm = np.zeros(nnz.value, dtype=REALTYPE)
        im_dm = np.zeros(nnz.value, dtype=REALTYPE)
        self._lib.negf_get_dm(self._href,
                byref(nnz), 
                byref(nrow),
                rowpnt,
                colind,
                re_dm,
                im_dm, 1)
        #Fix indexing
        rowpnt = rowpnt - 1
        colind = colind - 1
        dm = csr_matrix((re_dm + 1j*im_dm, colind, rowpnt), dtype='complex128')
        return dm



    def get_transmission(self):
        """
        Get a local copy of transmission from libnegf

        Returns:
            trans (ndarray): transmission for all possible
                lead pairs (2D array). Each row contains 
                the result for a lead pair (npair, values).
        """
        self._lib.negf_associate_transmission.argtypes = [
                self._href_type,
                POINTER(c_int * 2),
                POINTER(POINTER(c_double))
                ]
        tr_pointer = POINTER(c_double)()
        tr_shape = (c_int * 2)()
        self._lib.negf_associate_transmission(self._href,
                pointer(tr_shape),
                pointer(tr_pointer))
        tr_shape = (tr_shape[1] , tr_shape[0])
        trans = (np.ctypeslib.as_array(tr_pointer, shape=tr_shape)).copy()
        return trans


    def get_ldos(self):
        """
        Get a local copy of dos from libnegf

        Returns:
            ldos (ndarray): local DOS for all given orbital intervals
                (2D array). Each row contains the result for
                an interval (ninterval, values)
        """
        self._lib.negf_associate_ldos.argtypes = [
                self._href_type,
                POINTER(c_int * 2),
                POINTER(POINTER(c_double))
                ]
        ldos_pointer = POINTER(c_double)()
        ldos_shape = (c_int * 2)()
        self._lib.negf_associate_ldos(self._href,
                pointer(ldos_shape),
                pointer(ldos_pointer))
        ldos_shape = (ldos_shape[1] , ldos_shape[0])
        ldos = (np.ctypeslib.as_array(ldos_pointer, shape=ldos_shape)).copy()
        return ldos


    def set_ldos_intervals(self, istart, iend):
        """
        Define intervals for LDOS calculations

        Args:
            istart (int array): starting orbitals (fortran indexing)
            iend (int array): ending orbitals (fortran indexing)
        """
        nldos = istart.size
        self._lib.negf_init_ldos(self._href, c_int(nldos))
        self._lib.negf_set_ldos_intervals.argtypes = [
                self._href_type,
                c_int,
                ndpointer(c_int),
                ndpointer(c_int)]
        self._lib.negf_set_ldos_intervals(self._href,
                nldos,
                istart.astype(dtype=INTTYPE, copy=False),
                iend.astype(dtype=INTTYPE, copy=False))


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
      
