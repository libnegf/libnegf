import ctypes

class NEGF:

  def __init__(self, path = "./libnegf_x86_64.so"):
    
    #Initialize and store handler reference in self._href  
    self._lib = ctypes.cdll.LoadLibrary("./libnegf_x86_64.so") 
    self._handler_size = ctypes.c_int()
    self._lib.negf_gethandlersize(ctypes.byref(self._handler_size))
    self._handler = (ctypes.c_int * self._handler_size.value)()
    self._href = ctypes.pointer(self._handler)
    self._lib.negf_init_session(self._href)
    self._lib.negf_init(self._href)

    
  def __del__(self):
  
    #We only need to clean up the library
    self._lib.negf_destruct_libnegf(self._href)  

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
    self._lib.negf_solve_landauer(self._href)

  def write_tun_and_dos(self):
    """
    Write tunnelling and LDOS to file (for debugging)
    """
    self._lib.negf_write_tunneling_and_dos(self._href)
    
