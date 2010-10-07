module lib_param

  use precision, only : dp
  use globals
  use mat_def
  use structure

  implicit none
  private

  public :: Tnegf


  type Tnegf
     
     integer :: verbose

     character(LST) :: file_re_H
     character(LST) :: file_im_H
     character(LST) :: file_re_S
     character(LST) :: file_im_S
     character(LST) :: file_struct
     real(dp) :: mu_n
     real(dp) :: mu_p
     real(dp) :: Ec
     real(dp) :: Ev
     real(dp) :: DeltaEc
     real(dp) :: DeltaEv     
     real(dp) :: E
     real(dp) :: dos
     real(dp) :: eneconv

     type(z_CSR) :: H
     type(z_CSR) :: S
     type(z_CSR) :: Gr
     type(z_CSR) :: rho
     type(z_CSR) :: rho_eps
     logical    :: isSid

     type(TStruct_Info) :: str

     real(dp) :: delta
     real(dp) :: Emin
     real(dp) :: Emax
     real(dp) :: Estep
     real(dp) :: Temp

     integer :: Np_n(2)
     integer :: Np_p(2)
     integer :: n_kt
     integer :: n_poles
     real(dp) :: spin     

  end type Tnegf

end module lib_param
