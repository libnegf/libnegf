setup_test(
  f90Si_nin test.F90
  MPI
  SOURCES constants.F90 matconv.F90 readHS.F90
  LABELS long
)
unzip_test_data(f90Si_nin "hs.tar.xz")
