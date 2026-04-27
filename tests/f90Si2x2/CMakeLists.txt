setup_test(
  f90Si2x2 test.F90
  MPI
  SOURCES constants.F90 matconv.F90 readHS.F90
  LABELS long
)
unzip_test_data(f90Si2x2 "hs.tar.xz")
