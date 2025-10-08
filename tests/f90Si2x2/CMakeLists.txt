set(sources
	constants.F90
	readHS.F90
	matconv.F90
	test.F90)

setup_f90_mpi_test(f90Si2x2 "normal" ${sources}) 

unzip_test_data(f90Si2x2 "hs.tar.xz")

