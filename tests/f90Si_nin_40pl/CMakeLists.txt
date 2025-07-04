set(sources
	constants.F90
	readHS.F90
	matconv.F90
	test.F90)

setup_f90_mpi_test(f90Si_nin_40pl "long" ${sources})

unzip_test_data(f90Si_nin_40pl "hs.tar.xz")

