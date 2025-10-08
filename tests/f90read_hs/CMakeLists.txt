set(sources
	constants.F90
	readHS.F90
	matconv.F90
	test.F90)

setup_f90_test(f90read_hs "normal unit" ${sources})
