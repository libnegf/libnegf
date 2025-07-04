set(sources 
    main.f90
    random.f90)

setup_f90_test(testCUDA_decimation "normal cuda" ${sources})
