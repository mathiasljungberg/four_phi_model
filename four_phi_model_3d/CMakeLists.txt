cmake_minimum_required(VERSION 2.4)

PROJECT(Heff_fortran Fortran C CXX)

IF(CMAKE_Fortran_COMPILER MATCHES ifort)

  MESSAGE(STATUS  "ifort compiler found")

  SET(FLAGS "-O0 -warn all -g -CB")
  SET(INCDIR "-I/usr/include")
  SET(OMP_FLAGS "")
  SET(FFTW_LIB "-L/usr/lib -lfftw3 -lm")
  SET(MKLPATH "/opt/intelcompiler/mkl/lib/intel64")
  SET(LAPACK_LIB "-Wl,--start-group ${MKLPATH}/libmkl_intel_lp64.a ${MKLPATH}/libmkl_sequential.a ${MKLPATH}/libmkl_core.a -Wl,--end-group")

ELSE(CMAKE_Fortran_COMPILER MATCHES gfortran)

  MESSAGE(STATUS  "gfortran compiler found")

  IF(USER_OS STREQUAL mac)

    MESSAGE(STATUS  "option -DUSER_OS=mac, gfortran for mac")

    # gfortran/mac
    SET(FLAGS "-O0  -fbounds-check -ffree-line-length-none -framework vecLib")
    #SET(FLAGS "-O0  -fbounds-check -framework vecLib -DUSE_OWN_CZDOT")  
    SET(INCDIR "-I/usr/local/include -I/usr/include")
    SET(INCDIR "-I/usr/local/include -I/usr/include")
    SET(OMP_FLAGS "")
    SET(FFTW_LIB "-L/usr/local/lib -lfftw3 -lm")
    SET(LAPACK_LIB "")

  ELSE(USER_OS STREQUAL mac)

    MESSAGE(STATUS  "gfortran for Linux")

    # gfortran/linux
    SET(FLAGS "-O0  -fbounds-check -ffree-line-length-none -fmax-errors=5")
    SET(INCDIR "-I/usr/include")
    SET(OMP_FLAGS "")
    SET(FFTW_LIB "-L/usr/lib -lfftw3 -lm")
    SET(MKLPATH "")
    SET(LAPACK_LIB "-L/usr/lib -lblas -llapack")

  ENDIF(USER_OS STREQUAL mac)

ENDIF(CMAKE_Fortran_COMPILER MATCHES ifort)

SET(CMAKE_Fortran_FLAGS "${FLAGS}")
LINK_DIRECTORIES("/Users/mathias/work/fortran-utilities_build")
INCLUDE_DIRECTORIES("/Users/mathias/work/fortran-utilities_build")

#EXTERNAL_PROJECT_ADD("/Users/mathias/work/fortran-utilities_build/m_FFT.o")

SET(ALL_MODULES 
m_averages_func.F90             m_md_parameters.F90             m_symmetry.F90                  m_system_3d_scalar.F90
m_bandstructure.F90             m_system_3d.F90                 
hist_class.F90                  m_linalg.F90                    m_moment_fitting.F90            m_system_3d_mc.F90              parameters.F90
m_TCF.F90                       m_mc_parameters.F90             m_system_3d_mc_scalar.F90       
m_averages.F90                  m_mc_utils.F90                  m_strings.F90                   m_system_3d_md.F90
)

SET(EXT_OBJECTS "/Users/mathias/work/fortran-utilities_build/m_FFT.o" "/Users/mathias/work/fortran-utilities_build/sffteu.o")

SET(EXECUTABLES four_phi_model_3d bandstructure)

FOREACH(EXEC ${EXECUTABLES})
  ADD_EXECUTABLE(${EXEC} ${EXEC}.F90 ${EXT_OBJECTS} ${ALL_MODULES})
  TARGET_LINK_LIBRARIES(${EXEC} ${FLAGS} 
    ${LINK_ADDONS})
ENDFOREACH(EXEC)
