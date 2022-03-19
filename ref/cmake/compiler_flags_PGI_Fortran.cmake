####################################################################
# COMMON FLAGS
####################################################################
set( CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ")

####################################################################
# RELEASE FLAGS
####################################################################

set( CMAKE_Fortran_FLAGS_RELEASE "-O3 -mp")

####################################################################
# DEBUG FLAGS
####################################################################

set( CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} -O0 -g" )

####################################################################
# FLAGS FOR GPU
####################################################################

set( Fortran_GPU_FLAGS        "-acc -Minfo=accel -ta=tesla:cc60,cuda10.1" )
#set( Fortran_GPU_FLAGS        "" )
