####################################################################
# COMMON FLAGS
####################################################################
set( CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -g -traceback -Mnofma -Kieee" )

####################################################################
# RELEASE FLAGS
####################################################################

# Must turn off SIMD vectorization to get same results as for debug
set( CMAKE_Fortran_FLAGS_RELEASE "-fast -mp -Mnovect" )

####################################################################
# DEBUG FLAGS
####################################################################

set( CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} -O0 -Mbounds -Mchkptr -Mchkstk -Ktrap=fp" )
#set( CMAKE_Fortran_FLAGS_DEBUG "-O0" )

####################################################################
# FLAGS FOR GPU
####################################################################

set( Fortran_GPU_FLAGS        "" )
