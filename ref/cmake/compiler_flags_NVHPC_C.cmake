####################################################################
# COMMON FLAGS
####################################################################
set( CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -g")

####################################################################
# RELEASE FLAGS
####################################################################

set( CMAKE_C_FLAGS_RELEASE "-O3 -mp" )

####################################################################
# DEBUG FLAGS
####################################################################

set( CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -O0" )
