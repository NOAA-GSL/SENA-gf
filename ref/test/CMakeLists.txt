list( APPEND gf_test_output
  test_output/output.test
)
 
add_subdirectory( tools )

# Create test_output directory for test output and symlink all files
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/test_output)
foreach(FILENAME ${gf_test_output})
    execute_process( COMMAND ${CMAKE_COMMAND} -E create_symlink
           ${CMAKE_CURRENT_SOURCE_DIR}/../../test/${FILENAME}
           ${CMAKE_CURRENT_BINARY_DIR}/${FILENAME} )
endforeach()

set(COMPARE_SH ${CMAKE_BINARY_DIR}/bin/gf_compare.sh)

# Add unit tests here.

# CPU kernel
add_test(NAME cpu_kernel
         COMMAND bash -c "${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1 ${MPIEXEC_PREFLAGS} ../src/gf_kernel_cpu > test_output/gf_kernel_cpu.log")
add_test(NAME compare_cpu_kernel
           COMMAND ${COMPARE_SH} test_output/gf_kernel_cpu.log test_output/output.test)
if ( CMAKE_Fortran_COMPILER_ID STREQUAL "NVHPC" )
  set_property(TEST cpu_kernel PROPERTY ENVIRONMENT "FORTRANOPT=no_minus_zero")
endif()

# Tests when GPU is enabled
if ( ENABLE_GPU )

  # GPU kernel version 1
  add_test(NAME gpu_kernel
           COMMAND bash -c "${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1 ${MPIEXEC_PREFLAGS} ../src/gf_kernel_gpu > test_output/gf_kernel_gpu.log")
  add_test(NAME compare_gpu_kernel
           COMMAND ${COMPARE_SH} test_output/gf_kernel_gpu.log test_output/output.test)
  if ( CMAKE_Fortran_COMPILER_ID STREQUAL "NVHPC" )
    set_property(TEST gpu_kernel PROPERTY ENVIRONMENT "FORTRANOPT=no_minus_zero")
  endif()

endif()
