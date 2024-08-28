from gf_state import GFState
import numpy as np

s = GFState(rkind=np.float64)
s.print_state("Input state")
# Run driver
s.write_state("output_state_python.nc")
s.print_state("Output state")
 
