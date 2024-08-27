from gf_state import GFState
import numpy as np

s = GFState(rkind=np.float64)
s.print_state("Input state")
s.write_state("input_state_python.nc")
s.read_state("input_state_python.nc")
s.print_state("Output state")
 
