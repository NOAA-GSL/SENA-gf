from gf_state import GFState
from cu_gf_driver import cu_gf_driver_run

input_file = "input_state_0400.nc.baseline"
output_file = "test_output.nc"
errmsg = ""
errflg = 0

state = GFState()

state.read_state(input_file)

# Call GF driver
cu_gf_driver_run(state, errmsg, errflg)

state.write_state(output_file)
