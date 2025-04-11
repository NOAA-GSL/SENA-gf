from gf_state import GFState

input_file = "input_state_0400.nc"
output_file = "test_output.nc"

state = GFState()

state.read_state(input_file)
    

print(state.ntracer  )
print(state.im)
print(state.km)

state.write_state(output_file)
 
