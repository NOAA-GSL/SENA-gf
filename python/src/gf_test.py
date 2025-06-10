from gf_state import GFState
from cu_gf_driver import GFDriver

input_file = "input_state_0400.nc.baseline"
output_file = "test_output.nc"
errmsg = ""
errflg = 0

state = GFState()

for n in range(1,577):
# for n in range(76,77):
    print(f"Running GF test for state {n:04d}")

    input_file = f"data/input_state_{n:04d}.nc"
    state.read_state(input_file)

    # Call GF driver
    driver = GFDriver(state)
    driver.cu_gf_driver_run(state, errmsg, errflg)

    output_file = f"output_state_{n:04d}.nc"
    state.write_state(output_file)
