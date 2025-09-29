from gf_state import GFState
from cu_gf_driver import GFDriver

errmsg = ""
errflg = 0

first_step = 1
last_step = 576
# first_step = 576
# last_step = 576

first_input = f"data/input_state_{first_step:04d}.nc"
state = GFState()
state.read_state(first_input)
driver = GFDriver(state)


for n in range(first_step, last_step + 1):
    print(f"Running GF test for state {n:04d}")

    step_input = f"data/input_state_{n:04d}.nc"
    state.read_state(step_input)

    # Call GF driver
    # driver = GFDriver(state)
    driver.cu_gf_driver_run(state, errmsg, errflg)

    step_output = f"output_state_{n:04d}.nc"
    state.write_state(step_output)
