Overall the functions depend on the following R packages ("zoo","sf","TUWmodel","doParallel","foreach","hydroGOF")

-------------------------------------------------- PASS FUNCTIONS LIST --------------------------------------------------


PASS_GAR_TUWGAR.R: function to apply HydroPASS for the specific case of the TUW_GAR model

PASS_ME_functions.R: examples of objective functions to be passed to HydroPASS (examples working for the TUW_GAR application)

PASS.options.R: function to define options to be passed to HydroPASS

print_PASS: function to print HydroPASS outputs

.DT.app.R: random forest


-------------------------------------------------- ACCESSORY FUNCTIONS --------------------------------------------------


HU_topology.R: function to pre-calculate intersections between climate input pixels and hydrological units (HUs)

GAR_level.R: function to calculate the "Level" of each HU in the river network

Z_temp_pre.R: function to calculate the flow-chains between the HUs


-------------------------------------------------- SIMULATION FUNCTION --------------------------------------------------


TUWGAR_sim.R: function to simulate with TUWGAR over the entire GAR domain given parameters for all HUs, HUs shapes and connections and input. 



