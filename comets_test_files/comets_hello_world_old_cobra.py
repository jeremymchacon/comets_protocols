import cobra
import cometspy as c
import os
os.environ["GUROBI_COMETS_HOME"] = os.environ["GUROBI_HOME"]


import cobra.test
m = cobra.test.create_test_model("textbook")

m = c.model(m)
m.initial_pop = [0, 0, 0.01]

l = c.layout([m])
l.set_specific_metabolite("nh4_e", 1000)
l.set_specific_metabolite("pi_e", 1000)
l.set_specific_metabolite("glc__D_e", 1)

p = c.params()

sim = c.comets(l, p)
sim.run()
print(sim.run_output)
print("\n\ncheck hello_world_data.csv for the data")
sim.total_biomass.to_csv("hello_world_data.csv")
