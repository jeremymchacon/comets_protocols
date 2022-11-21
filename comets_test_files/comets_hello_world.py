import cobra
import cometspy as c
import os
os.environ["GUROBI_COMETS_HOME"] = os.environ["GUROBI_HOME"]

m = cobra.io.load_model("textbook")

# if line 6 doesnt work, you might have an older version of cobra, and instead do the following:
#
# import cobra.test
# m = cobra.test.create_test_model("textbook")

m = c.model(m)
m.initial_pop = [0, 0, 0.01]

l = c.layout([m])
l.set_specific_metabolite("nh4_e", 1000)
l.set_specific_metabolite("pi_e", 1000)
l.set_specific_metabolite("glc__D_e", 1)

p = c.params()

sim = c.comets(l, p)
sim.run()
sim.total_biomass.head()
