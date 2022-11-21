import cobra
import os
import numpy as np
from copy import copy
import re
import pandas as pd
import sys
import random
import math
import cometspy as c
pd.set_option('display.max_rows', None)
os.environ["GUROBI_COMETS_HOME"] = "/opt/gurobi952/linux64"

### THe command line will determine which seed and carbon source. 
### the carbon source will repeat (0 = csn, 1 = thr, 2 = tre, 3 = csn, 4 = thr, ...)
### for each repeat, the spatial seed will increase
###     python chosen4_multitoxin_commandline.py 1    # this would do threonine with spatial seed 0
###     python chosen4_multitoxin_commandline.py 4    # this would do threonine with spatial seed 1

### this expects that the four models are in a subdirectory called "models"

### There are many parameters to play with. I've placed especially relevant ones here, though one could adjust many othersw

model_dir = "./models/" # the relative directory where the models live

# toxin effect parameters
# first, the effects of F_toxin and G_toxin on strain B:
km_BF = 2.e-5 # km_xy is the km for toxin y on species x
km_BG = 2.e-5
hill_BF = 5
hill_BG = 5
km_FB = 2.e-5
km_FG = 2.e-5
hill_FB = 5
hill_FG = 5
km_GB = 1.e-3
hill_GB = 5
km_IB = 1
hill_IB = 5
#spatial + other parameters
sim_hours = 0.1
grid_size =[30, 30]
metabolite_diffusion = 5.e-6 # cm2 /s
metabolite_diffusion_h = metabolite_diffusion * 3600 # cm2 / h
toxin_diffusion = 5.e-7
space_width = 0.02
dt = 0.25 * (space_width * space_width) / metabolite_diffusion_h
max_cycles = int(sim_hours / dt)
founders_per_species = 10
initial_biomass = 1.e-10 # per founder, not per species in total
carbon_per_box = 4.e-7 # may wish to adjust this per carbon source to fix total carbon



rep = int(sys.argv[1])
carbon_sources = ["EX_csn_e", "EX_thr__L_e", "EX_tre_e"]
spatial_seed = math.floor(rep / len(carbon_sources))
carbon_source = carbon_sources[rep % len(carbon_sources)]
print(f"running spatial seed {spatial_seed} on carbon source {carbon_source}")

unlimited = ["EX_ca2_e","EX_cl_e","EX_co2_e","EX_cobalt2_e","EX_cu2_e",
             "EX_fe2_e","EX_fe3_e","EX_h_e","EX_h2o_e","EX_k_e","EX_mg2_e",
             "EX_mn2_e","EX_mobd_e","EX_na1_e","EX_nh4_e","EX_ni2_e",
             "EX_o2_e","EX_pi_e","EX_sel_e","EX_slnt_e","EX_so4_e",
             "EX_tungs_e","EX_zn2_e"]


### Load cobra models

model_shortnames = {"B" : "ST32113_grampos",
             "F" : "ST32123_grampos",
             "G" : "ST32124_grampos",
             "I" :"ST32133_grampos"}

model_ids = list(model_shortnames.values())
B = cobra.io.read_sbml_model(model_dir + model_shortnames["B"] + "_2022gapclosed.xml")
B.id = "B"
F = cobra.io.read_sbml_model(model_dir + model_shortnames["F"] + "_2022gapclosed.xml")
F.id = "F"
G = cobra.io.read_sbml_model(model_dir + model_shortnames["G"] + "_2022gapclosed.xml")
G.id = "G"
I = cobra.io.read_sbml_model(model_dir + model_shortnames["I"] + "_2022gapclosed.xml")
I.id = "I"

### add toxin production and exchanges to cobra models
### toxin production is always 1 mmol / gram cells born
from cobra import Metabolite, Reaction

B_toxin_c = Metabolite("B_toxin_c", compartment = "C_c")
B_toxin_e = Metabolite("B_toxin_e", compartment = "C_e")
B_objective = B.reactions.Growth
B_objective.add_metabolites({B_toxin_c: 1.})
B_toxin_tpp = Reaction("B_toxin_tpp", 
                       lower_bound = 0, upper_bound = 1000)
B_toxin_tpp.add_metabolites({B_toxin_c : -1,
                       B_toxin_e: 1})
B.add_reactions([B_toxin_tpp])
B.add_boundary(B_toxin_e, type = "exchange",
               lb = 0, ub = 1000)

F_toxin_c = Metabolite("F_toxin_c", compartment = "C_c")
F_toxin_e = Metabolite("F_toxin_e", compartment = "C_e")
F_objective = F.reactions.Growth
F_objective.add_metabolites({F_toxin_c: 1.})
F_toxin_tpp = Reaction("F_toxin_tpp", 
                       lower_bound = 0, upper_bound = 1000)
F_toxin_tpp.add_metabolites({F_toxin_c : -1,
                       F_toxin_e: 1})
F.add_reactions([F_toxin_tpp])
F.add_boundary(F_toxin_e, type = "exchange",
               lb = 0, ub = 1000)


G_toxin_c = Metabolite("G_toxin_c", compartment = "C_c")
G_toxin_e = Metabolite("G_toxin_e", compartment = "C_e")
G_objective = G.reactions.Growth
G_objective.add_metabolites({G_toxin_c: 1.})
G_toxin_tpp = Reaction("G_toxin_tpp", 
                       lower_bound = 0, upper_bound = 1000)
G_toxin_tpp.add_metabolites({G_toxin_c : -1,
                       G_toxin_e: 1})
G.add_reactions([G_toxin_tpp])
G.add_boundary(G_toxin_e, type = "exchange",
               lb = 0, ub = 1000)

B.add_boundary(F_toxin_e, type = "exchange",
               lb = 0, ub = 0)
B.add_boundary(G_toxin_e, type = "exchange",
               lb = 0, ub = 0)
F.add_boundary(B_toxin_e, type = "exchange",
               lb = 0, ub = 0)
F.add_boundary(G_toxin_e, type = "exchange",
               lb = 0, ub = 0)
G.add_boundary(B_toxin_e, type = "exchange",
               lb = 0, ub = 0)
G.add_boundary(F_toxin_e, type = "exchange",
               lb = 0, ub = 0)
I.add_boundary(B_toxin_e, type = "exchange",
               lb = 0, ub = 0)


#### Determine the minimum growth rate on the nutrient,
#### which will set the max growth rate for all species
cobra_medium = {key: 10. for key in unlimited}
cobra_medium
cobra_medium[carbon_source] = 10.

B.medium = {k:10. for k,v in cobra_medium.items() if k in [e.id for e in B.exchanges]}
B_gr = B.slim_optimize(0)
F.medium = {k:10. for k,v in cobra_medium.items() if k in [e.id for e in F.exchanges]}
F_gr = F.slim_optimize(0)
G.medium = {k:10. for k,v in cobra_medium.items() if k in [e.id for e in G.exchanges]}
G_gr = G.slim_optimize(0)
I.medium = {k:10. for k,v in cobra_medium.items() if k in [e.id for e in I.exchanges]}
I_gr = I.slim_optimize(0)
min_gr = min([B_gr, F_gr, G_gr, I_gr])

#### Make the comets models
B_comets = c.model(B)
B_comets.id = "B"
F_comets = c.model(F)
F_comets.id = "F"
G_comets = c.model(G)
G_comets.id = "G"
I_comets = c.model(I)
I_comets.id = "I"

B_comets.open_exchanges()
F_comets.open_exchanges()
G_comets.open_exchanges()
I_comets.open_exchanges()

#### Setup the toxin parameters
B_comets.add_multitoxin(2138, [229, 230], "ub", min_gr, [km_BF, km_BG], [hill_BF, hill_BG])
F_comets.add_multitoxin(2246, [265, 264], "ub", min_gr, [km_FB, km_FG], [hill_FB, hill_FG])
G_comets.add_multitoxin(2139, [234], "ub", min_gr, [km_GB], [hill_GB])
I_comets.add_multitoxin(2488, [275], "ub", min_gr, [km_IB], [hill_IB])

#### determine initial locations
def pick_unique_locations(width, height, n, edge_space = 0):
    locs = []
    while len(locs) < n:
        loc = (random.randrange(edge_space, width - edge_space),
               random.randrange(edge_space, height - edge_space))
        if loc not in locs:
            locs.append(loc)
    return(locs)

random.seed(spatial_seed)

total_founders = founders_per_species * 4
locs = pick_unique_locations(grid_size[0], grid_size[1], total_founders, 3)
B_locs = locs[0:founders_per_species]
F_locs = locs[founders_per_species:(founders_per_species*2)]
G_locs = locs[(founders_per_species*2):(founders_per_species*3)]
I_locs = locs[(founders_per_species*3):(founders_per_species*4)]

B_comets.initial_pop = [[loc[0], loc[1], initial_biomass] for loc in B_locs]
F_comets.initial_pop = [[loc[0], loc[1], initial_biomass] for loc in F_locs]
G_comets.initial_pop = [[loc[0], loc[1], initial_biomass] for loc in G_locs]
I_comets.initial_pop = [[loc[0], loc[1], initial_biomass] for loc in I_locs]

#### make the layout 
layout = c.layout([B_comets,F_comets,G_comets,I_comets])
layout.grid = grid_size
for met in unlimited:
    met = met[3:]
    if met in layout.media.metabolite.values:
        layout.set_specific_metabolite(met, 1000.)
layout.set_specific_metabolite(carbon_source[3:], carbon_per_box)

#### reduce the toxin diffusion constant to 5.e-7
layout.media.loc[layout.media.metabolite =="B_toxin_e", "diff_c"] = toxin_diffusion
layout.media.loc[layout.media.metabolite =="F_toxin_e", "diff_c"] = toxin_diffusion
layout.media.loc[layout.media.metabolite =="G_toxin_e", "diff_c"] = toxin_diffusion


#### make the parameters
params = c.params()
params.set_param("defaultDiffConst", metabolite_diffusion)
params.set_param("timeStep", dt) 
params.set_param("maxCycles", max_cycles)
params.set_param("spaceWidth", space_width) # cm, changing this will influence toxicity, need t play with
params.set_param("defaultKm", 0.000001)
params.set_param("maxSpaceBiomass", 1000)
params.set_param("minSpaceBiomass", 1.e-15)
params.set_param("writeMediaLog", False)
params.set_param("MediaLogRate", 1)
params.set_param("writeFluxLog", False)
params.set_param("FluxLogRate", 1)

#### run
sim = c.comets(layout, params)
sim.VERSION = "comets_multitoxin"
sim.set_classpath("bin", '/home/jeremy/Dropbox/work_related/harcombe_lab/segre/jars/comets_multitoxin.jar')
sim.run(delete_files = False)
sim.total_biomass.to_csv(f"total_biomass_rep_{rep}_carbon_{carbon_source}_seed_{spatial_seed}.csv")























