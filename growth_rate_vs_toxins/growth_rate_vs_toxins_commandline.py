import cometspy as c
import cobra
import math
import random
from cobra import Metabolite, Reaction, Model
import pandas as pd
import os
import sys

# note: this will change for your system, and may go away entirely on MSI
os.environ['GUROBI_COMETS_HOME'] = '/opt/gurobi900/linux64'


## Here is some other relevant treatment information--don't change these values without making sure it is noted in the output somehow
base_hours = 1 # for growth rate = 1, this will increase for slower growth rates (see below)
resistance_cost = 0. # 
production_cost = 0.02 # fraction of growth
toxin_coefficient = 1. # mmol toxin made per gram cells born
initial_biomass = 1.e-10
grid_size =[30, 30] # this was 100x100 in the actual sims!
metabolite_diff = 5.e-6
toxin_diff = 5.e-7

## figure out treatment information, and print that to stdout
replicate = int(sys.argv[1])

growth_rates = [0.125, 0.25, 0.5, 1.]
space_widths = [0.002, 0.004]


# iterate through these to figure out which to use, and the spatial seed
total_treatments = len(growth_rates) * len(space_widths)
spatial_seed = math.floor(replicate / total_treatments)
trt_num = 1 + replicate % total_treatments
print(f"treatment # = {trt_num}")
counter = 0
for gr in growth_rates:
    growth_rate = gr
    for sw in space_widths:
        counter += 1
        space_width = sw
        if counter == trt_num:
            break
    if counter == trt_num:
        break


print(f"rep = {replicate}, spatial_seed = {spatial_seed}, growthrate = {growth_rate}, space_width = {space_width}")

## I like to save results in their own directory
base = f"/home/jeremy/Dropbox/work_related/harcombe_lab/comets_protocols/growth_rate_vs_toxins/rep_{replicate}_spatial.seed_{spatial_seed}_growthrate_{growth_rate}_spacewidth_{space_width}/"
print(f"base directory : {base}")

try:
    os.mkdir(base)
except:
    print("Warning: base directory exists. data may be overwritten")

## I like to check if a simulation has already been run in this directory, and if so, not run another one
tb_file = base + "total.biomass.csv"
founder_file = base + "founder_locs.csv"

if os.path.isfile(tb_file):
    print("Total Biomass file exists, ending early")
    import sys
    sys.exit()
    
# If there is no total biomass file there, then do the simulation

random.seed(spatial_seed)    


# time step needs to shrink as D increases or dx decreases:
dt = 0.25 * ((space_width)**2) / metabolite_diff
# sim hours aka max cycles is hard to pre-determine, you may need to adjust later. 
hours = base_hours / growth_rate
max_cycles = int(hours / dt)
print(f"dt = {dt}, max_cycles = {max_cycles}")


# make params
p = c.params()
p.set_param("defaultVmax", 1.)
p.set_param("defaultKm", 0.000001) 
p.set_param('maxCycles', max_cycles)
p.set_param('timeStep', dt)
p.set_param('spaceWidth', space_width)
p.set_param('writeMediaLog', False)
p.set_param('writeBiomassLog', False)
p.set_param('minSpaceBiomass', 1.e-15)
p.set_param('defaultDiffConst', metabolite_diff)
p.set_param('flowDiffRate', 3.e-9)
p.set_param('writeFluxLog', False)
p.set_param("totalBiomassLogRate", 1)
p.set_param("numDiffPerStep", 1) # this SHOULD be fine since we determined dt based on the faster diff rate


def pick_unique_locations(width, height, n, edge_space = 0):
    locs = []
    while len(locs) < n:
        loc = (random.randrange(edge_space, width - edge_space),
               random.randrange(edge_space, height - edge_space))
        if loc not in locs:
            locs.append(loc)
    return(locs)

def make_RPS_cobra_models_biomass_cost(growth_rate = 1., toxin_cost = 0.02, resistance_cost = 0.0, toxin_prod = 1.):
    # this is different from make RPS cobra models because the toxin and resistance cost is put into the
    # biomass reaction, rather than a separate reaction. this is important because toxins will only
    # be created during growth.also, no need to multiply costs by growth rate. in all cases, they
    # vary directly with growth rate. finally, in all cases as much toxin is produced as carbon is taken
    # up, regardless of growth rate.
    #
    # These toy models use "carbon_c" to grow--that is the only resource
    #
    # the toxin is called "toxin_e"
    '''
    @toxin_prod:  a multiplier.  the multiple of mmol of toxin made per gram of carbon taken up
    '''
    carbon_e = Metabolite(id = "carbon_e",
               compartment = "e")
    carbon_c = Metabolite(id = "carbon_c",
                   compartment = "c")
    EX_carbon_e = Reaction(id = "EX_carbon_e",
                      lower_bound = -growth_rate, # growth rate
                      upper_bound = 1000.)
    EX_carbon_e.add_metabolites({carbon_e: -1})
    carbon_transfer = Reaction(id = "carbon_transfer",
                          lower_bound = 0.,
                          upper_bound = 1000.)
    carbon_transfer.add_metabolites({carbon_e: -1,
                                carbon_c: 1})
    Biomass = Reaction(id = "Biomass",
                  lower_bound = 0.,
                  upper_bound = 1000.)
    Biomass.add_metabolites({carbon_c: -1.})
    # make the toxicity-related metabolites and reactions
    toxin_c = Metabolite(id = "toxin_c", compartment = "c")
    toxin_e = Metabolite(id = "toxin_e", compartment = "e")

    EX_toxin_e = Reaction(id = "EX_toxin_e",
                         lower_bound = -1000.,
                         upper_bound = 1000.)
    EX_toxin_e.add_metabolites({toxin_e: -1})

    toxin_transfer = Reaction(id = "toxin_transfer",
                             lower_bound = -1000.,
                             upper_bound = 0.)
    toxin_transfer.add_metabolites({toxin_e: -1,
                                   toxin_c: 1})

    Biomass_producer = Reaction(id = "Biomass",
                  lower_bound = 0.,
                  upper_bound = 1000.)
    Biomass_producer.add_metabolites({carbon_c: -(1. + toxin_cost + resistance_cost),
                                     toxin_c: toxin_prod * (1. + toxin_cost + resistance_cost)})

    Biomass_resistant = Reaction(id = "Biomass",
                  lower_bound = 0.,
                  upper_bound = 1000.)
    Biomass_resistant.add_metabolites({carbon_c: -(1. + resistance_cost)})

    producer = Model("producer")
    producer.add_reactions([EX_carbon_e, carbon_transfer, Biomass_producer,
                        EX_toxin_e, toxin_transfer])
    producer.objective = Biomass_producer

    resistant = Model("resistant")
    resistant.add_reactions([EX_carbon_e, carbon_transfer, Biomass_resistant])
    resistant.objective = Biomass_resistant

    susceptible = Model("susceptible")
    susceptible.add_reactions([EX_carbon_e, carbon_transfer, Biomass,
                        EX_toxin_e])
    susceptible.objective = Biomass
    return((producer, resistant, susceptible))


producer, resistant, susceptible = make_RPS_cobra_models_biomass_cost(growth_rate = growth_rate, 
                                                                      toxin_cost = production_cost,
                                                                      resistance_cost = resistance_cost,
                                                                      toxin_prod = toxin_coefficient)
                                                                      
P = c.model(producer)
R = c.model(resistant)
S = c.model(susceptible)

P.obj_style = "MAX_OBJECTIVE_MIN_TOTAL"
R.obj_style = "MAX_OBJECTIVE_MIN_TOTAL"
S.obj_style = "MAX_OBJECTIVE_MIN_TOTAL"


### Fix COMETS thinking biomass reactions with only reactants are exchanges
R.reactions.loc[R.reactions.REACTION_NAMES == "Biomass", "EXCH"] = False
index = R.reactions.loc[R.reactions.REACTION_NAMES == "Biomass", "EXCH_IND"].values[0]
R.reactions.loc[R.reactions.REACTION_NAMES == "Biomass", "EXCH_IND"] = 0
R.reactions.loc[R.reactions.EXCH_IND > index, "EXCH_IND"] -= 1

S.reactions.loc[S.reactions.REACTION_NAMES == "Biomass", "EXCH"] = False
index = S.reactions.loc[S.reactions.REACTION_NAMES == "Biomass", "EXCH_IND"].values[0]
S.reactions.loc[S.reactions.REACTION_NAMES == "Biomass", "EXCH_IND"] = 0
S.reactions.loc[S.reactions.EXCH_IND > index, "EXCH_IND"] -= 1

### setup the toxicity
biomass_id = S.reactions.loc[S.reactions.REACTION_NAMES == "Biomass", "ID"].values[0]
toxin_exch_id = S.reactions.loc[S.reactions.REACTION_NAMES == "EX_toxin_e", "EXCH_IND"].values[0]
gr_absense_of_toxin = growth_rate
toxin_conc_where_effect_starts = 0.
slope_of_toxin_effect = -growth_rate
toxin_conc_where_effect_saturates = 1.
S.add_signal(biomass_id, toxin_exch_id, 'ub', 'bounded_linear', 
             parms = [gr_absense_of_toxin, 
                      toxin_conc_where_effect_starts, 
                      slope_of_toxin_effect,
                     toxin_conc_where_effect_saturates])


# HERE IS WHERE THE NUMBER OF FOUNDERS IS SET
locs = pick_unique_locations(grid_size[0], grid_size[1], 50, 3)
producer_locs = locs[0:5]
resistant_locs = locs[5:10]
susceptible_locs = locs[10:]

S.initial_pop = [[loc[0], loc[1], initial_biomass] for loc in susceptible_locs]
P.initial_pop = [[loc[0], loc[1], initial_biomass] for loc in producer_locs]
R.initial_pop = [[loc[0], loc[1], initial_biomass] for loc in resistant_locs]


## Make the layout
l = c.layout([P, R, S])
l.grid = grid_size
l.set_specific_metabolite("carbon_e", 4.e-7)
l.media.loc[l.media.metabolite == "toxin_e", "diff_c"] = toxin_diff


# Saving initial founder locs
spatial_data = pd.DataFrame(columns = ["strain", "x", "y"])
producer_data = pd.DataFrame({"strain" : "producer",
                             "x" : [x[0] for x in producer_locs],
                             "y" : [x[1] for x in producer_locs]})
resistant_data = pd.DataFrame({"strain" : "resistant",
                             "x" : [x[0] for x in resistant_locs],
                             "y" : [x[1] for x in resistant_locs]})
susceptible_data = pd.DataFrame({"strain" : "susceptible",
                             "x" : [x[0] for x in susceptible_locs],
                             "y" : [x[1] for x in susceptible_locs]})
spatial_data = spatial_data.append(producer_data, 
                                   ignore_index = True).append(resistant_data, 
                                                               ignore_index = True).append(susceptible_data, ignore_index = True)
spatial_data.to_csv(founder_file)


## run the model, then save total biomass

sim = c.comets(l, p)
sim.working_dir = base
sim.VERSION = 'comets_multitoxin.jar'
sim.set_classpath('bin', '/home/jeremy/Dropbox/work_related/harcombe_lab/segre/jars/comets_multitoxin.jar')
try:
    sim.run(delete_files = False)
except:
    print(sim.run_output)
sim.total_biomass.to_csv(tb_file)