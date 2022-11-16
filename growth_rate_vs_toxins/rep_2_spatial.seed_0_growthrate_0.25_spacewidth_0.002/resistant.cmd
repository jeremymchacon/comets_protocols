SMATRIX  2  3
    1   1   -1
    1   2   -1
    2   2   1
    2   3   -1.0
//
BOUNDS 0 1000
    1   -0.25   1000.0
//
OBJECTIVE
    3
//
METABOLITE_NAMES
    carbon_e
    carbon_c
//
REACTION_NAMES
    EX_carbon_e
    carbon_transfer
    Biomass
//
EXCHANGE_REACTIONS
 1
//
OBJECTIVE_STYLE
MAX_OBJECTIVE_MIN_TOTAL
//
OPTIMIZER GUROBI
//
