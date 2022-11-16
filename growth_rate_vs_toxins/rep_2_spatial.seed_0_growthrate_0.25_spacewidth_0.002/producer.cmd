SMATRIX  4  5
    1   1   -1
    1   2   -1
    2   2   1
    2   3   -1.02
    3   3   1.02
    3   5   1
    4   4   -1
    4   5   -1
//
BOUNDS 0 1000
    1   -0.25   1000.0
    4   -1000.0   1000.0
    5   -1000.0   0.0
//
OBJECTIVE
    3
//
METABOLITE_NAMES
    carbon_e
    carbon_c
    toxin_c
    toxin_e
//
REACTION_NAMES
    EX_carbon_e
    carbon_transfer
    Biomass
    EX_toxin_e
    toxin_transfer
//
EXCHANGE_REACTIONS
 1 4
//
OBJECTIVE_STYLE
MAX_OBJECTIVE_MIN_TOTAL
//
OPTIMIZER GUROBI
//
