SMATRIX  3  4
    1   1   -1
    1   2   -1
    2   2   1
    2   3   -1.0
    3   4   -1
//
BOUNDS 0 1000
    1   -0.25   1000.0
    4   -1000.0   1000.0
//
OBJECTIVE
    3
//
METABOLITE_NAMES
    carbon_e
    carbon_c
    toxin_e
//
REACTION_NAMES
    EX_carbon_e
    carbon_transfer
    Biomass
    EX_toxin_e
//
EXCHANGE_REACTIONS
 1 4
//
MET_REACTION_SIGNAL
3 2 ub bounded_linear 0.25 0.0 -0.25 1.0
//
OBJECTIVE_STYLE
MAX_OBJECTIVE_MIN_TOTAL
//
OPTIMIZER GUROBI
//
