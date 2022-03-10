import gurobipy as gp
import random

from utils.mip import build_model, build_model_continuous, summarize_model_results

# TODO: define Node class for branch and bound solver. 
#   index/id, 
#   forced_0/1_x/y variable names, 
#   F_hpr, 
#   etc.
class Node():
    pass


if __name__ == '__main__':
    random.seed(2021) # For reproducibility

    NUM_ROWS = 4
    NUM_COLS = 4
    NUM_DISTRICTS = 4
    NUM_SUBDISTRICTS = 2 * NUM_DISTRICTS

    NUM_UNITS = NUM_ROWS * NUM_COLS
    UNITS_PER_DISTRICT = NUM_UNITS // NUM_DISTRICTS
    UNITS_PER_SUBDISTRICT = NUM_UNITS // NUM_SUBDISTRICTS

    NUM_DEFINER_SUPPORTERS = int(NUM_UNITS * 0.50)
    NUM_COMBINER_SUPPORTERS = NUM_UNITS - NUM_DEFINER_SUPPORTERS

    assert(NUM_UNITS // NUM_SUBDISTRICTS == NUM_UNITS / NUM_SUBDISTRICTS) # Require exact population balance

    print("Definer has", NUM_DEFINER_SUPPORTERS, "supporters.")

    unit_weights = [-1] * NUM_COMBINER_SUPPORTERS + [1] * NUM_DEFINER_SUPPORTERS
    random.shuffle(unit_weights)
    unit_ids = [i for i in range(NUM_UNITS)]

    # Print voter distribution
    index = 0
    while index < NUM_UNITS:
        if index % NUM_COLS == 0:
            print()
        
        weight = unit_weights[index]
        output_char = 1 if weight == 1 else 0 # Use 0 to represent -1 since it's a single character
        print(output_char, end='')

        index += 1
    print()
    
    try:
        # Build model
        m = build_model(num_rows=NUM_ROWS, num_cols=NUM_COLS, \
            num_districts=NUM_DISTRICTS, unit_ids=unit_ids, unit_weights=unit_weights)

        # Optimize model
        m.optimize()
        summarize_model_results(m, unit_ids=unit_ids, unit_weights=unit_weights, num_cols=NUM_COLS)

        print("\n*** RESET and RESOLVE ***\n")
        m.reset()
        m.addConstr(m.getVarByName("x[0,4]") == 0, "forced_zero_x_{}_{}".format(0, 4))
        
        m.optimize()
        summarize_model_results(m, unit_ids=unit_ids, unit_weights=unit_weights, num_cols=NUM_COLS)

        print("\n*** SOLVE CONTINUOUS RELAXATION ***\n")
        m_cont = build_model_continuous(num_rows=NUM_ROWS, num_cols=NUM_COLS, \
            num_districts=NUM_DISTRICTS, unit_ids=unit_ids, unit_weights=unit_weights)
        
        m_cont.optimize()
        summarize_model_results(m_cont, unit_ids=unit_ids, unit_weights=unit_weights, num_cols=NUM_COLS)

    except gp.GurobiError as e:
        print('Error code ' + str(e.errno) + ': ' + str(e))

    except AttributeError as e:
        print('Encountered an attribute error.')
