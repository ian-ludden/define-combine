import pyomo.environ as pe
import pao

SOLVER_NAME = "pao.pyomo.MIBS"

def build_sample_model_one():
    M = pe.ConcreteModel()

    M.x = pe.Var(bounds=(2,6))
    M.y = pe.Var()

    M.L = pao.pyomo.SubModel(fixed=[M.x,M.y])
    M.L.z = pe.Var(bounds=(0,None), )

    M.o = pe.Objective(expr=M.x + 3*M.L.z, sense=pe.minimize)
    M.c = pe.Constraint(expr= M.x + M.y == 10)

    M.L.o = pe.Objective(expr=M.L.z, sense=pe.maximize)
    M.L.c1 = pe.Constraint(expr= M.x + M.L.z <= 8)
    M.L.c2 = pe.Constraint(expr= M.x + 4*M.L.z >= 8)
    M.L.c3 = pe.Constraint(expr= M.x + 2*M.L.z <= 13)
    
    return M


def build_sample_model_two():
    M = pe.ConcreteModel()
    M.x = pe.Var()#bounds=(0,None))#domain=pe.NonNegativeIntegers)

    M.L = pao.pyomo.SubModel(fixed=[M.x])
    M.L.y = pe.Var()#bounds=(0,None))#domain=pe.NonNegativeIntegers)

    M.o = pe.Objective(expr=M.x + 10*M.L.y, sense=pe.maximize)
    
    M.L.o = pe.Objective(expr=-1*M.L.y, sense=pe.maximize)
    M.L.c1 = pe.Constraint(expr= -25*M.x + 20*M.L.y <= 30)
    M.L.c2 = pe.Constraint(expr= 1*M.x + 2*M.L.y <= 10)
    M.L.c3 = pe.Constraint(expr= 2*M.x - 1*M.L.y <= 15)
    M.L.c4 = pe.Constraint(expr= 2*M.x + 10*M.L.y >= 15)

    return M


def build_sample_model_three():
    M = pe.ConcreteModel()
    M.x = pe.Var(bounds=(0,None))
    M.y = pe.Var(bounds=(0,None))

    M.o = pe.Objective(expr=M.x - 4*M.y)

    M.L = pao.pyomo.SubModel(fixed=M.x)
    M.L.o = pe.Objective(expr=M.y)
    M.L.c1 = pe.Constraint(expr=   -M.x -   M.y <= -3)
    M.L.c2 = pe.Constraint(expr= -2*M.x +   M.y <=  0)
    M.L.c3 = pe.Constraint(expr=  2*M.x +   M.y <= 12)
    M.L.c4 = pe.Constraint(expr=  3*M.x - 2*M.y <=  4)

    return M
    

if __name__ == '__main__':
    # M = build_sample_model_one()
    # opt = pao.Solver(SOLVER_NAME)
    # results = opt.solve(M)
    # print(results)
    # print(M.x.value, M.y.value, M.L.z.value)

    # M2 = build_sample_model_two()
    # opt = pao.Solver(SOLVER_NAME)
    # results2 = opt.solve(M2)
    # print(results2)
    # print(M2.x.value, M2.L.y.value)

    M = pe.ConcreteModel("Name of model.")
    M.x = pe.Var(bounds=(0,None))
    M.y = pe.Var(bounds=(0,None))

    M.o = pe.Objective(expr=M.x - 4*M.y)

    M.L = pao.pyomo.SubModel(fixed=M.x)
    M.L.o = pe.Objective(expr=M.y)
    M.L.c1 = pe.Constraint(expr=   -M.x -   M.y <= -3)
    M.L.c2 = pe.Constraint(expr= -2*M.x +   M.y <=  0)
    M.L.c3 = pe.Constraint(expr=  2*M.x +   M.y <= 12)
    M.L.c4 = pe.Constraint(expr=  3*M.x - 2*M.y <=  4)

    print(M)
    for v in [M.x, M.y]:
        print(v.value)

    M.pprint()
    
    
    with pao.Solver(SOLVER_NAME) as solver:
        print(type(solver), solver)
        results3 = solver.solve(M)
        print(results3.solver.termination_condition)
        print(results3)
        print("Optimal termination?", results3.check_optimal_termination())
        print("x = ", M.x.value, ", y = ", M.y.value, sep='')

