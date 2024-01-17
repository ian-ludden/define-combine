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
    M.x = pe.Var(domain=pe.NonNegativeIntegers)

    M.L = pao.pyomo.SubModel(fixed=[M.x])
    M.L.y = pe.Var(domain=pe.NonNegativeIntegers)

    M.o = pe.Objective(expr=M.x + 10*M.L.y, sense=pe.maximize)
    
    M.L.o = pe.Objective(expr=-1*M.L.y, sense=pe.maximize)
    M.L.c1 = pe.Constraint(expr= -25*M.x + 20*M.L.y <= 30)
    M.L.c2 = pe.Constraint(expr= 1*M.x + 2*M.L.y <= 10)
    M.L.c3 = pe.Constraint(expr= 2*M.x - 1*M.L.y <= 15)
    M.L.c4 = pe.Constraint(expr= 2*M.x + 10*M.L.y >= 15)

    return M


if __name__ == '__main__':
    M = build_sample_model_one()
    opt = pao.Solver(SOLVER_NAME)
    # results = opt.solve(M)
    # print(M.x.value, M.y.value, M.L.z.value)

    M = build_sample_model_two()
    opt = pao.Solver(SOLVER_NAME)
    # results = opt.solve(M)
    # print(M.x.value, M.L.y.value)

