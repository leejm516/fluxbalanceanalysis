# MILP solve

using BenchmarkTools
using MathOptInterface
using CPLEX
using MAT

const MOI = MathOptInterface

vars = matread("./MILPProblem_ijr.mat")
milpPr = vars["MILPpr"]

S = milpPr["A"];
b = vec(milpPr["b"]);
lb = vec(milpPr["lb"]);
ub = vec(milpPr["ub"]);
c = vec(milpPr["c"]);
csense = vec(milpPr["csense"]);
vartype = vec(milpPr["vartype"]);

num_constraints, num_variables = size(S);

@btime begin
    model = CPLEX.Optimizer()

    # Solver parameters
    cpxParams = Dict("CPX_PARAM_SCRIND" => 0,
                    "CPXPARAM_Simplex_Tolerances_Feasibility" => 1e-9,
                    "CPXPARAM_Simplex_Tolerances_Optimality" => 1e-9,
                    "CPXPARAM_Network_Tolerances_Optimality" => 1e-9,
                    "CPXPARAM_Network_Tolerances_Feasibility" => 1e-9,
                    "CPXPARAM_Barrier_ConvergeTol" => 1e-9,
                    "CPXPARAM_MIP_Tolerances_MIPGap" => 1e-12,
                    "CPXPARAM_MIP_Tolerances_AbsMIPGap" => 1e-12,
                    "CPXPARAM_MIP_Tolerances_Integrality" => 1e-12)
    model.params = cpxParams
    MOI.empty!(model)


    # Create the variables in the problem.
    x = MOI.add_variables(model, num_variables)

    # lb, ub, and Binary variables
    for i = 1:num_variables
        if vartype[i] == "B"
            MOI.add_constraint(model, MOI.SingleVariable(x[i]), MOI.ZeroOne())
        else
            MOI.add_constraint(model, MOI.SingleVariable(x[i]), MOI.GreaterThan(lb[i]));
            MOI.add_constraint(model, MOI.SingleVariable(x[i]), MOI.LessThan(ub[i]));
        end
    end


    # Objective function
    objective_function = MOI.ScalarAffineFunction(MOI.ScalarAffineTerm.(c, x), 0.0)
    MOI.set(model, MOI.ObjectiveFunction{MOI.ScalarAffineFunction{Float64}}(),
            objective_function)
    MOI.set(model, MOI.ObjectiveSense(), MOI.MAX_SENSE)


    beq_idx = findall(x -> x == "E", csense)
    beq_len = length(beq_idx)

    bl_idx = findall(x -> x == "L", csense)
    bl_len = length(bl_idx)

    bg_idx = findall(x -> x == "G", csense)
    bg_len = length(bg_idx)

    # Mass balance
    massbal_function = MOI.VectorAffineTerm.(1:beq_len, MOI.ScalarAffineTerm.(S[beq_idx, :],
                        reshape(x, 1, num_variables)))
    massbal_cons = MOI.VectorAffineFunction(vec(massbal_function), zeros(beq_len))
    MOI.add_constraint(model, massbal_cons, MOI.Zeros(beq_len))

    # loopless
    bl_function = MOI.VectorAffineTerm.(1:bl_len, MOI.ScalarAffineTerm.(S[bl_idx, :],
                    reshape(x, 1, num_variables)))
    bl_cons = MOI.VectorAffineFunction(vec(bl_function), (-1.0)*b[bl_idx])
    MOI.add_constraint(model, bl_cons, MOI.Nonpositives(bl_len))

    bg_function = MOI.VectorAffineTerm.(1:bg_len, MOI.ScalarAffineTerm.(S[bg_idx, :],
                    reshape(x, 1, num_variables)))
    bg_cons = MOI.VectorAffineFunction(vec(bg_function), (-1.0)*b[bg_idx])
    MOI.add_constraint(model, bg_cons, MOI.Nonnegatives(bg_len))
end


# All set!
MOI.optimize!(model)
