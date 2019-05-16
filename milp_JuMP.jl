# MILP solve

using BenchmarkTools
using JuMP
using CPLEX
using MAT

vars = matread("./MILPProblem_ijr.mat")
milpPr = vars["MILPpr"]

S = milpPr["A"];
b = vec(milpPr["b"]);
lb = vec(milpPr["lb"]);
ub = vec(milpPr["ub"]);
c = vec(milpPr["c"]);
csense = vec(milpPr["csense"]);
vartype = vec(milpPr["vartype"]);

m, n = size(S);

@btime begin
    MOImodel = CPLEX.Optimizer()
    cpxParams = Dict("CPX_PARAM_SCRIND" => 0,
                    "CPXPARAM_Simplex_Tolerances_Feasibility" => 1e-9,
                    "CPXPARAM_Simplex_Tolerances_Optimality" => 1e-9,
                    "CPXPARAM_Network_Tolerances_Optimality" => 1e-9,
                    "CPXPARAM_Network_Tolerances_Feasibility" => 1e-9,
                    "CPXPARAM_Barrier_ConvergeTol" => 1e-9,
                    "CPXPARAM_MIP_Tolerances_MIPGap" => 1e-12,
                    "CPXPARAM_MIP_Tolerances_AbsMIPGap" => 1e-12,
                    "CPXPARAM_MIP_Tolerances_Integrality" => 1e-12)
    MOImodel.params = cpxParams
    MOI.empty!(MOImodel)

    model = direct_model(MOImodel)


    bin_idx = findall(x -> x =="B", vartype)
    con_idx = findall(x -> x =="C", vartype)
    flux_idx, g_idx = con_idx[1:(bin_idx[1]-1)], con_idx[bin_idx[1]:end]

    @variable(model, lb[flux_idx[i]] <= v[i=1:length(flux_idx)] <= ub[flux_idx[i]]) # Fluxes
    @variable(model, lc[i =1:length(bin_idx)], Bin) # a variables
    @variable(model, lb[g_idx[i]] <= lg[i=1:length(g_idx)] <= ub[g_idx[i]]) # G variables

    beq_idx = findall(x -> x == "E", csense)
    bg_idx = findall(x -> x == "G", csense)
    bl_idx = findall(x -> x == "L", csense)

    @constraint(model, S[beq_idx, :]*vcat(v, lc, lg) .== 0) # Mass balance
    @constraint(model, S[bl_idx, :]*vcat(v, lc, lg) - b[bl_idx] .<= 0) # loopless
    @constraint(model, S[bg_idx, :]*vcat(v, lc, lg) - b[bg_idx] .>= 0) # loopless

    # Set the objective function.
    if milpPr["osense"] < 0.0
        @objective(model, Max, v[findfirst(x->x !=0, c)])
    else
        @objective(model, Min, v[findfirst(x->x !=0, c)])
    end
end

    optimize!(model)
