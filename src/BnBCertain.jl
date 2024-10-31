module BnBCertain
    using Revise, LinearAlgebra
    using DAQP, ASCertain, CertificationTools
    using JuMP, Gurobi, CDDLib, GLPK          
    using Polyhedra, Plots, PGFPlotsX, LaTeXStrings
    #using StatsBase, DataFrames, CSV, JLD2 #Debugger, 

    export MPLIQP, MIQP, generate_mpMIQP, 
           BNBSettings, BNBRegion,
           param_sampling, param_location, compute_centers,
           bnb_cert, #bnb_cert_hu, 
           compare_on_cert    

    @enum BinState LOWER UPPER FREE 
    @enum BNBState Root InfeasCut DomCut IntFeasCut Branch PartialDomCut Unbounded
    @enum SortStrategy DepthFirst BreadthFirst BestFirst #DF BrF BF #
    @enum BranchStrategy FirstBin MostInfeasBin MinBin MaxBin LastBin #branching methods: in order, randomely, most infeasible branching ...
           
    include("types.jl")
    include("utilsMI.jl")
    include("bnb_cert.jl")  
    #include("bnb_cert_hu.jl") 
    include("bnb_func_cert.jl")
    include("bnb_on_hu.jl") #
    include("compare_on_cert.jl")
    include("find_solutionQLP.jl")
    include("MIP_online.jl")
end