# multi-parametric mixed-integer linear/quadratic program
mutable struct MPMIQP
  H::Matrix{Float64}
  f::Vector{Float64}
  f_theta::Matrix{Float64}
  H_theta::Matrix{Float64}
  A::Matrix{Float64}
  b::Vector{Float64}
  W::Matrix{Float64} 
  bin_ind::Vector{Int64}
  bounds_table::Vector{Int64}
  senses::Vector{Cint} 

  MPMIQP()=new()
  MPMIQP(H,f,f_theta,H_theta,A,b,W,bin_ind,bounds_table=Int64[],senses=Cint[])= new(H,f,f_theta,H_theta,A,b,W,bin_ind,bounds_table,senses)
end

mutable struct MIQP
    H::Matrix{Float64}
    f::Vector{Float64}
    A::Matrix{Float64}
    b::Vector{Float64}
    bin_ind::Vector{Int64}
    bounds_table::Vector{Int64}
    senses::Vector{Cint} 

    MIQP()=new()
    MIQP(H,f,A,b,bin_ind,bounds_table,senses) =new(H,f,A,b,bin_ind,bounds_table,senses)   
end
  
function MIQP(mpMIQP::MPMIQP, theta::Vector{Float64})
 return MIQP(mpMIQP.H,mpMIQP.f+mpMIQP.f_theta*theta,mpMIQP.A,mpMIQP.b+mpMIQP.W*theta,mpMIQP.bin_ind,mpMIQP.bounds_table,mpMIQP.senses) 
end

# parametric BNB node with relaxation
mutable struct NodeMP
    # relaxation data
    H::Matrix{Float64}
    f::Vector{Float64}
    f_theta::Matrix{Float64}
    H_theta::Matrix{Float64} 
    A::Matrix{Float64}
    b::Vector{Float64}
    W::Matrix{Float64} 
    bounds_table::Vector{Int64}
    senses::Vector{Cint} #Vector{ConstraintType}
    bin_ind::Vector{Int64}
    #bnb data
    bin_values::Vector{Int64}
    bin_indices::Matrix{Int64}
    bin_states::Vector{BinState} 
    n::UInt16
    nb::UInt16
    bin_ind_all::Vector{Int64}
    AS::Vector{Int64} 
    priority::Float64 
    #branchrule::BranchStrategy
    level::Int32
    extra_cost::Float64 
    extra_cost_param::Matrix{Float64}
    priority_BF::Matrix{Float64}
    feasible::Bool #when warm-starting 
    

    NodeMP()=new()
    NodeMP(H,f,f_theta,H_theta,A,b,W,bounds_table,senses,bin_ind,bin_values,bin_indices,bin_states,n,nb,bin_ind_all,AS,priority,level,extra_cost,extra_cost_param,priority_BF,feasible)= new(H,f,f_theta,H_theta,A,b,W,bounds_table,senses,bin_ind,bin_values,bin_indices,bin_states,n,nb,bin_ind_all,AS,priority,level,extra_cost,extra_cost_param,priority_BF,feasible)
  end

# node in online BNB
mutable struct Node
    H::Matrix{Float64}
    f::Vector{Float64}
    A::Matrix{Float64}
    b::Vector{Float64} 
    bounds_table::Vector{Int64}
    senses::Vector{Cint} #Vector{ConstraintType}    
    bin_ind::Vector{Int64}

    bin_values::Vector{Int64}
    bin_indices::Matrix{Int64}
    bin_states::Vector{BinState} 
    n::UInt16
    nb::UInt16
    bin_ind_all::Vector{Int64}
    AS::Vector{Int64} 
    priority::Float64 
    #branchrule::BranchStrategy
    level::Int32
    extra_cost::Float64
    feasible::Bool  #when warm-starting 
  
    Node()=new()
    Node(H,f,A,b,bounds_table,senses,bin_ind,bin_values,bin_indices,bin_states,n,nb,bin_ind_all,AS,priority,level,extra_cost,feasible) =new(H,f,A,b,bounds_table,senses,bin_ind,bin_values,bin_indices,bin_states,n,nb,bin_ind_all,AS,priority,level,extra_cost,feasible)    
end

abstract type AbstractRegion end
# parametric solution (F th + G) and value function (th A th + B th + C)
mutable struct ParamSol    
    F_theta::Matrix{Float64}
    G_theta::Vector{Float64}
    A_theta::Matrix{Float64}
    B_theta::Matrix{Float64} 
    C_theta::Float64
end

mutable struct ComplexityMeasure    
  iter::Int64
  #iter_avg::Float64
  iter_hu::Int64
  node::Int64
  decomp_node::Int64
  tot_gen_nodes::Int64
  loop::Int64
  flops::Float64
  time::Float64
  iter_tot::Int64
  iter_hu_tot::Int64
end

mutable struct BNBRegion <:AbstractRegion
    Ath::Matrix{Float64}
    bth::Vector{Float64} 
    iter::Int64
    #Jbar::Vector{Matrix{Float64}}  #xbar::Vector{Float64}
    AS::Vector{Int64}
    sol_bar::Vector{ParamSol} 
    T::Vector{NodeMP}
    #kappa::Vector{Int64} #kappa::Dict{Any,Any}
    kappa::ComplexityMeasure
    seq_nodes::Matrix{Int64}
    #seq_status::Vector{BNBState}
    tree_status::Matrix{Any}
    tot_nodes::Int64 
    #decomp_nodes::Int32 
    #tot_gen_nodes::Int32
    #iter_hu::Int32
    #iter_avg::Float32
    #tree_avg::Float32
    #ind::Vector{Int64}
    info_status::Matrix{Any} #[k j level bnbstate binaryvalues] 
end


Base.@kwdef mutable struct BNBSettings
  sorting_method::SortStrategy = DepthFirst # "depth_first", "breadth_first", "best_first"
  branching_rule::BranchStrategy = FirstBin 
  bound_tot_nodes::Int64 = 1000;   # limit on total nodes in BnB
  chosen_solver::String = "GLPK" #"Gurobi"
  daqp_solver_method::Bool = true;
  warm_start::Bool = false;
  cut_int_feas::Bool = true;
  problem_type_QP::Bool = false; #QP: true; LP:false
  LPcompare_call_polyhedra::Bool = true; #false; # call polyhedron package to compare LPs
  compute_minimal_rep::Bool = true; #false; # find minimal reperesntation of polyhedrons
  remove_redundancy::Bool = true; #Use certification tool to remove redundancy
  remove_empty_reg::Bool = true; #remove empty polyhedrons after calling certifictaion code
  store_final_regions::Bool = false;
  online_gurobi::Bool = false; # find optimal online solution with gurobi
  check_dom_cut::Bool = true; #not ignore dominance cut 
  check_int_feas_cut::Bool = true; #not ignore integer-feasible cut 
  add_extra_cost::Bool = true; #to J at each step
  MInfB_check_chebyC::Bool  = true; #used in most infeas approach to calculate distace from chebycenter when all are binary  
  #____________heuristic____________
  apply_heuristic::Bool = false;
  hu_trigger::Vector{Int64} = [1]; #trigger heuristic function, defult: just after the root node
  add_iter_hu::Bool = false;    # Add up iteration from heuristic to total iteration number
  bound_loop_hu::Int16 = 50;    # iteration loop limit
  bound_iter_hu::Int16 = 500;   # simplex iteration number limit
  bound_time_hu::Int16 = 5000;  # time loop limit
  arith_mean::Bool = false;     # take arithmatic mean or geometric mean, f.e. from online time
  #______heuristic type_________
  heuris_start::Bool = false;
  heuris_start_type::String = "FP" #"RENS" # "Div"
  heuris_imp::Bool = false;
  heuris_imp_type::String = "LB" #"RINS" # ""
  #______Feasibility Pump
  #heuristic_FP::Bool = false;   # Feasibility Pump heuristics
  heuristic_obj_FP::Bool = true;# Objective Feasibility Pump heuristics
  obj_FP_alpha::Float32 = 1;    # Objective Feasibility Pump heuristics alpha: scale between delta and c
  obj_FP_phi::Float32 = 0.9;    # Objective Feasibility Pump heuristics:scale alpha
  nbin_flip_FP::Float16 = .5;   # ratio of number of binaries to be flipped in case of stalling
  FP_cycle_flip::Bool = true;   # Flip variables if cycle detected in OFP
  FP_cycl_abort::Bool = true;   # Abort FP algorithm if cycle detected in OFP
  FP_cycl_perturb::Vector{Float64} = zeros(8); #if not abort, purturb the solution
  #______Diving
  #heuristic_Div::Bool = false; # Diving heuristics
  heuris_Div_br0::Bool = true; # Diving heuristics:branch on xi = 0, false: branch on 1
  #______Rounding: Rens
  #heuristic_RENS::Bool = false; # RENS heuristics
  heuris_rens_ratio::Float16 = 0.5; #ratio of fixed relaxed binary variables
  #______Local branching
  #heuristic_LB::Bool = false; # Local branching heuristics
  heuris_start_LB::Bool = false; # Call Local branching heuristics after start heuristic
  neigh_size_hu::Int8 = 2;# k: Neighberhood size
  bound_node_hu::Int16 = 20; # node limit in bnb algorithm
  #____________Suboptimality settings___________
  subopt_sol::Bool = false; #true; 
  #___epsilont-suboptimality___
  subopt_epsilon::Bool = false; #relax dom cut
  tol_domcut_abs::Float64 = 1e-6;
  tol_domcut_rel::Float64 = 1e-6;
  tol_abs_param::Bool = false; #absoulte allowance parameter is parametric?
  tol_domcut_abs_param::Matrix{Float64} = [1e-6.*Matrix(I,2,2) 1e-6.*ones(2)] #eps = aJ: [Beps Ceps],dim: 1*(nth + 1) for MILP || [Aeps Beps' Ceps], dim:(nth)*(nth + 2) for MIQ
  #___T-cut method___
  subopt_Tcut::Bool = false; # #generated-node cut
  bound_tot_gen_nodes::Float64 = Inf # bounds on total generated nodes
  bound_decomp_nodes::Float64 = Inf # bounds on total generated nodes #Int32
  #___M-cut method___
  subopt_Mcut::Bool = false; # #active-node cut
  bound_active_nodes::Int32 = 10000 # bounds on stored active nodes
  find_subopt_err_MILP::Bool = true; # find suoptimal error by solving an MILP in each region, If false:solve relaxation instead
  #________ cert opts:__________
  eps_primal::Float64 = 1e-6
  eps_dual::Float64 = 0
  eps_zero::Float64 = 1e-10
  eps_gap::Float64  = 1e-6
  verbose::Int8 = 1 #2 #0
  iter_limit::Int64  = 1e3 #2
  storage_level::Int8 = 1
  max_constraints::Int64 = 1000
  delta_lam::Float64 = 0
  delta_mu::Float64 = 0
  delta_alpha::Float64 =0
  rm_callbacks::Vector{Function} = Function[] 
  add_callbacks::Vector{Function} = Function[]
  termination_callbacks::Vector{Function} = Function[]
end

mutable struct RegHeuristicFP  
  Ath::Matrix{Float64}
  bth::Vector{Float64} 
  F_h::Matrix{Float64}
  G_h::Vector{Float64}
  F_theta::Matrix{Float64}
  G_theta::Vector{Float64}
  A_theta::Matrix{Float64}
  B_theta::Matrix{Float64} 
  C_theta::Float64
  int_feas_found::Bool
  iter::Int64 
  bin_ind::Vector{Int64}
  loop::Int64 
  neigh_size::Int64 #k=neighberhood size: for local branching improvement heuristic
  exec::Bool #for local branching improvement heuristic
end

mutable struct RegHeuristic  
  Ath::Matrix{Float64}
  bth::Vector{Float64} 
  #ub_th::Vector{Float64} 
  #lb_th::Vector{Float64} 
  sol::ParamSol
  iter::Int64 
  loop::Int64 
  int_feas_found::Bool
  alphaa::Float64 #used for feasibility pump method
  G_rp::Vector{Float64}#used for feasibility pump method, previous rounded sol
  T::Vector{NodeMP} #used for diving method
  bin_ind::Vector{Int64}#used for diving method
  AS::Vector{Int64}
  lam::Matrix{Float64} #used for diving method
  neigh_size::Int64 #used for local branching #k=neighberhood size: for improvement heuristic
  exec::Bool #used for local branching improvement heuristic
end
