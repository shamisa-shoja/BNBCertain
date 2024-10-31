## _____________________________________________________________________
import Pkg; 
using Pkg
Pkg.activate(".")

include("BnBCertain.jl")
using .BnBCertain
using StatsBase, DataFrames, CSV, JLD2 

## _____________________________________________________________________
# Define the problem 
n = 8; nb = 4; nth = 2; m = n+8; bnd= 1; example = "rand"; # random example
prob, P_theta = BnBCertain.generate_mpMIQP(n,m,nth,nb,bnd,example); prob.senses = zeros(Cint,m+2*nb);

# Setting
opts = BNBSettings(); 
opts.problem_type_QP = false; #MILP or MIQP?
(!opts.problem_type_QP) && (prob.H = zeros(n,n); prob.f_theta = zeros(n,nth)); # For MILP case!
strategies = [BnBCertain.DepthFirst true] 

## _____________________________________________________________________
# Run certification with default settings
parts, max_iter, max_node, max_decomp_nodes, n_fin_reg = bnb_cert(prob, P_theta, opts)
println("worst-case iteration number= $(max_iter), worst-case number of nodes = $(max_node)")
