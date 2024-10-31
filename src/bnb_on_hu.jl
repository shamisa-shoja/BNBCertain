#online branch and bound algorithm
#mainprob = MIQP(prob, th_sam) #for test
function bnb_on_hu(mainprob::MIQP, opts::BNBSettings) 
  #n = length(mainprob.f);cont_ind = setdiff([1:n;],mainprob.bin_ind); bin_ind = mainprob.bin_ind; 

  #initialization
  root_node = initialize_node(mainprob);
  T = Node[]; 
  push!(T,root_node) 

  Jbar  = Inf;  xbar = NaN; 
  seq_nodes = 5*ones(Int64,length(mainprob.bin_ind)); 
  tree_status = [Root 0]; info_status = [0 0 0 Root root_node.bin_values']
  tot_nodes = 0; total_nodes = 2^(length(mainprob.bin_ind)+1)-1;
  iter = 0; out_loop = 0; decomp_nodes = 0; tot_gen_nodes = 1; iter_hu = 0; iter_tot = 0; iter_hu_tot = 0; #eps_dom_abs = opts.tol_domcut_abs; eps_dom_rel = opts.tol_domcut_rel;
  (!opts.subopt_sol) && (opts.subopt_epsilon = false; opts.subopt_Mcut = false; opts.bound_decomp_nodes = Inf) #turn off suboptimality options

  #____________________________

  while(!isempty(T) && (tot_nodes < opts.bound_tot_nodes) && (decomp_nodes < opts.bound_decomp_nodes))  
    out_loop+=1;    
    node = pop!(T); #pop new problem
    AS_parent = opts.warm_start ? deepcopy(node.AS) : Int64[]
    #_______________________________
    infeas_cut = false;
    for i = 1:size(node.A,1)        
      if(all(node.A[i,:] .== 0) && node.b[i] <0) 
        # THE NODE IS INFEASIBLE, Cut without solving
        infeas_cut = true; flag = InfeasCut; #J = Inf; x = NaN; AS_opt = Int64[];
        tot_nodes += 1; seq_nodes = [seq_nodes node.bin_values];
        tree_status = [tree_status; flag node.level];  info_status = [info_status; out_loop 0 node.level flag node.bin_values']; #[out_loop j level bnbstate binaryvalues]
        break  
      end
    end
    #_______________________________
    if (!infeas_cut)        
     ## solve the relaxation    
      if opts.daqp_solver_method
        if opts.problem_type_QP #MIQP         
          x,lam,AS_opt,J,iter_opt,ASs = DAQP.daqp_jl(node,AS_parent);  #DAQP.daqp_jl(node.H,node.f,node.A,node.b,node.sense,AS_parent); 
        else #MILP
          x,lam,exitflag,AS_opt,iter_opt  = ASCertain.dsimplex(node.f,node.A,node.b,AS_parent) 
          #@time (flag_g, xg, Jg, ASg) = QLPSolve(node, opts) #Gurobi
          #________________
          if (exitflag == 2); (!isempty(x)) ? (J = -lam'*node.b[AS_opt]) : (J = 0) #J = node.f'*x #optimal
          elseif (exitflag == 3); J = Inf #infeasible 
          else; J = -Inf; end #unbounded
          #____________
        end #MILP or MIQP
      else
        opts.chosen_solver == "GLPK"
        @time (exitflag1, x1, J1, AS_opt1) = QLPSolve(node,opts); iter_opt = 1; #number of nodes
      end #solver
  
      ## update bnb data
      J = J + node.extra_cost; #(opts.add_extra_cost) && (J = J + node.extra_cost);
      iter = iter + iter_opt;
      iter_tot = iter_tot + iter_opt; iter_hu_tot = iter_hu_tot + iter_opt;
      tot_nodes += 1; seq_nodes = [seq_nodes node.bin_values];
      (opts.verbose>=2)&&print("\r>> #proc./tot. nodes: $(tot_nodes)/$(total_nodes) | Stack T: $(length(T))| iter: $(iter)|");

      #________________________________________
      #EVALUATE CUT conditions
      T,xbar,Jbar,iter,iter_hu,iter_tot,iter_hu_tot,seq_nodes,tree_status,info_status,decomp_nodes,tot_gen_nodes = eval_cut_cond_on_hu(x,J,exitflag,T,xbar,Jbar,node,mainprob,AS_opt,lam,iter,iter_hu,iter_tot,iter_hu_tot,seq_nodes,tree_status,info_status,tot_gen_nodes,decomp_nodes,out_loop,opts)
    end #if !infeas_cut
  end # while T !empty
  #____________________________
  (opts.verbose>=2)&&print("\r>> ======= BnB terminate ========= | total iter: $(iter)|");
  return xbar,Jbar, iter, iter_hu, iter_tot, tot_nodes, seq_nodes, tree_status, info_status, decomp_nodes, tot_gen_nodes, iter_hu_tot #, [], [] #seq_status #seq_nodes[:,2:end]
end #function

## ___________________________________________________________________________
############################################################################
   ############################ Functions ################################
############################################################################
function initialize_node(mainprob::MIQP)
  n = size(mainprob.f,1);
  m = size(mainprob.b,1);
  nb = size(mainprob.bin_ind,1);
  bin_ind = mainprob.bin_ind;#cont_ind = setdiff([1:n;],mainprob.bin_ind);
  # Add 0≤xb ≤ 1 constraints 
  eye =Matrix(1.0I,n,n)
  A_bin = [-eye[mainprob.bin_ind,:];eye[mainprob.bin_ind,:]];
  b_bin = [zeros(nb);ones(nb)];
  A = [mainprob.A; A_bin]; 
  b = [mainprob.b; b_bin]; 
  # Specifiy constraint coupeling and types
  nb = length(bin_ind);
  bounds_table = collect(1:m+2*nb); #bounds_table = []; 
  senses = zeros(Cint,m+2*nb);#repeat([INEQUALITY],m+2*nb); 
  nb = length(bin_ind);
  bin_values = -ones(Int64,nb);
  bin_indices = Matrix{Int64}(undef,nb,2); 
  bin_indices[:,1] = (m+1):(m+nb);
  bin_indices[:,2] = (m+nb+1):(m+2*nb);
  bin_states = repeat([FREE],nb); 
  bin_ind_all = deepcopy(bin_ind); #all binary variables
  AS= Int64[];
  priority = 0;
  level = 0;
  extra_cost = 0;
  feasible = true; #when warm-starting

  root_node = Node(mainprob.H,mainprob.f,A,b,bounds_table,senses,mainprob.bin_ind,bin_values,bin_indices,bin_states,n,nb,bin_ind_all,AS,priority,level,extra_cost,feasible);
  return root_node
end
#_________________________________________________________
function isIntfeas(node, AS) #::MIQP,::Vector{Int64}
  if typeof(AS) != BitVector
    AS1 = falses(size(node.b,1));
    AS1[AS] .= true;
    AS = AS1;
  end
  # find if the solution is integer feasible
  AS_bin = [AS[node.bin_indices[:,1]]; AS[node.bin_indices[:,2]]];
  leaf_node = isempty(node.bin_ind);
  int_feas = leaf_node || (length(findall(AS_bin)) >= length(node.bin_ind));
  return int_feas 
end  
#_________________________________________________________
function update_xbar(x,n,bin_ind,bin_ind_curr,bin_values)
  cont_ind = setdiff([1:n;],bin_ind); 
  x_tot = zeros(n); 
  x_bin = x[bin_ind_curr]; 
  x_tot[cont_ind] = x[cont_ind];  

  if size(x,1) == n
    x_tot[bin_ind] = round.(x_bin)
  else   
    bin_fixed = findall(bin_values .> -1)
    bin_free = findall(bin_values .== -1)
    x_tot[bin_ind[bin_fixed]] = bin_values[bin_fixed]       
    if !isempty(bin_ind[bin_free])
      x_tot[bin_ind[bin_free]] = round.(x_bin);
    end
  end
  return x_tot
end
#_________________________________________________________
function branchingScore(x,bin_ind, branchrule)
  # bin_ind=node.bin_ind; branchrule =opts.branching_rule; # for debug
  xbin = x[bin_ind];
  if branchrule == FirstBin
    branch_ind = 1; 
  elseif branchrule == LastBin
    branch_ind = size(bin_ind,1);
  elseif branchrule == MinBin
    branch_ind = argmin(xbin); #For integer case: # (abs.(xbin-round.(xbin))); 
  elseif branchrule == MaxBin
    branch_ind = argmax(xbin); #For integer case: # (abs.(xbin-round.(xbin)));
  elseif branchrule == MostInfeasBin #closest to 0.5
    branch_ind = argmax(abs.(0.5.-abs.(xbin.-0.5))); #argmin(abs.(xbin .- 0.5)); #For integer case: # (abs.(xbin-round.(xbin) .- 0.5)); 
  else
    branch_ind = 1;  #default: in order
  end
  return branch_ind
end 

#________________________________________________________________________
#_________________________________________________________
function ordering_cost(level,sorting,J)
  if sorting == DepthFirst
    priority_order = 1 / (level+1); 
  elseif sorting == BreadthFirst
    priority_order = level+1;
  elseif sorting == BestFirst
    priority_order = J;
  else
    priority_order = 1 / (level+1); #default:"depth_first"
  end
  return priority_order
end
#_________________________________________________________
function branching(prob, br_ind::Int64, priority_order, AS, lam, opts) #::MIQP,::Vector{Int64}
  n = length(prob.f);
  #indices = [collect(1:this-1), collect(this+1:n)] #[1:this-1, this+1:n];
  branch_ind = prob.bin_ind[br_ind]
  indices = [1:n;];
  indices = setdiff(indices, branch_ind);
  prob1 = deepcopy(prob);#fieldnames(typeof(prob)) #see the field of the struct
  prob0 = deepcopy(prob);

  H11 = prob.H[indices, indices];
  H12 = prob.H[indices, branch_ind];
  H22 = prob.H[branch_ind, branch_ind];
  prob0.H = H11; prob1.H = H11;

  f1  = prob.f[indices];
  f2  = prob.f[branch_ind]; 
  prob0.f = f1;
  prob1.f = f1[:]+H12[:];

  if isdefined(prob,:f_theta)
    f_theta = prob.f_theta[indices,:];
    prob0.f_theta = f_theta; prob1.f_theta = f_theta; 
  end

  m = size(prob.A,1);
  rm_bin_indices = prob.bin_indices[br_ind, :];
  const_ind = setdiff([1:m;],rm_bin_indices); #[bin_indices[:,1]; bin_indices[:,2]] #indices; 
  A1 = prob.A[:,indices]; # remove the column for removed xb
  A = A1[const_ind,:];    # remove rows for binary constraints of xb
  prob0.A = A; prob1.A = A; 

  b = prob.b[const_ind];
  prob0.b = b; prob1.b = b - prob.A[const_ind,branch_ind];

  if isdefined(prob,:W)
    W = prob.W[const_ind,:];
    prob0.W = W; prob1.W = W;
  end

  prob0.extra_cost = prob.extra_cost; 
  prob1.extra_cost = prob1.extra_cost + .5*H22 + f2; #add constant terms in value function

  #_________Add just current extra cost____________
  #prob0.extra_cost = 0;
  #prob1.extra_cost = .5*H22 + f2; #add current constant terms in value function

  if isdefined(prob,:extra_cost_param)
    extra_cost_param = prob.f_theta[branch_ind,:];#add parametric terms in value function
    #_________Add just current extra cost____________
    #prob0.extra_cost_param = []; prob1.extra_cost_param = []; #Matrix{Float64}(I,nth,nth) 
    prob0.extra_cost_param = zeros(size(extra_cost_param,1),1);
    prob1.extra_cost_param = zeros(size(extra_cost_param,1),1)+ extra_cost_param; 
  end

  bin_ind = [prob.bin_ind[1: br_ind-1]; prob.bin_ind[br_ind+1:end].-1]; #remove fixed binary var
  prob0.bin_ind = bin_ind; prob1.bin_ind = bin_ind; 

  bin_indices = prob.bin_indices[1:end .!= br_ind, :];
  bin_indices[:,1] = [prob.bin_indices[1: br_ind-1, 1]; prob.bin_indices[br_ind+1:end, 1].-1]; 
  bin_indices[:,2] = [prob.bin_indices[1: br_ind-1, 2].-1; prob.bin_indices[br_ind+1:end, 2].-2]; 
  prob0.bin_indices = bin_indices; prob1.bin_indices = bin_indices; 

  prob0.bin_states[br_ind] = LOWER;
  prob1.bin_states[br_ind] = UPPER;
  n_const = size(prob0.b,1);
  prob0.bounds_table = collect(1:n_const);
  prob0.senses = zeros(Cint,n_const);
  prob1.bounds_table = collect(1:n_const);
  prob1.senses = zeros(Cint,n_const);

  if (opts.warm_start && !isempty(AS))
    AS_updated0, n0_feas, AS_updated1, n1_feas = basis_recovery(prob,prob0, AS,br_ind,lam) #use (half of) simplex method to update AS
    prob0.AS = AS_updated0; prob0.feasible = n0_feas;
    prob1.AS = AS_updated1; prob1.feasible = n1_feas;
  else
    prob0.AS = copy(AS); 
    prob1.AS = copy(AS);
  end

  #bin_fixed = findall(prob.bin_values .> -1);
  bin_free = findall(prob.bin_values .== -1);
  branch = bin_free[br_ind];
  prob0.bin_values[branch] = 0;
  prob1.bin_values[branch] = 1;

  prob0.level = prob0.level + 1;  
  prob1.level = prob1.level + 1;
  prob0.priority = priority_order;
  prob1.priority = priority_order;
 
  return prob0,prob1
end
#_________________________________________________________
function basis_recovery(node,node0, AS_parent,br_ind,lam)
  #br_ind = node.bin_ind[br_indx];
  m = node.bin_indices[1,1]-1; #ineq constraints without binary constraints
  nb = length(node.bin_ind); n = size(node0.f,1);
  AS_ineq = AS_parent[AS_parent .<= m];  #ActiveSet of original inequality constraints  
  bin_ind0 = m .< AS_parent .< (m+nb+1);        
  AS_b_lb = AS_parent[bin_ind0];         #AS of lower bound of left binary constraints
  AS_b_ub = AS_parent[AS_parent .> m+nb];
  #remove index of branchvar and branchvar+nb from AS
  AS_bin_mod = [AS_b_lb[AS_b_lb .< m+br_ind]; AS_b_lb[AS_b_lb .> m+br_ind].-1; AS_b_ub[AS_b_ub .< m+nb+br_ind].-1; AS_b_ub[AS_b_ub .> m+nb+br_ind].-2];
  AS_updated = [AS_ineq; AS_bin_mod];

  if (length(AS_updated) > n)
    #apply (half of) simplex method to remove one index from AS (in simplex: lenght(AS) = n)
    n0_feas = true; n1_feas = true; 
    fix_id_0,fix_id_1 = node.bin_indices[br_ind,:];
    #___________________________
    AS0 = fix_warm_AS(node,AS_parent,lam,fix_id_0);  
    if !isempty(AS0)
      AS_b_lb0 = AS0[m .< AS0 .< (m+nb+1)];        #AS of lower bound of left binary constraints
      AS_b_ub0 = AS0[AS0 .> m+nb];#AS of lower bound of left binary constraints
      AS_updated0 = [AS0[AS0 .<= m]; AS_b_lb0[AS_b_lb0 .< (m+br_ind)]; AS_b_lb0[AS_b_lb0 .> (m+br_ind)].-1; AS_b_ub0[AS_b_ub0 .< m+nb+br_ind].-1; AS_b_ub0[AS_b_ub0 .> m+nb+br_ind].-2];
    else
      AS_updated0 = Int64[]; n0_feas = false
    end
    #___________________________
    AS1 = fix_warm_AS(node,AS_parent,lam,fix_id_1);
    if !isempty(AS1)
      AS_b_lb1 = AS1[m .< AS1 .< (m+nb+1)];         #AS of lower bound of left binary constraints
      AS_b_ub1 = AS1[AS1 .> m+nb];#AS of lower bound of left binary constraints
      AS_updated1 = [AS1[AS1 .<= m]; AS_b_lb1[AS_b_lb1 .< (m+br_ind)]; AS_b_lb1[AS_b_lb1 .> (m+br_ind)].-1; AS_b_ub1[AS_b_ub1 .< m+nb+br_ind].-1; AS_b_ub1[AS_b_ub1 .> m+nb+br_ind].-2];
    else
      AS_updated1 = Int64[]; n1_feas = false
    end
    #___________________________
    return AS_updated0, n0_feas, AS_updated1, n1_feas
    #___________________________
  else #no need to call simplex method
    return AS_updated, true, AS_updated, true
  end
end
#_________________________________________________________

function appending(T,childs,opts)
  if isempty(childs)    
    return T; # no new node to be pushed in T
  end
  #____________________________________
  
  if isempty(T) #T empty, no need to find a place, insert in the end
    sortedT = childs;
    if (opts.subopt_Mcut && (length(sortedT) > max_act_nod))
      ind_start = length(sortedT)- max_act_nod + 1;
      sortedT = [sortedT[ind_start:end];]; #sortedT = [sortedT[1:max_act_nod];]; #bests are in the end       
    end
    return sortedT;
  end
  #____________________________________
  #ordering nodes: higher priority (smaller rho): put in the end of T (first popped)
  sorting = opts.sorting_method; 
  max_act_nod = opts.bound_active_nodes; #for active-node subopt cut
  #_______________________
  if sorting == DepthFirst     
    sortedT = [T[1:end]; childs[1:end]];  #return sortedT
    #____________________________________
  elseif sorting == BreadthFirst    
    sortedT = [childs[1:end];  T[1:end]]; #return sortedT
    #________________________________
  else #if sorting == BestFirst #and other node selection strategies
    sortedT = NodeMP[];   
    prior_nodes = zeros(length(T));#node0 = childs[end];
    for i = length(T):-1:1     
      prior_nodes[i] = T[i].priority
    end
    ind1 = findall((childs[end].priority.- prior_nodes) .<= opts.eps_zero);
    ind = isempty(ind1) ? 0 : ind1[end];
    sortedT = [T[1:ind]; childs[1:end]; T[ind+1:end]]   
  end #if best_first
  if (opts.subopt_Mcut && (length(sortedT) > max_act_nod))
    ind_start = length(sortedT)- max_act_nod + 1;
    sortedT = [sortedT[ind_start:end];]; #sortedT = [sortedT[1:max_act_nod];]; #bests are in the end       
  end
  return sortedT; 
end # function

function eval_cut_cond_on_hu(x,J,exitflag,T,xbar,Jbar,node,mainprob,AS_opt,lam,iter,iter_hu,iter_tot,iter_hu_tot,seq_nodes,tree_status,info_status,tot_gen_nodes,decomp_nodes,out_loop,opts)
  ## Evaluate 3 cut conditions and check heuristics
  n= node.n; #bin_ind = deepcopy(node.bin_ind_all); cont_ind = setdiff([1:n;],node.bin_ind_all); 
  #________________________________________________
    dom_cut = false; int_feas = false; flag = BNBState[]; int_sol_found_hu = false;
    if (exitflag == 3) #string(exitflag) == "INFEASIBLE"  
        dom_cut = true; flag = InfeasCut;
    elseif (exitflag != 2) #string(exitflag) == "UNBOUNDED"  
        dom_cut = true; flag = Unbounded; #check if cut unbounded solution
        AS_opt = Int64[]; 
    end 
    #____________________________
    (!dom_cut) && (int_feas = isIntfeas(node, AS_opt)) #is integer feasible
    #____________________________
    ## Apply heuristics
    trigger_hu = any(out_loop .== opts.hu_trigger)
    if (opts.apply_heuristic && !dom_cut && !int_feas  && trigger_hu) #start heuristic
        xbar,Jbar,iter,iter_hu,iter_tot,iter_hu_tot = heuris_start_on(x,J,xbar,Jbar,node,mainprob,AS_opt,lam,iter,iter_hu,iter_tot,iter_hu_tot, opts)
    end #if heuristic   
    #_____________________________
    ## Evaluate Cut conditions
    if !dom_cut # && opts.check_dom_cut)
        dom_cut = (opts.subopt_epsilon) ? (opts.tol_domcut_abs + J >= Jbar) : (J >= Jbar) #(1+ eps_dom_rel)*# J - Jbar
        (dom_cut) && (flag = DomCut)
    end 
    #_______________________________
    if(dom_cut) #-opts.eps_zero Dominance cut 
        # Dominance and Infeasibility Cut conditions
        tree_status = [tree_status; flag node.level]; info_status = [info_status; out_loop 0 node.level flag node.bin_values']; #[out_loop j level bnbstate binaryvalues] 
        #continue;
        #______________________________________
        return T, xbar, Jbar, iter, iter_hu, iter_tot, iter_hu_tot, seq_nodes, tree_status, info_status, decomp_nodes, tot_gen_nodes
    end    
    #___________________________________
    # better solution
    if((int_feas) || isempty(node.bin_ind))  #for leaf node with unbounded case 
        Jbar = J;
        xbar = update_xbar(x,n,node.bin_ind_all,node.bin_ind,node.bin_values); #update solution to the real dimention n 
        tree_status = [tree_status; IntFeasCut node.level]; #seq_status = [seq_status; "int cut"];
        info_status = [info_status; out_loop 0 node.level IntFeasCut node.bin_values']; #[out_loop j level bnbstate binaryvalues] 
        #___________________________
        #improvement heuristic: local branching
        if (opts.apply_heuristic && opts.heuris_imp) # && opts.heuristic_LB) #&& heu_first_call
            if opts.heuris_imp_type == "LB" #
                int_sol_foundLB, x_huLB, J_hLB, iter_hLB, lLB = heuristics_LB(xbar, Jbar, mainprob, node, opts) 
                if int_sol_foundLB
                    (J_hLB < Jbar) && (xbar = deepcopy(x_huLB); Jbar = deepcopy(J_huLB); iter_hu = iter_hu + iter_hLB); #Jbar = node.f'*xbar; #J_h
                    (opts.add_iter_hu) && (iter = iter + iter_hLB);
                end #if int_sol_foundLB
            end # if LB
        end #if improvement heuristic
        #____________________________
    else #do branching   
        ##TODO: uncooment function in selecting br_ind
        br_ind = branchingScore(x, node.bin_ind, opts.branching_rule); #allow for different branching rules #1; #
        priority_order = ordering_cost(node.level,opts.sorting_method,J);         
        node0,node1 = branching(node, br_ind, priority_order, AS_opt, lam, opts);
        child_nodes = Node[];
        (node1.feasible) && push!(child_nodes, node1); #push if node is feasible by warm-starting
        (node0.feasible) && push!(child_nodes, node0); #length(child_nodes) # node0.AS 
        curr_act_nod = length(T);
        T = appending(T, child_nodes, opts); #length(T) #push!(T,node1,node0);
        tot_gen_nodes = tot_gen_nodes + length(T) - curr_act_nod;  #for GN-cut subopt method #maybe some child nodes are cut away     
        decomp_nodes = decomp_nodes + 1;
        tree_status = [tree_status; Branch node.level];  info_status = [info_status; out_loop 0 node.level Branch node.bin_values']; #[out_loop j level bnbstate binaryvalues] #seq_status = [seq_status; "branch"];
    end #if int feas
    #______________________________________
    return T, xbar, Jbar, iter, iter_hu, iter_tot, iter_hu_tot, seq_nodes, tree_status, info_status, decomp_nodes, tot_gen_nodes
end

function heuris_start_on(x,J,xbar,Jbar,node,mainprob,AS_opt,lam,iter,iter_hu,iter_tot,iter_hu_tot, opts)
    int_sol_found_hu = false;
    #______ Feasibility Pump_______
    if (opts.heuris_start_type == "FP") #(opts.heuristic_FP)
        int_sol_found_hu, x_hu, iter_h, l_h, J_h, AS_h = heuristics_FP(x, node, AS_opt, opts) 
        (int_sol_found_hu) && (xbar = deepcopy(x_hu); Jbar = node.f'*xbar; iter_hu_tot = iter_hu_tot + iter_h); #Jbar = node.f'*xbar; 
        (int_sol_found_hu && opts.add_iter_hu) && (iter = iter + iter_h); 
        iter_hu = iter_hu + iter_h; iter_tot = iter_tot + iter_h; # maximum(node.A*x_hu -node.b)
        #___________________
    elseif (opts.heuris_start_type == "Div") #(opts.heuristic_Div)
        #______ Diving_______
        int_sol_found_hu, x_huD, iter_hD, lD, J_hD, AS_hD = heuristics_div(x, J, Jbar, node, AS_opt, lam, opts) 
        (int_sol_found_hu) && (xbar = deepcopy(x_huD); Jbar =  node.f'*xbar;iter_hu_tot = iter_hu_tot + iter_hD); #Jbar = node.f'*xbar; #J_h
        (int_sol_found_hu && opts.add_iter_hu) && (iter = iter + iter_hD); 
        iter_hu = iter_hu + iter_hD; iter_tot = iter_tot + iter_hD;# maximum(node.A*x_huD -node.b)
    elseif (opts.heuris_start_type == "RENS") #(opts.heuristic_RENS)
        #______ RENS_______
        int_sol_found_hu, x_huR, J_hR, iter_hR = heuristics_rens(x, J, node, mainprob, opts) 
        (int_sol_found_hu) && (xbar = deepcopy(x_huR); Jbar =  node.f'*xbar; iter_hu_tot = iter_hu_tot + iter_hR); #Jbar = node.f'*xbar; #J_h
        (int_sol_found_hu && opts.add_iter_hu) && (iter = iter + iter_hR); 
        iter_hu = iter_hu + iter_hR;  iter_tot = iter_tot + iter_hR; 
    end #if heuristic_FP 
    #___________________________improve the solution________________
    #improvement heuristic: local branching
    if (int_sol_found_hu && opts.heuris_imp) # && opts.heuristic_LB) #&& heu_first_call
        if opts.heuris_imp_type == "LB" #
            int_sol_foundLB, x_huLB, J_hLB, iter_hLB, lLB = heuristics_LB(xbar, Jbar, mainprob, node, opts) 
            if int_sol_foundLB
                (J_hLB < Jbar) && (xbar = deepcopy(x_huLB); Jbar = deepcopy(J_huLB); iter_hu = iter_hu + iter_hLB); #Jbar = node.f'*xbar; #J_h
                (opts.add_iter_hu) && (iter = iter + iter_hLB);
            end #if int_sol_foundLB
        end # if LB
    end #if improvement heuristic
    #_______________________________
    return xbar,Jbar,iter,iter_hu,iter_tot,iter_hu_tot
end