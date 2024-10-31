#certification branch and bound algorithm for MILP/MIQPs

function bnb_cert(mainprob::MPMIQP, P_theta, opts) 
  ## ____________________initialization______________________
  #bnb data
  n = length(mainprob.f); nth = size(mainprob.W,2); nb = length(mainprob.bin_ind);
  ub_th = P_theta.ub;  lb_th = P_theta.lb; 
  #__________
  (!opts.subopt_sol) && (opts.subopt_epsilon = false; opts.subopt_Mcut = false; opts.bound_decomp_nodes = Inf) #turn off suboptimality options
  #_________

  root_node = initialize_nodeMP(mainprob);
  T = NodeMP[];
  push!(T,root_node) 

  # Push initial region with its data to stack S
  S = BNBRegion[];
  final_partition = BNBRegion[];

  region0 = initialize_region(n, nb, P_theta, T);
  push!(S,region0);

  #Setup workspace for daqp...
  ws = ASCertain.setup_workspace(P_theta,1000); # TODO larger?
  certopts = CertSettings(); certopts.verbose=0; certopts.storage_level=2; certopts.max_constraints = 20000 #10000
  outer_loop = 0; max_iter = 0; max_size_tree = 0; max_decomp_nodes = 0; max_iter_hu = 0; max_iter_tot = 0; max_iter_hu_tot=0;
  #_________________________________________________________
  ## _______________________Main Loop ______________________
  while(!isempty(S)) 
    outer_loop+=1; 
    region = pop!(S);
    max_size_tree = max(max_size_tree,region.tot_nodes) #worst case size of tree 
    max_iter = max(max_iter,region.iter) #worst case accumulated iterations #region.kappa[1]  
    max_decomp_nodes = max(max_decomp_nodes,region.kappa.decomp_node) #worst number of decomposing nodes (branching) in the tree 
    max_iter_hu = max(max_iter_hu,region.kappa.iter_hu) #
    max_iter_tot = max(max_iter_tot,region.kappa.iter_tot) # total heuristics and bnb computations
    max_iter_hu_tot = max(max_iter_hu_tot,region.kappa.iter_hu_tot) #total heuristics and bnb computations

    ## ___________________________Terminated Region_______________________
    if (isempty(region.T) || (region.tot_nodes >= opts.bound_tot_nodes) || (region.kappa.decomp_node >= opts.bound_decomp_nodes)) #(region.tot_gen_nodes > opts.bound_dtot_gen_nodes)
      # region terminated
      ind_fin_reg = length(final_partition); #index of current region in final partitioning
      (opts.verbose>=2)&&println("#outer_loop: $(outer_loop), flag = PushtoF, #reg_F: $(ind_fin_reg+1)")
      region.kappa.loop = region.kappa.loop +1; region.kappa.time = round(time() - region.kappa.time);
      push!(final_partition, region)
    else 
      ## ___________________________Processing Region_______________________
      #pop new problem
      node = pop!(region.T); #node_reg = deepcopy(node); #for debug
      AS_parent = opts.warm_start ? deepcopy(node.AS) : Int64[]; #AS = findall(prob.senses.&(ACTIVE+IMMUTABLE).!=0) ∪ node.AS
      #_______________________________
      infeas_prob, region_infeas= isprobInfeasMP(node,region,outer_loop) #chech if the problem formulation is infeasible
      #_______________________________
      if (!infeas_prob)
        ## define dual problem
        (opts.problem_type_QP) ? (dualprob = DualCertProblem(node)) : (dualprob = ASCertain.DualLPCertProblem(node.f, node.A, [node.W'; node.b'], size(node.W,2), length(node.f), node.bounds_table, node.senses));
        ASCertain.reset_workspace(ws)
        #_______________________________
        #Remove Redundancy  
        Arth = region.Ath; brth = region.bth;  #work on the same region  
        if (opts.compute_minimal_rep) 
          n_cons = size(region.bth,1)
          if (n_cons >= 100)
            try
              p = polyhedron(hrep(Array(region.Ath'), region.bth), CDDLib.Library()) #, lib); #, CDDLib.Library() #plot(p,color = "red", alpha = 0.1)
              #p = polyhedron(hrep(Array(region.Ath'), region.bth), lib); #, Use Gurobi
              removehredundancy!(p); Pn = MixedMatHRep(p); 
              Arth = deepcopy(Pn.A); brth = deepcopy(Pn.b); Arth = Array(Arth'); 
              Pn = []; p = []; #clear up memory
            catch

            end
          end
        else
          #____________________
          Arth,brth = ASCertain.remove_redundant(region.Ath,region.bth) #USe certification tools to remove redundancy
        end
        reg_th = (A = deepcopy(Arth), b = deepcopy(brth), ub=ub_th, lb=lb_th); Arth = []; brth = []; #clear up memory
        ## ____________________________ Certification ____________________________
        ##
        partition,_ = certify(dualprob,reg_th,AS_parent,certopts);  # certify MILP/MIQP 
        sol_part_ld = primal_solution(node,dualprob,partition,opts.problem_type_QP) #lower dimention solution#extraxt primal solution and obj function from dual solution
        (opts.verbose>=2)&&print("\r>> #$(outer_loop) | Stack T: $(length(region.T)) | Stack S: $(length(S)) |Stack F: $(length(final_partition)) | partition: $(length(partition)) | max iter: $(max_iter) |");
        sol_part = update_sol_to_entire(sol_part_ld, mainprob, n, nth, node.bin_ind_all, node.bin_values,opts.problem_type_QP) #increase dim of sol to n
        sol_part_ld = []; #clear up memory
        #plot_partition(partition,inds=[1,2], axis=[P_theta.lb[1],P_theta.ub[1],P_theta.lb[2],P_theta.ub[2]]) #for debug       
        ## ___________________________Process on each Region _______________________
        ## evaluate cut conditions in each single region
        inner_loop= 0;
        for part in partition #inner_loop in eachindex(partition) # inner_loop= 1:length(partition) 
          inner_loop +=1; # part = partition[inner_loop] #part.iter
            #____________________________
            if (opts.remove_empty_reg) #(opts.compute_minimal_rep) 
                #remove empty regions
                Pn = polyhedron(hrep(Array(part.Ath'), part.bth), CDDLib.Library());#npoints(Pn); #polyhedron and plot
                empty_reg = false;
                try; empty_reg =  (isempty(Pn) || (dim(Pn) != nth) || (fulldim(Pn) != nth)) #(npoints(Pn) == 0)
                catch; empty_reg = true; end #
                if empty_reg
                continue; # empty or lower dimentional region
                end
                Pn = []; #clear up memory
            end
            #___________________________
            S = eval_cut_cond(part,region,sol_part[inner_loop],node,inner_loop,outer_loop,S,opts); #length(S)
         #______________________________
        end # for each region
      else
        push!(S,region_infeas) #if infeasibile problem
      end # if !infeas_prob
    end # if T empty
  end #main loop
  ## _________________________________________________________
  (opts.verbose>=2)&&print("\r>> ** Cert BnB terminate **| The worst acc iter: $(max_iter)|"); 
   return final_partition, max_iter, max_size_tree, max_decomp_nodes, length(final_partition), max_iter_hu #, info_cert_all
end #function
## ________________________________________________________________________________________________
## ___________________________________________________________________________
function eval_cut_cond(part,region,sol_partt,node,inner_loop,outer_loop,S,opts)
    #______________________________________
    new_kappa = deepcopy(region.kappa);
    new_kappa.iter = part.iter + region.iter; new_kappa.node = new_kappa.node+1;  new_kappa.loop = new_kappa.loop+1; new_kappa.iter_tot = part.iter + new_kappa.iter_tot; new_kappa.iter_hu_tot = part.iter + new_kappa.iter_hu_tot;
    new_region = BNBRegion(deepcopy(part.Ath),deepcopy(part.bth),part.iter + region.iter, part.AS, deepcopy(region.sol_bar),deepcopy(region.T),
        new_kappa,[region.seq_nodes node.bin_values], deepcopy(region.tree_status),region.tot_nodes+1, deepcopy(region.info_status));  #(part.iter + region.iter_avg)/2, 
    #___________________________________________________________________
    ##handle infeasible and unbounded solution immediately
    dom_cut = false;  int_feas = false; flag = BNBState[];
    if string(part.state) == "INFEASIBLE"  
         dom_cut = true; flag = InfeasCut;
    elseif string(part.state) == "UNBOUNDED"  
         dom_cut = true; flag = Unbounded; #check if cut unbounded solution
    end
    #______________________________
    # Dominance and Infeasibility Cut conditions
    if !dom_cut #&&(opts.check_dom_cut)
        dom_cut, flag, reg_part, reg_part_prune =  isdominated(part,sol_partt,new_region.sol_bar, opts)
        (dom_cut) && (flag = DomCut)     
    end
    #_________________Dominance Cut________________      
    if (dom_cut)
        # Dominance cut (includes infeasibility cut)
        (opts.verbose>=2)&&println("#outer_loop: $(outer_loop), flag = $(flag), #inner_loop: $(inner_loop) / $(length(partition))")  
        new_region.tree_status = [new_region.tree_status; flag node.level]; 
        new_region.info_status = [new_region.info_status; outer_loop inner_loop node.level flag node.bin_values']; #[outer_loop inner_loop level bnbstate binaryvalues] 
        (flag == Unbounded) && (new_region.AS = Int64[]); 
        push!(S,new_region)
        return S; #continue; 
    end # if dominated
    #____________________________
    if flag == PartialDomCut
        # partially prune part of a region
        (opts.verbose>=2)&&println("#outer_loop: $(outer_loop), flag = $(PartialDomCut), #inner_loop: $(inner_loop) / $(length(partition))")  
        new_region_prune = deepcopy(new_region);
        new_region_prune.tree_status = [new_region_prune.tree_status; DomCut node.level];
        new_region_prune.info_status = [new_region_prune.info_status; outer_loop inner_loop node.level flag node.bin_values']; #[outer_loop inner_loop level bnbstate binaryvalues]           
        new_region_prune.Ath = reg_part_prune.Ath; #shrink region
        new_region_prune.bth = reg_part_prune.bth;
        push!(S,new_region_prune)
        # update current region
        new_region.Ath = reg_part.Ath; #shrink region
        new_region.bth = reg_part.bth;
    end # if partially dominated
    #_________________________________________________________
    # Integer Feasibility
    int_feas = isIntfeasMP(node, part.AS); 
    #_________________________________________________________
    #____________________ Integer Feasiblity Cut _____________
    if (int_feas || isempty(node.bin_ind))
        flag = IntFeasCut; (opts.verbose>=2)&&println("#outer_loop: $(outer_loop), flag = $(flag), #inner_loop: $(inner_loop) / $(length(partition))") 
        new_region.tree_status = [new_region.tree_status; flag node.level]; 
        new_region.info_status = [new_region.info_status; outer_loop inner_loop node.level flag node.bin_values']; 
        (!opts.problem_type_QP) && (new_region.sol_bar = ParamSol[])  #MILP #just store one sol_bar #erase stored upper bound 
        push!(new_region.sol_bar,sol_partt)   #update best-known solution         
        push!(S,new_region) #push in the end
        return S
    end #else 
    #____________________Branching _______________
    flag = Branch; (opts.verbose>=2)&&println("#outer_loop: $(outer_loop), flag = $(flag), #inner_loop: $(inner_loop) / $(length(partition))")
    priority_order = ordering_costMP(node.level,opts.sorting_method,sol_partt);   
    br_ind, reg_branch, flag_br_once= branchingScoreMP(sol_partt, node.bin_ind_all, node.bin_values, opts.branching_rule, part, opts); # #find ind of branching value #1; #          
    #_______________________  
    branch_ind = br_ind; #println(outer_loop); #for debug # Just one binary variable to branch on
    node0,node1 = branchingMP(node,branch_ind,priority_order,part.AS, part.Lam, opts);
    new_region.tree_status = [new_region.tree_status; flag node.level]; 
    new_region.info_status = [new_region.info_status; outer_loop inner_loop node.level flag node.bin_values']; #[outer_loop inner_loop level bnbstate binaryvalues] 
    child_nodes = NodeMP[];
    (node1.feasible) && push!(child_nodes, node1); #push if node is feasible by warm-starting
    (node0.feasible) && push!(child_nodes, node0); 
    S = appendingMP(new_region, child_nodes, S, opts); 
    return S
    #_____________________
end