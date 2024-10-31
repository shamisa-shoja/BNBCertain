## _____________________________________________________________________________________________________
############################ Functions called in BnB Certification code ################################
##_____________________________________________________________________________________________________
function initialize_nodeMP(prob::MPMIQP)
    n = size(prob.f,1);
    m = size(prob.b,1);
    nb = size(prob.bin_ind,1);
    nth = size(prob.W,2);
    bin_ind = prob.bin_ind;  #cont_ind = setdiff([1:n;],prob.bin_ind);
    # Add 0≤xb ≤ 1 constraints 
    eye =Matrix(1.0I,n,n)
    A_bin = [-eye[prob.bin_ind,:];eye[prob.bin_ind,:]];
    b_bin = [zeros(nb);ones(nb)];
    A = [prob.A; A_bin]; 
    b = [prob.b; b_bin]; 
    W = [prob.W; zeros(2*nb,size(prob.W,2))]; 
    # Keep track on which constraints in new QP corresponds to binary relaxations
    bounds_table = collect(1:m+2*nb); #bounds_table = []; 
    senses = zeros(Cint,m+2*nb);#senses = repeat([INEQUALITY],m+2*nb); 
    bin_values = -ones(Int64,nb);
    bin_indices = Matrix{Int64}(undef,nb,2); 
    bin_indices[:,1] = (m+1):(m+nb);
    bin_indices[:,2] = (m+nb+1):(m+2*nb);
    bin_states = repeat([FREE],nb);  # Specifiy constraint coupeling and types
    bin_ind_all = deepcopy(bin_ind); #all binary variables
    AS= Int64[];
    priority = 0;
    level = 0;
    extra_cost = 0;
    extra_cost_param = zeros(nth,1);
    priority_BF = zeros(1,nth+1)
    feasible = true;
  
    root_node = NodeMP(prob.H,prob.f,prob.f_theta,prob.H_theta,A,b,W,bounds_table,senses,bin_ind,bin_values,bin_indices,bin_states,n,nb,bin_ind_all,AS,priority,level,extra_cost,extra_cost_param,priority_BF,feasible) 
    return root_node
end
## _______________________________________________________
function initialize_region(n, nb, P_theta, T)
    nth = size(P_theta.ub,1);
    # Pack initial region in the form A'θ ≤ b
    Ath = [Matrix{Float64}(I,nth,nth) Matrix{Float64}(-I,nth,nth) P_theta.A];
    bth = [P_theta.ub;-P_theta.lb;P_theta.b];
    iter0 = 0;
    AS = Int64[];
    sol_bar = ParamSol[];  
    F_theta = zeros(n,nth);
    G_theta = NaN*ones(n); #solution
    A_theta = zeros(nth,nth);
    B_theta = zeros(1,nth);
    C_theta = Inf; #value function
    sol_bar0 = ParamSol(F_theta, G_theta, A_theta, B_theta, C_theta);
    push!(sol_bar,sol_bar0)
    t_s = time();
    kappa0 = ComplexityMeasure(0,0,0,0,1,0,0,t_s,0,0) #kappa0 = [0]; #Dict();
    seq_nodes = 5*ones(Int64,nb,1);
    tree_status = [Root 0] 
    tot_nodes = 0; 
    #decomp_nodes = 0; tot_gen_nodes = 1; iter_hu = 0; iter_avg = 0; tree_avg = 0; ind = [0; 0]; #[k;j], where in the stack it's located
    info_status = [0 0 0 Root -ones(Int64,1,nb)]; 
    
    region0 = BNBRegion(Ath, bth, iter0, AS, sol_bar, T, kappa0, seq_nodes, tree_status, tot_nodes, info_status); #ind); #, size_avg, seq_status,
    return region0
end  

## ___________________________________________________________________________
function isprobInfeasMP(node,region,outer_loop)
  infeas_prob = false; region_infeas = [];
  for i = 1:size(node.A,1)        
    if(all(node.A[i,:] .== 0) && all(node.W[i,:] .== 0) && node.b[i] <0) #if(all(node.A[i,:] == 0) && all(node.b[1:end-1,i] == 0) && node.b[end,i] <0)
      # THE NODE IS INFEASIBLE, Cut without solving
      infeas_prob = true; flag = InfeasCut;
      region_infeas = deepcopy(region); 
      region_infeas.seq_nodes = [region_infeas.seq_nodes node.bin_values]; 
      region_infeas.tot_nodes += 1;region_infeas.AS = Int64[];
      region_infeas.tree_status = [region_infeas.tree_status; flag node.level];  
      region_infeas.info_status = [region_infeas.info_status; outer_loop 0 node.level flag node.bin_values']; #[outer_loop inner_loop level bnbstate binaryvalues] 
      #push!(S,region_infeas)
      break  
    end
  end
  return infeas_prob, region_infeas
end
## _______________________________________________________
function primal_solution(node,dualprob,partition,problem_type_QP)
    sol_part = ParamSol[]; 
    if problem_type_QP
      #_________________ MIQP sol________________
      # primal solution from dual solution --> MIQP
      n = size(node.f); nth = size(node.W,2); 
      unb_sol_F = -node.H\node.f_theta;
      unb_sol_G = -node.H\node.f;
      R = cholesky(node.H);
      Fv = (R.U')\node.f_theta;
      Gv = (R.U')\node.f;
      vQ = 0.5*(Fv')*Fv; 
      vR = Gv'*Fv; 
      vS = 0.5*(Gv')*Gv;
  
      for part in partition #p in 1:length(partition) #part = partition[p]
        # first check if optimal or unbounded or infeasible
        if string(part.state) != "UNBOUNDED" # "OPTIMAL" #Int(partition[p].state == 1/2?)   
          AS = part.AS;    
          if !isempty(AS) #optimal sol
            F_lam = part.Lam[1:end-1,:]
            F_lam = F_lam'[:,:]; #transpose!
            G_lam =part.Lam[end,:];
            M = dualprob.M[AS,:]; M = M'[:,:];
            Fu = M*F_lam; 
            Gu = M*G_lam; 
            uQ = 0.5*(Fu')*Fu;
            uR = Gu'[:,:]*Fu;
            uS = 0.5*(Gu')*Gu;
  
            ##primal function: F_theta th + G_theta
            F_theta =  -R.U\(M*F_lam) + unb_sol_F;
            G_theta = -R.U\(M*G_lam) + unb_sol_G;
          else #infeasile 
            uQ = 0;
            uR = 0;
            uS = Inf; 
            #primal function: F_theta th + G_theta
            F_theta =  unb_sol_F;
            G_theta = unb_sol_G;
          end
          #Objective function: th' A_theta th + B_theta th+ C_theta
          A_theta = uQ.-vQ; 
          B_theta = uR.-vR; 
          C_theta = uS-vS; 
          #sum up extra cost from fixing a binary variable and reducing the size of the problem
          if opts.add_extra_cost
            #B_theta = B_theta + node.extra_cost_param';
            C_theta = C_theta + node.extra_cost; # ; # 
          end
          sol_one_part = ParamSol(F_theta,G_theta,A_theta,B_theta,C_theta)
        else #if string(partition[p].state) == "UNBOUNDED" #Int(partition[p].state == 5)      
            #unbounded
            sol_one_part = ParamSol(zeros(n,nth),NaN*ones(n),zeros(nth,nth),zeros(1,nth),-Inf)
        end  #if optimal or unbounded or infeasible
        push!(sol_part,sol_one_part)
      end #for p
      #______________________
    else
        #_________________ MILP sol________________
        # primal solution from dual solution --> MILP
        n = length(node.f); nth = size(node.W,2);
        A_theta = zeros(nth,nth); #MILP: J = B_theta * th + C_theta
  
        for part in partition #p in 1:length(partition) 
          # first check if optimal or unbounded or infeasible
           if string(part.state) == "OPTIMAL" #"UNBOUNDED" # 
              AS = part.AS;    
              F_theta = node.A[AS,:]\node.W[AS,:];
              G_theta = node.A[AS,:]\node.b[AS];
              B_theta = -part.Lam*node.W[AS,:]; #[end,:]
              C_theta1 = part.Lam*node.b[AS];
              C_theta = -C_theta1[1];
              #sum up extra cost from fixing a binary variable and reducing the size of the problem
              if (false) #opts.add_extra_cost # TODO: check #( false) # added below
                #B_theta = B_theta + node.extra_cost_param'; # TODO : enable
                C_theta = C_theta + node.extra_cost; # ; # 
              end
              sol_one_part = ParamSol(F_theta,G_theta,A_theta,B_theta,C_theta)
          
            elseif string(part.state) == "INFEASIBLE"      
            #infeasible
              sol_one_part = ParamSol(zeros(n,nth),NaN*ones(n),zeros(nth,nth),zeros(1,nth),Inf)
            else #if string(partition[p].state) == "UNBOUNDED"  
            #unbounded
            sol_one_part = ParamSol(zeros(n,nth),NaN*ones(n),zeros(nth,nth),zeros(1,nth),-Inf)
          end  #if optimal or unbounded or infeasible
          push!(sol_part,sol_one_part)
        end #for p
    end #if MIQP
  
    return sol_part
end
## ___________________________________________________________________
## ____________________________________________________________________
function update_sol_to_entire(sol_part, mainprob, n, nth, bin_ind, bin_values,problem_type_QP)
    sol_updated = ParamSol[]; #bin_values = node.bin_values; sol_part = deepcopy(sol_part_ld) #for debug
    for sol in sol_part #sol = sol_part[1] #for debug
        F_th = sol.F_theta;
        G_th = sol.G_theta;
        bin_fixed = findall(bin_values .> -1) #bin_free = findall(bin_values .== -1)
        ints = bin_ind[bin_fixed]
        if (!isempty(ints))
          #for j in eachindex(ints) # = 1:length(ints)
          #  F_theta = [F_th[1:ints[j]-1,:]; zeros(1,nth); F_th[ints[j]:end,:]];
          #  G_theta = [G_th[1:ints[j]-1,:]; bin_values[j]; G_th[ints[j]:end,:]];
          #end
          #TODO : fix below for when branching is not in order
          F_theta = [F_th[1:ints[1]-1,:]; zeros(length(ints),nth); F_th[ints[1]:end,:]];
          G_theta = vec([G_th[1:ints[1]-1,:]; bin_values[bin_fixed]; G_th[ints[1]:end,:]]);
  
          if !isinf(sol.C_theta)
            if problem_type_QP#for MIQP #opts.
              A_theta = 0.5*F_theta'*mainprob.H*F_theta + mainprob.f_theta'* F_theta; 
              B_theta = (G_theta'*mainprob.H + mainprob.f')*F_theta + G_theta'* mainprob.f_theta;
              C_theta = 0.5*G_theta'*mainprob.H*G_theta+ mainprob.f'*G_theta;
            else
              A_theta = zeros(nth,nth);
              B_theta = mainprob.f'*F_theta;
              C_theta = (mainprob.f'*G_theta)[1];
            end
          else
            A_theta = sol.A_theta; B_theta = sol.B_theta; C_theta = sol.C_theta;
          end
          sol_reg_up = ParamSol(F_theta,G_theta,A_theta,B_theta,C_theta)
          push!(sol_updated,sol_reg_up)
        else
          push!(sol_updated,sol)
        end
    end
    return sol_updated
end
## ______________________________________________________________________
  function isdominated(reg,sol,sol_bar,opts)    #::Region
    #reg =deepcopy(partition[j]); sol = deepcopy(sol_part[j]); sol_bar = deepcopy(new_region.sol_bar);    #for debug     
    if (sol.C_theta == Inf)  #infeasible
      #__________invoke infeasibility cut____________
      return true, InfeasCut, [], []  
    elseif (sol.C_theta == -Inf)  #unbounded
        #__________dismiss infeasibility cut____________
        return false, Unbounded, [], []    
    elseif isinf(sol_bar[end].C_theta) 
      #__________dissmiss dominance cut____________
      return false, BNBState[], [], [] #Jbar is inf --> not cut
    else 
      dom_cut = false; flag = BNBState[];
      reg_J_better = []; reg_J_worse = []; nth = size(sol_bar[1].B_theta,2);
      #__________compare J and Jbar ____________ 
      for i in eachindex(sol_bar) # = 1: length(sol_bar) #for conservative miqp
        # Check if J-Jbar ≥ 0 ∀θ in reg      
        #Ath = reg.Ath'; #bth = reg.bth; # in region: Ath θ <= bth
        # Delta J = J - Jbar = θ A_delta θ + B_delta' θ + C_delta  
        if (opts.subopt_epsilon) 
          if opts.tol_abs_param  
            A_eps = opts.tol_domcut_abs_param[:,1:end-1]; B_eps = opts.tol_domcut_abs_param[:,end]; 
          else
            A_eps = zeros(nth,nth); B_eps = zeros(nth);
          end
          #___________ Suboptimality ____________
          A_delta = zeros(nth,nth);
          B_delta = vec(sol.B_theta - sol_bar[i].B_theta); #(1+ opts.tol_domcut_rel).*
          C_delta= (sol.C_theta + opts.tol_domcut_abs) - sol_bar[i].C_theta; #(1+ opts.tol_domcut_rel).*
          eq_A = (norm(A_delta) <= opts.eps_zero ) #) || (norm(sol.A_theta - sol_bar[i].A_theta) <= opts.eps_zero); #optimality or subopt.
          eq_B = (norm(B_delta) <= opts.eps_zero ) #|| (norm(sol.B_theta - sol_bar[i].B_theta) <= opts.eps_zero);
          non_param_diff = (eq_B) && (eq_A); #Diff J is non parametric
          #Better_J =  (non_param_diff)  && (C_delta <= -opts.eps_zero) #|| (sol.C_theta - sol_bar[i].C_theta <= -opts.eps_zero)) #obvoiusly better
        else
          #___________ Optimality ____________
          A1 = sol.A_theta - sol_bar[i].A_theta; A_delta = (A1 + A1')/2;
          B_delta = vec(sol.B_theta - sol_bar[i].B_theta);
          C_delta= sol.C_theta - sol_bar[i].C_theta;
          non_param_diff = (norm(B_delta) <= opts.eps_zero ) && (norm(A_delta) <= opts.eps_zero); #Diff J is non parametric
          #Better_J =  (non_param_diff)  && (C_delta <= -opts.eps_zero) #obvoiusly better
        end
        #______________________ 
        Better_J =  (non_param_diff)  && (C_delta <= -opts.eps_zero) #obvoiusly better
        equal_worse_J =  (non_param_diff)  && (C_delta >= 0 ) #-opts.eps_zero) #equal or worse 
        
        if (equal_worse_J)
          #__________invoke dominance cut____________
          flag = DomCut; dom_cut = true; #holds in entire region # J >= Jbar --> cut
          reg_J_better = []; reg_J_worse = []; #reg; #only need it in partially cut    
        elseif Better_J
          #__________dismiss infeasibility cut____________
          dom_cut = false; flag = BNBState[];
          reg_J_better = []; ##only need it in partially cut 
          reg_J_worse = [];
        else
          #__________invoke partially dominance cut?____________
          dom_cut, flag, reg_J_better, reg_J_worse  = compare_value_functions(reg,A_delta,B_delta,C_delta, opts)
        end #if equal J 
        if (flag == DomCut); break; end  #J is worse
      end # for i
      return dom_cut, flag, reg_J_better, reg_J_worse 
    end # if feasible 
end #function
  
## ________________________________________________________________________
function compare_value_functions(reg,A_delta,B_delta,C_delta, opts) 
    dom_cut= false; flag = BNBState[]; 
    nth = size(reg.Ath,1);
  
    if opts.problem_type_QP    
      # compare QPs
        QP_diff = MIQP(A_delta,B_delta,reg.Ath,reg.bth,[],[],[]);
        # find min of J-Jbar
        flag_QP_diff, ~ , J_QP_diff_min, ~ = QLPSolve(QP_diff, opts);
        J_QP_diff_min = J_QP_diff_min + C_delta; 
  
        if (string(flag_QP_diff) == "OPTIMAL" &&  (J_QP_diff_min >= 0 ))
            #__________invoke dominance cut____________
            # no better solution --> cut entire region    
            dom_cut = true; flag = DomCut; 
            return dom_cut, flag, [], [] #reg
        else
            #__________dismiss dominance cut____________      
            return dom_cut, flag, [], [] #reg, [] # better solution --> maybe partly!
        end # if no QP_diff solution
      #____________________________________________ 
      #_______________ #compare LPs _______________  
    else   
      #_______________________________
      if !opts.LPcompare_call_polyhedra  
        # form regions with better and worse J
  
        reg_J_better = (Ath = [reg.Ath B_delta], bth = [reg.bth; -C_delta - opts.eps_zero]) # delta J <= 0  #J_diff = B_diff*theta+C_diff <= 0
        reg_J_worse = (Ath = [reg.Ath -B_delta], bth =  [reg.bth; C_delta]) #+ opts.eps_zero # delta J >= 0 # J_diff = B_diff*theta+C_diff >= 0  #';   
        # solve an LP to find if J improves in part of a region
        #opts.chosen_solver = "Gurobi"; #"GLPK";  #opts.problem_type_QP = false; #for LP
        prob_J_better = MIQP(zeros(nth,nth), ones(nth), reg_J_better.Ath', reg_J_better.bth, [],[],[]);  #dummy problem   
        flag_J_better, ~ , J_better, ~ = QLPSolve(prob_J_better, opts); #th_diff_min
        #(x_better,lam_better,AS_better,J_better,iter_opt_better,~) = DAQP.daqp_jl(prob_J_better,Int64[])
        
        if (string(flag_J_better) != "OPTIMAL")  #infeasible
          #__________invoke dominance cut____________
          # no better solution --> cut entire region
          dom_cut = true; flag = DomCut;
          return dom_cut, flag, [], [] #reg
        else
          # better solution
          prob_J_worse = MIQP(zeros(nth,nth),ones(nth),reg_J_worse.Ath',reg_J_worse.bth,[],[],[]);     
          flag_J_worse, ~ , J_worse, ~ = QLPSolve(prob_J_worse, opts); #th_diff_min
          #(x_worse,lam_worse,AS_worse,J_worse,iter_opt_worse,~) = DAQP.daqp_jl(prob_J_worse,Int64[])
          
          if (string(flag_J_worse) == "OPTIMAL") 
            #__________invoke partially dominance cut____________
            # partially better solution --> cut partial of region
            (opts.remove_redundancy) ? ((Arth,brth) = ASCertain.remove_redundant(reg_J_worse.Ath,reg_J_worse.bth); reg_J_worsee = (Ath = Arth, bth = brth) ) : (reg_J_worsee = reg_J_worse)
            dom_cut = false; flag = PartialDomCut;
            return dom_cut, flag, reg_J_better, reg_J_worsee
          else
            #__________dissmiss dominance cut____________
            # better solution in entire region
            dom_cut = false; flag = BNBState[];
            return dom_cut, flag, [], [] #reg, []
          end
        end # if no QP_diff solution 
        #____________________________________
      else
        ## ________________ Polyhedron_____________________
        #using Polyhedra, CDDLib
        # first check if region is feasible: # call polyhedron package to compare LPs 
        #part_feas = polyhedron(hrep(Array(reg.Ath'), reg.bth)); #, CDDLib.Library() 
        #if (npoints(part_feas) == 0) #&& (
        #  return true, InfeasCut, [], [] #region is empty
        #end
        #____ if region is feasible, check if decompose it or not
        reg_J_better = (Ath = [reg.Ath B_delta], bth = [reg.bth; -C_delta - opts.eps_zero]) # delta J <= 0  #J_diff = B_diff*theta+C_diff <= 0
        reg_J_worse = (Ath = [reg.Ath -B_delta], bth =  [reg.bth; C_delta]) # delta J >= 0 # J_diff = B_diff*theta+C_diff >= 0  #';       
        
        P_better = polyhedron(hrep(Array(reg_J_better.Ath'), reg_J_better.bth), CDDLib.Library())
        #plot(P_better, color="red", alpha=0.2) #for 2D
        #m = Polyhedra.Mesh(P_better); using MeshCat; vis = Visualizer(); setobject!(vis, m); IJuliaCell(vis)
        flag_J_better = !isempty(P_better) #(npoints(P_better) != 0)
        flag_J_better_full = (dim(P_better) == nth);
        #___________________
        if (flag_J_better && flag_J_better_full)
          #_______________dismiss cut condition yet_______________
          # better solution
          # see if it is better in the whole region
          P_worse = polyhedron(hrep(Array(reg_J_worse.Ath'), reg_J_worse.bth), CDDLib.Library()) # B_delta * th <= -C_delta + eps : -B_delta * th <= C_delta - eps
          flag_J_worse = !isempty(P_worse) #(npoints(P_worse) != 0)
          flag_J_worse_full = (dim(P_worse) == nth);
          #plot!(P_worse, color="green", alpha=0.2) #for 2D
          #m = Polyhedra.Mesh(P_worse); using MeshCat; vis = Visualizer(); setobject!(vis, m); IJuliaCell(vis)
          #partial_cut = flag_J_worse && flag_J_worse_full
          if (flag_J_worse && flag_J_worse_full)
             # partially better solution --> cut partial of region
            if (opts.compute_minimal_rep)  #TODO: uncomment (false) #
              #remoive redundant constraints: find minimal repersentation           
              #______________
              removehredundancy!(P_better)
              p_better = MixedMatHRep(P_better)      
              A_b = p_better.A; #b_b = p_better.b; 
              reg_J_b = (Ath = Array(A_b'), bth =  p_better.b);
              #______________
              removehredundancy!(P_worse)
              p_worse = MixedMatHRep(P_worse)      
              A_w = p_worse.A; #b_w = p_worse.b; 
              reg_J_w = (Ath = Array(A_w'), bth =  p_worse.b);    
              #___________________
              #_________________________________________              
              dom_cut = false; flag = PartialDomCut;
              return dom_cut, flag, reg_J_b, reg_J_w
            else
              #_________________________________________              
              dom_cut = false; flag = PartialDomCut;
              return dom_cut, flag, reg_J_better, reg_J_worse
            end
          else
            #__________dissmiss dominance cut____________
            # better solution in entire region
            dom_cut = false; flag = BNBState[];
            return dom_cut, flag, [], [] #reg, []
          end
        #____________________________________________
        else
          #__________invoke dominance cut____________
          # no better solution --> cut entire region
          dom_cut = true; flag = DomCut;
          return dom_cut, flag, [], [] #reg
        end
        #_____________________________        
      end #if call Polyhedron package
    end # if QP
end
## ________________________________________________________________________
function isIntfeasMP(node, AS) 
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
## ________________________________________________________________________
function branchingScoreMP(sol, bin_ind_all, bin_val, branchrule, reg, opts)
    #sol = sol_part[j]; bin_ind_all = bin_ind; branchrule=opts.branching_rule; reg = partition[j]; bin_val = node.bin_values; #for debug
    indb = findall(bin_val .==-1); #remove binary variables that have already been fixed
    bin_ind = bin_ind_all[indb];
    nb = size(bin_ind,1); #n = size(sol.G_theta,1);  nth = size(reg.Ath,1);
    reg_branch = []; flag_br_once = true; check_chebyC = opts.MInfB_check_chebyC; reg_br = []; #outputs
    #reg_br = (Ath = [], bth = [], bin_indx = []); push!(reg_branch,reg_br);
    if branchrule == FirstBin
      branch_ind = 1; 
    elseif branchrule == LastBin
      branch_ind = nb;
    elseif branchrule == MinBin
      branch_ind = argmin(xbin); #For integer case: # (abs.(xbin-round.(xbin))); 
    elseif branchrule == MaxBin
      branch_ind = argmax(xbin); #For integer case: # (abs.(xbin-round.(xbin)));
    elseif branchrule == MostInfeasBin #closest to 0.5
      lib = DefaultLibrary{Float64}(()->Gurobi.Optimizer()); #lib = DefaultLibrary{Float64}(()->GLPK.Optimizer());
      F = sol.F_theta; G = sol.G_theta; # x = F \theta + G
      Fb = F[bin_ind,:]; Gb = G[bin_ind]; # xbin = x[bin_ind];   
      #branch_ind = argmin(abs.(xbin .- 0.5)); #For integer case: # (abs.(xbin-round.(xbin) .- 0.5));
      #____________
      #first check if the solution in binary variables is not-parametric
      nonparam_sol = iszero(Fb);
      if nonparam_sol
        branch_ind = argmax(0.5.-abs.(Gb.-0.5)); #the one closest to 0.5 #TODO:fix like below
        #____________
      else
        #find which binary variable is not binary
        notbin = ones(nb); eps_zero = 1e-6;#isbin = zeros(nb);
        for j = 1 : nb
          notbin[j] = !(all(abs.(Fb[j,:]) .<= eps_zero) && (abs(Gb[j]) <= eps_zero || abs(Gb[j]-1) <= eps_zero)) #if(all(node.A[i,:] == 0) && all(node.b[1:end-1,i] == 0) && node.b[end,i] <0)
        end
        notbin_ind = bin_ind[notbin .!=0]; nbb = length(notbin_ind);
        #______________
        if (nbb == 1) 
          #just one (or none!) relaxed binary variable left that is not binary
          branch_ind = findall(!iszero,notbin)[1]; #bin_ind[notbin_ind[1]];
        elseif (nbb == 0) 
          # all are binary --> similar to onilne choose first one
          branch_ind = 1; #all are binary, branch on the first (default)
          #_______________
          #=  #Look at the chebychevCenter and calculate distance to 0.5    
          if check_chebyC         
            #solver = GLPK.Optimizer; lib = DefaultLibrary{Float64}(solver); ##using GLPK;
            lib = DefaultLibrary{Float64}(()->GLPK.Optimizer());
            P = polyhedron(hrep(Array(reg.Ath'),reg.bth),lib);
            th_ch,_  = chebyshevcenter(P,proper=false); #r_ch
            xb_ch = Fb*th_ch + Gb; #Fb[notbin_ind,:]*th_ch + Gb[notbin_ind,:]; #x_ch = F*th_ch + G; xb_ch = x_ch[bin_ind];
            branch_ind = argmax(0.5.-abs.(xb_ch.-0.5)); #.-0.5); 
          else
            branch_ind = 1; #all are binary, branch on the first (default)
          end #
          =# #____________
        else
          #add polyhedral cuts and possibly partition
               #____________
          ##compare score function of binary variables to partition Theta to pick the one closest to 0.5 
          H_comp  = deepcopy(reg.Ath);  K_comp = deepcopy(reg.bth); 
          ##_______________
          ##=
          if nbb == 2
            ind1 = notbin_ind[1]; ind2 = notbin_ind[2];
            #_______________________
            #_____________________
            #first check if x2-0.5 >= 0 ? then consider two following cases
            reg_feas2p, _, _= isRegFeasible([H_comp -F[ind2,:]], [K_comp; G[ind2]-0.5], opts); #x2 >= 0.5? #H_m2p, K_m2p 
            reg_feas2n, _, _ = isRegFeasible([H_comp F[ind2,:]], [K_comp; -G[ind2]+0.5], opts); #x2 <= 0.5?
            reg_feas1p, _, _ = isRegFeasible([H_comp -F[ind1,:]], [K_comp; G[ind1]-0.5], opts); #x1 >= 0.5?
            reg_feas1n, _, _ = isRegFeasible([H_comp F[ind1,:]], [K_comp; -G[ind1]+0.5], opts); #x1 <= 0.5?
            #_____________________
            ##x1 >= 0.5, x2 >= 0.5
            if (reg_feas1p && reg_feas2p) 
              reg_feas_pp, H_pp, K_pp = isRegFeasible([H_comp -F[ind1,:] -F[ind2,:] (F[ind1,:] - F[ind2,:])], [K_comp; G[ind1]-0.5; G[ind2]-0.5; G[ind2]-G[ind1]], opts);
              (reg_feas_pp) && (reg_pp = (Ath = H_pp, bth = K_pp, bin_indx = findall(bin_ind .== ind1)[1]); push!(reg_branch, reg_pp) ); #); #no partitioning! # x1 is closer to 0.5, pick it)
              reg_feas_ppc, H_ppc, K_ppc = isRegFeasible([H_comp -F[ind1,:] -F[ind2,:]  -(F[ind1,:] - F[ind2,:])], [K_comp; G[ind1]-0.5; G[ind2]-0.5; -(G[ind2]-G[ind1])], opts)  
              (reg_feas_ppc) && (reg_ppc = (Ath = H_ppc, bth = K_ppc, bin_indx = findall(bin_ind .== ind2)[1]); push!(reg_branch, reg_ppc) ); #); #no partitioning! # x1 is closer to 0.5, pick it)
            end
            #_____________________
            ##x1 <= 0.5, x2 >= 0.5 
            if (reg_feas1n && reg_feas2p) 
              reg_feas_np, H_np, K_np = isRegFeasible([H_comp F[ind1,:] -F[ind2,:] (-F[ind1,:] - F[ind2,:])], [K_comp; -G[ind1]+0.5; G[ind2]-0.5; (G[ind2]+G[ind1]-1)], opts);
              (reg_feas_np) && (reg_np = (Ath = H_np, bth = K_np, bin_indx = findall(bin_ind .== ind1)[1]); push!(reg_branch, reg_np) ); #); #no partitioning! # x1 is closer to 0.5, pick it)
              reg_feas_npc, H_npc, K_npc = isRegFeasible([H_comp F[ind1,:] -F[ind2,:]  -(-F[ind1,:] - F[ind2,:])], [K_comp; -G[ind1]+0.5; G[ind2]-0.5; -(G[ind2]+G[ind1]-1)], opts)  
              (reg_feas_npc) && (reg_npc = (Ath = H_npc, bth = K_npc, bin_indx = findall(bin_ind .== ind2)[1]); push!(reg_branch, reg_npc) ); #); #no partitioning! # x1 is closer to 0.5, pick it)
            end
            #_____________________
            ##x1 >= 0.5, x2 <= 0.5
            if (reg_feas1p && reg_feas2n) 
              reg_feas_pn, H_pn, K_pn = isRegFeasible([H_comp -F[ind1,:] F[ind2,:] (F[ind1,:] + F[ind2,:])], [K_comp; G[ind1]-0.5; -G[ind2]+0.5; (-G[ind2]-G[ind1]+1)], opts);
              (reg_feas_pn) && (reg_pn = (Ath = H_pn, bth = K_pn, bin_indx = findall(bin_ind .== ind1)[1]); push!(reg_branch, reg_pn) ); #); #no partitioning! # x1 is closer to 0.5, pick it)
              reg_feas_pnc, H_pnc, K_pnc = isRegFeasible([H_comp -F[ind1,:] F[ind2,:]  -(F[ind1,:] + F[ind2,:])], [K_comp; G[ind1]-0.5; -G[ind2]+0.5; -(-G[ind2]-G[ind1]+1)], opts)  
              (reg_feas_pnc) && (reg_pnc = (Ath = H_pnc, bth = K_pnc, bin_indx = findall(bin_ind .== ind2)[1]); push!(reg_branch, reg_pnc) ); #); #no partitioning! # x1 is closer to 0.5, pick it)
            end
            #_____________________
            ##x1 <= 0.5, x2 <= 0.5
            if (reg_feas1n && reg_feas2n) 
              reg_feas_nn, H_nn, K_nn = isRegFeasible([H_comp -F[ind1,:] F[ind2,:] (F[ind1,:] + F[ind2,:])], [K_comp; G[ind1]-0.5; -G[ind2]+0.5; (-G[ind2]-G[ind1]+1)], opts);
              (reg_feas_nn) && (reg_nn = (Ath = H_nn, bth = K_nn, bin_indx = findall(bin_ind .== ind1)[1]); push!(reg_branch, reg_nn) ); #); #no partitioning! # x1 is closer to 0.5, pick it)
              reg_feas_nnc, H_nnc, K_nnc = isRegFeasible([H_comp -F[ind1,:] F[ind2,:]  -(F[ind1,:] + F[ind2,:])], [K_comp; G[ind1]-0.5; -G[ind2]+0.5; -(-G[ind2]-G[ind1]+1)], opts)  
              (reg_feas_nnc) && (reg_nnc = (Ath = H_nnc, bth = K_nnc, bin_indx = findall(bin_ind .== ind2)[1]); push!(reg_branch, reg_nnc) ); #); #no partitioning! # x1 is closer to 0.5, pick it)
            end
            #_____________________
            if !(reg_feas1p || reg_feas2p || reg_feas1n || reg_feas2n)
              # Non of the constraints hoolds: Look at the chebychevCenter and calculate distance to 0.5
              P = polyhedron(hrep(Array(reg.Ath'),reg.bth),lib);
              th_ch,_  = chebyshevcenter(P,proper=false); #r_ch
              xb_ch = Fb*th_ch + Gb; #Fb[notbin_ind,:]*th_ch + Gb[notbin_ind,:]; #x_ch = F*th_ch + G; xb_ch = x_ch[bin_ind];
              branch_ind = argmax(0.5.-abs.(xb_ch.-0.5)); #.-0.5);
            end
          else # 
            #________________
            for i = 1 : nbb
              #H = reg.Ath; K = reg.bth; #H_aug = deepcopy(H); K_aug = deepcopy(K); 
              Fi, Gi = score_MInfB(F[notbin_ind,:],G[notbin_ind],nbb,i);
              H_aug = [reg.Ath Fi]; K_aug = [reg.bth; Gi]; #H_aug = [reg.Ath Fi[:, (4*i+1):end]]; K_aug = [reg.bth; Gi[(4*i+1):end]];
              reg_feas, H_r, K_r = isRegFeasible(H_aug, K_aug, opts)
              if reg_feas
                reg_i = (Ath = H_r, bth = K_r, bin_indx = findall(bin_ind .== i));
                push!(reg_branch, reg_i); #push!(reg_br, P_r); #, P_r
                H_comp = [H_comp -Fi]; K_aug = [K_comp; -Gi .- eps_zero]; #complement of the region
              end
              (i==nbb)&&(break)
            end
          end #if
          #____________
          #____________
          if length(reg_branch) > 1
            flag_br_once = false;
            branch_ind = length(reg_branch); #-1 #number of variables to branch on
            #___________________
            add_reg_comp = false; # true; # #check to see if some parts are not considered
            if add_reg_comp
              #Pu = deepcopy(reg_br[1]); for k = 2 : length(reg_branch); Pu = convexhull(Pu, reg_br[k]); end
              reg_feas_comp, H_compp, K_compp = isRegFeasible(H_comp, K_comp, opts); #check to see if some parts are not considered
              #________
              if reg_feas_comp
                #Some parts of region is left
                P_comp = polyhedron(hrep(Array(H_compp'),K_compp),lib);
                #removehredundancy!(P_comp1); P_comp = MixedMatHRep(P_comp1);   #find minimal repersentation  #already done in isReg function
                th_chc,_  = chebyshevcenter(P_comp,proper=false); #r_ch
                xb_chc = Fb*th_chc + Gb; #Fb[notbin_ind,:]*th_ch + Gb[notbin_ind,:]; #x_ch = F*th_ch + G; xb_ch = x_ch[bin_ind];
                br_indc = argmax(0.5.-abs.(xb_chc.-0.5)); #.-0.5);
                reg_c = (Ath = H_compp, bth = K_compp, bin_indx = br_indc);
                push!(reg_branch, reg_c); #push!(reg_br, P_r); #, P_r
              end
            end
            #_________
          elseif length(reg_branch) == 1
            branch_ind = reg_branch[end].bin_indx;
          else
            #___________ All halfplanes are redundanct
            if check_chebyC
              #Look at the chebychevCenter and calculate distance to 0.5
              P = polyhedron(hrep(Array(reg.Ath'),reg.bth),lib);
              th_ch,_  = chebyshevcenter(P,proper=false); #r_ch
              xb_ch = Fb*th_ch + Gb; #Fb[notbin_ind,:]*th_ch + Gb[notbin_ind,:]; #x_ch = F*th_ch + G; xb_ch = x_ch[bin_ind];
              branch_ind = argmax(0.5.-abs.(xb_ch.-0.5)); #.-0.5); 
            else
              branch_ind = findall(notbin .== 1)[1]; #TODO: fix: how to define it? maybe using chebycenter?
            end #cheby
          end
        end #if just one not binary
      end #if not parametric
    else
      branch_ind = 1;  #default: in order
    end
    return branch_ind, reg_branch, flag_br_once
end 
## ________________________________________________________________________
function score_MInfB(F,G,nb,i) #Most infeasible branching score function
Fi = zeros(nth, 4*nb); Gi = zeros(4*nb); #k= 0;#(nb-1)        
    for k = 1 : nb # while (k < nb);  k+=1;
      (k==i)&&(continue);
      Fik = [F[i,:] - F[k,:]  F[i,:] + F[k,:]  -F[i,:] - F[k,:]  -F[i,:] + F[k,:]];  #when finding max s
      Gik = [ G[k] - G[i]; -G[k] - G[i]+1; G[k] + G[i]-1; -G[k] + G[i]];
      inds = (k-1)*4 +1; 
      Fi[:,inds:k*4] = Fik; Gi[inds:k*4] = Gik;
      #push!(Fi, Fik); push!(Gi, Gik); #Fi = [Fi Fik]; Gi = [Gi; Gik]; #(k == nb) & (return Fi, Gi) 
    end
    return Fi, Gi
end
## ________________________________________________________________________
function isRegFeasible(H,K,opts)
    #H = [reg_r.Ath Array((node.A*reg_r.sol.F_theta-node.W)')]; K= [reg_r.bth; node.b-node.A*reg_r.sol.G_theta]
    #H_aug = deepcopy(H); K_aug= deepcopy(K); H= [H -F_h[indb,:]]; K = [K; G_h[indb]-0.5-eps_zero]
    #P = polyhedron(hrep(Array(H'),K)); plot(Pn, color="black", alpha=0.5) # npoints(P); fulldim_P = (fulldim(P) == nth); 
    ##for debug
    lib = DefaultLibrary{Float64}(()->Gurobi.Optimizer()); #lib = DefaultLibrary{Float64}(()->CDD.Optimizer())
    #lib = DefaultLibrary{Float64}(()->GLPK.Optimizer()); #OutputFlag = 0
    #with_optimizer(Gurobi.Optimizer; OutputFlag=0) #, OutputFlag=0
    reg_feas = false; H_r = H; K_r = K; #P = [];
    nth = size(H,1);
    #_______________________________
    if !opts.LPcompare_call_polyhedra  
      # solve an LP to find if region is feasible
      #opts.chosen_solver = "Gurobi"; #"GLPK";  #opts.problem_type_QP = false; #for LP
      prob_feas = MIQP(zeros(nth,nth), ones(nth), Array(H'), K, [],[],[]);  #dummy problem   
      flag_feas, _ , _, _ = QLPSolve(prob_feas, opts); #th_diff_min
      #(x_better,lam_better,AS_better,J_better,iter_opt_better,~) = DAQP.daqp_jl(prob_feas,Int64[])
      
      if (string(flag_feas) != "OPTIMAL")  #infeasible
        reg_feas = false; H_r = []; K_r = []; #region is infeasible, clean the space
      else
        reg_feas = true;     
        H_r, K_r = ASCertain.remove_redundant(H,K); #USe certification tools to remove redundancy
      end 
      #____________________________________
    else
      ## ________________ Polyhedron_____________________
      #using Polyhedra, CDDLib
      #Pn = polyhedron(hrep(Array(H'), K), lib) # plot(Pn, color="black", alpha=0.5)
      Pn = polyhedron(hrep(Array(H'), K), CDDLib.Library())
      #flag_feas = !isempty(Pn) #(npoints(Pn) != 0)
      #reg_feas_full = (dim(Pn) == nth);npoints(Pn);
      #removehredundancy!(Pn); Pnn = MixedMatHRep(p); dim(Pnn) #Arth = Pn.A; brth = Pn.b
      #plot(Pn, color="blue", alpha=0.2) # m = Polyhedra.Mesh(Pn); using MeshCat; vis = Visualizer(); setobject!(vis, m); IJuliaCell(vis)
      fulldim_Pn = true;
      try
        fulldim_Pn =  (npoints(Pn) != 0) && (dim(Pn) == nth) && (fulldim(Pn) == nth)
      catch
        fulldim_Pn = false;
      end #
      #fulldim_Pn =  (!isempty(Pn)) && (dim(Pn) == nth) && (fulldim(Pn) == nth) #(npoints(Pn) != 0) && 
  
      #___________________
      if fulldim_Pn #(flag_feas && reg_feas_full)
        reg_feas = true; 
          if (opts.compute_minimal_rep) 
            #remoive redundant constraints       
            removehredundancy!(Pn)
            P = MixedMatHRep(Pn);   #find minimal repersentation       
            H_r = P.A; K_r = P.b; H_r = Array(H_r');
            #______________
          end    
      else
        reg_feas = false; H_r = []; K_r = []; #region is infeasible, clean the space
      end    
      #__________
    end #if call Polyhedron package
    #_____________________________   
    return reg_feas, H_r, K_r #, P   
end 
## _________________________________________________________
function ordering_costMP(level,sorting,J)
    if sorting == DepthFirst
      priority_order = 1 / (level+1); 
    elseif sorting == BreadthFirst
      priority_order = level+ 1; #1.;
    elseif sorting == BestFirst
      priority_order = [J.B_theta J.C_theta]; #parametric priority: lower bound
    else
      priority_order = 1 / (level+1); #default:"depth_first"
    end
    return priority_order
end
## __________________________________________________________________________
function branchingMP(prob, br_ind::Int64, priority_order, AS, lam, opts) #::MIQP,::Vector{Int64}
    #AS = partition[j].AS; lam =partition[j].Lam
    n = length(prob.f); #fieldnames(typeof(prob)) #see the field of the struct
    branch_ind = prob.bin_ind[br_ind]
    indices = [1:n;];
    indices = setdiff(indices, branch_ind);
    prob1 = deepcopy(prob);
    prob0 = deepcopy(prob);
  
    # reduce new child nodes
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
    const_ind = setdiff([1:m;],rm_bin_indices); 
    A1 = prob.A[:,indices]; # remove the column for removed xb
    A = A1[const_ind,:];    # remove rows for binary constraints of xb
    prob0.A = A; prob1.A = A; 
  
    b = prob.b[const_ind];
    prob0.b = b; prob1.b = b - prob.A[const_ind,branch_ind];
  
    if isdefined(prob,:W)
      W = prob.W[const_ind,:];
      prob0.W = W; prob1.W = W;
    end
    #prob1.b + prob1.W*th_min #for debug, compare with online
  
    prob0.extra_cost = prob.extra_cost; 
    prob1.extra_cost = prob1.extra_cost + .5*H22 + f2; #add constant terms in value function
    #prob1.extra_cost = prob1.extra_cost + f2; #add constant terms in value function
    #_________Add just current extra cost____________
    #prob0.extra_cost = 0;
    #prob1.extra_cost = .5*H22 + f2; #add current constant terms in value function
  
    if isdefined(prob,:extra_cost_param)
      extra_cost_param = prob.f_theta[branch_ind,:];#add parametric terms in value function
      #prob0.extra_cost_param = prob.extra_cost_param; 
      #prob1.extra_cost_param = prob1.extra_cost_param + extra_cost_param; 
      #_________Add just current extra cost____________
      #println(extra_cost_param); #readline();
      #prob0.extra_cost_param = []; prob1.extra_cost_param = []; 
      prob0.extra_cost_param = zeros(size(extra_cost_param,1),1);#
      prob1.extra_cost_param = zeros(size(extra_cost_param,1),1)+ extra_cost_param; 
    end
  
    bin_ind = [prob.bin_ind[1: br_ind-1]; prob.bin_ind[br_ind+1:end].-1]; #remove fixed binary var
    prob0.bin_ind = bin_ind; prob1.bin_ind = bin_ind; 
  
    bin_indices = prob.bin_indices[1:end .!= br_ind, :];
    #bin_indices[:,1] = bin_indices[:,1].-1; bin_indices[:,2] = bin_indices[:,2].-2; #Only for FirstBin branching rule
    bin_indices[:,1] = [prob.bin_indices[1: br_ind-1, 1]; prob.bin_indices[br_ind+1:end, 1].-1]; #bin_indices[:,1].-1; 
    bin_indices[:,2] = [prob.bin_indices[1: br_ind-1, 2].-1; prob.bin_indices[br_ind+1:end, 2].-2]; #bin_indices[:,2].-2;
    prob0.bin_indices = bin_indices; prob1.bin_indices = bin_indices; 
  
    prob0.bin_states[br_ind] = LOWER;
    prob1.bin_states[br_ind] = UPPER;
    #prob0.senses[prob.bin_indices[br_ind,1]]=EQUALITY;
    #prob1.senses[prob.bin_indices[br_ind,2]]=EQUALITY;
    #prob1.sense[prob.bin_indices[br_ind,1]]|=(DAQP.ACTIVE+DAQP.IMMUTABLE); #correct format
    #prob0.sense[prob.bin_indices[br_ind,2]]|=(DAQP.ACTIVE+DAQP.IMMUTABLE);# commented since prob is reduced
    n_const = size(prob0.b,1);
    prob0.bounds_table = collect(1:n_const);
    prob0.senses = zeros(Cint,n_const);
    prob1.bounds_table = collect(1:n_const);
    prob1.senses = zeros(Cint,n_const);
    
    if (opts.warm_start && !isempty(AS))
      AS_updated0, n0_feas, AS_updated1, n1_feas = basis_recoveryMP(prob,prob0, AS,br_ind,lam) #use (half of) simplex method to update AS
      prob0.AS = AS_updated0; prob0.feasible = n0_feas;
      prob1.AS = AS_updated1; prob1.feasible = n1_feas;
    else
      prob0.AS = copy(AS); 
      prob1.AS = copy(AS);
    end
  
    bin_free = findall(prob.bin_values .== -1);#bin_fixed = findall(prob.bin_values .> -1);
    branch = bin_free[br_ind];
    prob0.bin_values[branch] = 0;
    prob1.bin_values[branch] = 1;
  
    prob0.level = prob0.level + 1;  
    prob1.level = prob1.level + 1;
  
    if (opts.sorting_method == BestFirst) # BF: parametric order
      prob0.priority_BF = priority_order; prob1.priority_BF = priority_order; # BF
      prob0.priority += 1; prob1.priority += 1; # default: BrF? 
    else
      #typeof(priority_order) == Float64  #DF and BrF
      prob0.priority = priority_order; prob1.priority = priority_order;
    end
   
  return prob0,prob1
end
## _________________________________________________________
function basis_recoveryMP(node,node0, AS_parent,br_ind,lam)
    #br_ind = node.bin_ind[br_indx];
    m = node.bin_indices[1,1]-1; #ineq constraints without binary constraints
    nb = length(node.bin_ind); n = size(node0.f,1);
    AS_ineq = AS_parent[AS_parent .<= m];  #ActiveSet of original inequality constraints  
    bin_ind0 = m .< AS_parent .< (m+nb+1); #AS_parent[AS_parent .> m]        
    AS_b_lb = AS_parent[bin_ind0];         #AS of lower bound of left binary constraints
    AS_b_ub = AS_parent[AS_parent .> m+nb];#AS of lower bound of left binary constraints
    #remove index of branchvar and branchvar+nb from AS
    AS_bin_mod = [AS_b_lb[AS_b_lb .< m+br_ind]; AS_b_lb[AS_b_lb .> m+br_ind].-1; AS_b_ub[AS_b_ub .< m+nb+br_ind].-1; AS_b_ub[AS_b_ub .> m+nb+br_ind].-2];
    AS_updated = [AS_ineq; AS_bin_mod];
    #AS_bin_mod = [AS_parent[m .< AS_parent .< (m+nb+1)].-1; AS_parent[(m+nb+1) .<= AS_parent].-2];
    #AS_updated = [AS_parent[AS_parent .<= m]; AS_bin_mod]; 
  
    if (length(AS_updated) > n)
      #apply (half of) simplex method to remove one index from AS (in simplex: lenght(AS) = n)
      n0_feas = true; n1_feas = true; 
      fix_id_0,fix_id_1 = node.bin_indices[br_ind,:];
      #___________________________
      AS0 = fix_warm_AS(node,AS_parent,lam,fix_id_0);  
      if !isempty(AS0)
        #AS_updated0 = [AS0[AS0 .<= m]; AS_bin_mod];
        #AS_updated0 = [AS0[AS0 .<= m]; AS0[m .< AS0 .< (m+nb+1)].-1; AS0[(m+nb+1) .<= AS0].-2];   
        AS_b_lb0 = AS0[m .< AS0 .< (m+nb+1)];         #AS of lower bound of left binary constraints
        AS_b_ub0 = AS0[AS0 .> m+nb];#AS of lower bound of left binary constraints
        AS_updated0 = [AS0[AS0 .<= m]; AS_b_lb0[AS_b_lb0 .< (m+br_ind)]; AS_b_lb0[AS_b_lb0 .> (m+br_ind)].-1; AS_b_ub0[AS_b_ub0 .< m+nb+br_ind].-1; AS_b_ub0[AS_b_ub0 .> m+nb+br_ind].-2];
      else
        AS_updated0 = Int64[]; n0_feas = false
      end
      #___________________________
      AS1 = fix_warm_AS(node,AS_parent,lam,fix_id_1);
      if !isempty(AS1)
        AS_b_lb1 = AS1[m .< AS1 .< (m+nb+1)]; #(m+nb+br_ind)        #AS of lower bound of left binary constraints
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
## _________________________________________________________
function appendingMP(reg,childs,S,opts)
    (S == Nothing) && (S = BNBRegion[])
    if isempty(childs)
      #no new node to be pushed in T
      push!(S,reg)
      return S;
    end
    #____________________________________
    tot_gen_nod = deepcopy(reg.kappa.tot_gen_nodes); #for generated-node subopt cut
    decomp_nod = deepcopy(reg.kappa.decomp_node); #for generated-node subopt cut
    max_act_nod = opts.bound_active_nodes; #for active-node subopt cut
    #______________________________
    if isempty(reg.T) #T empty, no need to find a place, insert in the end
      #push!(reg.T,node1,node0);
      #for child_node in childs
      #  push!(reg.T,child_node)
      #end
      sortedT = childs; #[reg.T[1:end]; childs[1:end]];
      if (opts.subopt_Mcut && (length(sortedT) > max_act_nod))
        #reg.gen_nodes = tot_gen_nod + length(sortedT) - curr_act_nod; 
        ind_start = length(sortedT)- max_act_nod + 1;
        sortedT = [sortedT[ind_start:end];]; #sortedT = [sortedT[1:max_act_nod];];  #bests are in the end     
      #else; reg.gen_nodes = tot_gen_nod + length(childs);
      end
      reg.T = sortedT;   
      reg.kappa.decomp_node = decomp_nod + 1;
      reg.kappa.tot_gen_nodes = tot_gen_nod + length(childs);
      push!(S,reg)
      return S;
    end
    #____________________________________
    #ordering nodes: higher priority (smaller rho): put in the end of T (first popped)
    sorting = opts.sorting_method; 
    curr_act_nod = length(reg.T);
    #______________________________
    if sorting == DepthFirst
      #for child_node in childs
      #  push!(reg.T,child_node); end; #push!(reg.T,node0,node1); #DF (put in the end)
      sortedT = [reg.T[1:end]; childs[1:end]];
      if (opts.subopt_Mcut && (length(sortedT) > max_act_nod))
        #reg.gen_nodes = tot_gen_nod + length(sortedT) - curr_act_nod; 
        ind_start = length(sortedT)- max_act_nod + 1;
        sortedT = [sortedT[ind_start:end];]; #sortedT = [sortedT[1:max_act_nod];];  #bests are in the end     
      end
      reg.kappa.decomp_node = decomp_nod + 1;
      reg.kappa.tot_gen_nodes = tot_gen_nod + length(sortedT) - curr_act_nod; 
      reg.T = sortedT;   
      push!(S,reg)
      return S
      #____________________________________
    elseif sorting == BreadthFirst  
      #for child_node in childs
      #  pushfirst!(reg.T,child_node); end #pushfirst!(reg.T,node0,node1); #BrF (put in the beginning)
      sortedT = [childs[1:end];  reg.T[1:end]]; 
      if (opts.subopt_Mcut && (length(sortedT) > max_act_nod))     
        ind_start = length(sortedT)- max_act_nod + 1;
        sortedT = [sortedT[ind_start:end];]; #sortedT = [sortedT[1:max_act_nod];];       
      end
      reg.kappa.decomp_node = decomp_nod + 1;
      reg.kappa.tot_gen_nodes = tot_gen_nod + length(sortedT) - curr_act_nod; 
      #reg.gen_nodes = tot_gen_nod + length(sortedT) - curr_act_nod; 
      reg.T = sortedT;
      #reg.gen_nodes = reg.gen_nodes + length(childs);
      push!(S,reg)
      return S
    end   
    #______________________________________
    if sorting == BestFirst
      node0 = childs[end];
      B_theta_currpr = node0.priority_BF[:,1:end-1]; #current parametric nodes priority
      C_theta_currpr = node0.priority_BF[1,end];
  
      if isinf(C_theta_currpr)
        #pushfirst!(reg.T,node1,node0); #worst priority: put in the beginning (pop from end)
        sortedT = [childs[1:end];  reg.T[1:end]];
        if (opts.subopt_Mcut && (length(sortedT) > max_act_nod))
          #reg.gen_nodes = tot_gen_nod + length(sortedT) - curr_act_nod; 
          ind_start = length(sortedT)- max_act_nod + 1;
          sortedT = [sortedT[ind_start:end];]; 
        end
        reg.kappa.decomp_node = decomp_nod + 1;
        reg.kappa.tot_gen_nodes = tot_gen_nod + length(sortedT) - curr_act_nod; 
        reg.T = sortedT;
        push!(S,reg)
        return S;
      end 
      #________________________________
      nth = size(node0.W,2); sortedT = NodeMP[];     
      ordering_done = false; # rho[i] - curr_rho >= 0 --> the place to insert
      for i = length(reg.T):-1:1 
        if isinf(reg.T[i].priority_BF[1,end]) #found a worst one
          sortedT = [reg.T[1:i]; childs[1:end]; reg.T[i+1:end]]
          if (opts.subopt_Mcut && (length(sortedT) > max_act_nod))
            ind_start = length(sortedT)- max_act_nod + 1;
            sortedT = [sortedT[ind_start:end];];   
          end
          reg.kappa.decomp_node = decomp_nod + 1;
          reg.kappa.tot_gen_nodes = tot_gen_nod + length(sortedT) - curr_act_nod; 
          reg.T = sortedT; ordering_done = true;
          push!(S,reg)
          return S;
          #____________________________
        else 
          B_delta = vec(B_theta_currpr - reg.T[i].priority_BF[:,1:end-1]);
          C_delta= C_theta_currpr - reg.T[i].priority_BF[1,end]; 
          #equal_Js = (norm(B_delta) <= opts.eps_zero) && (abs(C_delta) <= opts.eps_zero) #&& (norm(A_delta) <= opts.eps_zero)
          higher_prior_1 = (norm(B_delta) <= opts.eps_zero) && (C_delta <= opts.eps_zero) #&& (norm(A_delta) <= opts.eps_zero)
          if higher_prior_1 # (equal_Js || Better_J) #equal: best-depth, put it here first
            #________ put entire region here ____________
            # priorities are equal or better, put here
            sortedT = [reg.T[1:i]; childs[1:end]; reg.T[i+1:end]]
            if (opts.subopt_Mcut && (length(sortedT) > max_act_nod))
              #reg.gen_nodes = tot_gen_nod + length(sortedT) - curr_act_nod; 
              ind_start = length(sortedT)- max_act_nod + 1;
              sortedT = [sortedT[ind_start:end];]; #sortedT = [sortedT[1:max_act_nod];];  
            #else; reg.gen_nodes = tot_gen_nod + length(childs);
            end
            reg.kappa.decomp_node = decomp_nod + 1;
            reg.kappa.tot_gen_nodes = tot_gen_nod + length(sortedT) - curr_act_nod; 
            #reg.gen_nodes = tot_gen_nod + length(sortedT) - curr_act_nod; 
            reg.T = sortedT; ordering_done = true;
            push!(S,reg)
            return S;
            #____________________________
          else
            # compare priorities
            higher_prior, flag_partially, reg_J_better, reg_J_worse  = compare_value_functions(reg,zeros(nth,nth),-B_delta,-C_delta, opts)        
            if higher_prior #priority is higher (smaller) in entire region: rho[i] >= curr_rho 
              # push nodes entirely here
              sortedT = [reg.T[1:i]; childs[1:end]; reg.T[i+1:end]]
              if (opts.subopt_Mcut && (length(sortedT) > max_act_nod))
                ind_start = length(sortedT)- max_act_nod + 1;
                sortedT = [sortedT[ind_start:end];]; 
              end
              reg.kappa.decomp_node = decomp_nod + 1;
              reg.kappa.tot_gen_nodes = tot_gen_nod + length(sortedT) - curr_act_nod;  
              reg.T = sortedT; ordering_done = true;
              push!(S,reg)
              return S;
              #_________________________
            elseif flag_partially == PartialDomCut #partially better
              #put part of region where rho[i] >= curr_rho here
              reg_higher_pr = deepcopy(reg);
              reg_higher_pr.Ath = reg_J_worse.Ath; reg_higher_pr.bth = reg_J_worse.bth; #update region
              sortedT = [reg.T[1:i]; childs[1:end]; reg.T[i+1:end]]
              if (opts.subopt_Mcut && (length(sortedT) > max_act_nod))              
                ind_start = length(sortedT)- max_act_nod + 1;
                sortedT = [sortedT[ind_start:end];]; #sortedT = [sortedT[1:max_act_nod];];  
              end
              reg_higher_pr.kappa.decomp_node = decomp_nod + 1;
              reg_higher_pr.kappa.tot_gen_nodes = tot_gen_nod + length(sortedT) - curr_act_nod;  
              reg_higher_pr.T = sortedT
              push!(S,reg_higher_pr)
              #___________
              reg.Ath = reg_J_better.Ath; reg.bth = reg_J_better.bth; #update current region to the remaining region
              #sleep(3*60) #wait 3 min #wait_for_key("press any key to continue")
            end #if found place
          end #if equal prioritySS
        end #if inf priority
        #____________________________
        if (isempty(reg.bth)  || (norm(reg.bth) <= opts.eps_zero) ) #isempty(reg)) #
          ordering_done = true;
          return S
        end
      end # for
      #____________________________
      if !ordering_done
        #pushfirst!(reg.T,node1,node0); #worst priority: put in the beginning
        sortedT = [childs[1:end]; reg.T[1:end]];
        if (opts.subopt_Mcut && (length(sortedT) > max_act_nod))
          ind_start = length(sortedT)- max_act_nod + 1;
          sortedT = [sortedT[ind_start:end];]; #sortedT = [sortedT[1:max_act_nod];];     
        end
        reg.kappa.decomp_node = decomp_nod + 1;
        reg.kappa.tot_gen_nodes = tot_gen_nod + length(sortedT) - curr_act_nod; 
        #reg.gen_nodes = tot_gen_nod + length(sortedT) - curr_act_nod; 
        reg.T = sortedT;
        push!(S,reg)
        return S 
      end
      #________________________________
    else #room for other node selection strategies
      sortedT = NodeMP[];
      prior_nodes = zeros(length(reg.T));
      for i = length(reg.T):-1:1      
        prior_nodes[i] = reg.T[i].priority
      end
      ind1 = findall(node0.priority .- prior_nodes .<= opts.eps_zero);
      ind = isempty(ind1) ? 1 : ind1[end];
      sortedT = [reg.T[1:ind]; childs[1:end]; reg.T[ind+1:end]];
      reg.kappa.decomp_node = decomp_nod + 1;
      reg.kappa.tot_gen_nodes = tot_gen_nod + length(sortedT) - curr_act_nod; 
      reg.T = sortedT;
      push!(S,reg)
    return S;   
    end #if best_first
end # function
## ______________________________________________________

function update_xbar_cert(sol_part, mainprob, n, nth, cont_ind, bin_ind_curr, bin_ind, bin_values,problem_type_QP)
  sol_updated = ParamSol[]; #bin_values = node.bin_values; sol_part = deepcopy(sol_part_ld) #for debug
  if size(sol_part[1].G_theta,1) == n
    # root node, solution is updated
    sol_part[1].F_theta[bin_ind,:] .= 0;
    sol_part[1].G_theta[bin_ind] = round.(Int, sol_part[1].G_theta[bin_ind]);
    return sol_part
    #____________________
  else 
      for sol in sol_part #sol = sol_part[1] #for debug
          F_tot = zeros(n,nth); G_tot = zeros(n);
          F_tot[cont_ind,:] = sol.F_theta[cont_ind,:]; G_tot[cont_ind] = sol.G_theta[cont_ind]; 
          Gb = sol.G_theta[bin_ind_curr]; #Fb = sol.F_theta[bin_ind_curr,:]; 

          bin_fixed = findall(bin_values .> -1)
          bin_free = findall(bin_values .== -1)
          F_tot[bin_ind[bin_fixed], :].= 0; #zeros(size(bin_ind[bin_fixed],1), nth)
          G_tot[bin_ind[bin_fixed], :]= bin_values[bin_fixed]; #
          if !isempty(bin_ind[bin_free])
            F_tot[bin_ind[bin_free]] .= 0;
            G_tot[bin_ind[bin_free]] = round.(Int, Gb);
          end
          F_theta = deepcopy(F_tot); G_theta = deepcopy(G_tot);
          #_____________________________________________
          if !isinf(sol.C_theta)
            if problem_type_QP #for MIQP
                A_theta = 0.5*F_theta'*mainprob.H*F_theta + mainprob.f_theta'* F_theta; 
                B_theta = (G_theta'*mainprob.H + mainprob.f')*F_theta + G_theta'* mainprob.f_theta;
                C_theta = 0.5*G_theta'*mainprob.H*G_theta+ mainprob.f'*G_theta;
            else
                A_theta = zeros(nth,nth);
                B_theta = mainprob.f'*F_theta;
                C_theta = (mainprob.f'*G_theta)[1];
            end
          else
            A_theta = sol.A_theta; B_theta = sol.B_theta; C_theta = sol.C_theta;
          end
          sol_reg_up = ParamSol(F_theta,G_theta,A_theta,B_theta,C_theta)
          push!(sol_updated,sol_reg_up)
      end #for
    return sol_updated
  end #if
end 


## Extra functions
#=
## _________________________________________________________________
function update_sol(sol,n,nth,cont_ind,bin_ind,bin_ind_curr,bin_values)
    #(sol_part[j],n,nth,cont_ind,bin_ind,node.bin_ind,node.bin_values); bin_values = node.bin_values; bin_ind_curr = node.bin_ind; #for debug
    F_theta = zeros(n,nth);
    G_theta = zeros(n); 
    F_theta[cont_ind,:] = sol.F_theta[cont_ind,:];
    G_theta[cont_ind,:] = sol.G_theta[cont_ind];
    F_bin = sol.F_theta[bin_ind_curr,:];
    G_bin = sol.G_theta[bin_ind_curr];
  
    if size(sol.G_theta,1) == n
      G_theta[bin_ind] = round.(G_bin)
    else   
      bin_fixed = findall(bin_values .> -1)
      bin_free = findall(bin_values .== -1)
      #F_theta[bin_ind[bin_fixed],:] = zeros(length(bin_fixed),nth); #already done when initializing!
      G_theta[bin_ind[bin_fixed]] = bin_values[bin_fixed];    
      if !isempty(bin_ind[bin_free])
        #F_theta[bin_ind[bin_free],:] = zeros(length(bin_free),nth); #already done when initializing!
        G_theta[bin_ind[bin_free]] = round.(G_bin);
      end
      # ading extra cost to objective function here
    end
    updated_sol = ParamSol(F_theta,G_theta,sol.A_theta,sol.B_theta,sol.C_theta);
 return updated_sol
end
=#

  