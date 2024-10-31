## Generate random mpMIQP
function generate_mpMIQP(n,m,nth,nb,bnd,example) 
  ## Random Example
  if example == "rand"
    M = randn(n,n)
    H = M*M'
    f = randn(n)
    f_theta = randn(n,nth)
    H_theta = zeros(0,0);
    A = randn(m,n)   
    b = 2 .+ rand(m); 
    W = randn(m,nth); 
    bin_ind = collect(n-nb+1:n); 

  elseif example == "LP_specific" 
    #LP example from reference
    f = [-3; -2; 10; 5];
    A = [1 0 0 0; 0 1 0 0; 1 1 0 0; 1 2 0 0;1 0 -20 0; 0 1 0 -20;
         1 -1 0 0; 0 0 -1 -1;-1 0 0 0;0 -1 0 0];
    b = [10.; 10; 20; 12; 0; 0; -4; -1; 0; 0];
    W = [1 2 0; -1 1 0; 0 -1 0;1 0 -1; 0 0 0; 0 0 0; 0 0 1; 0 0 0;  0 0 0; 0 0 0];
    n = size(f,1); nth = size(W,2);
    H = zeros(n,n);
    f_theta = zeros(n,nth);
    H_theta = zeros(0,0);
    bin_ind = [n-1;n];
    nb = length(bin_ind);
    m = size(A,1); 
  elseif example == "rand_matlab"   
    # read saved matlab file
    vars = matread("data_Nx8_Nb4_Ncon10_Nth3_ub4_Nv6.mat") #read all data
    probm = vars["prob"]; #get(vars,"prob")
    f = vec(probm["f"]); 
    A = probm["A"]; b = vec(probm["b"]);         W = probm["W"];
    m = size(A,1);  n = size(f,1);          nth = size(W,2);
    H =zeros(n,n);  f_theta = zeros(n,nth); H_theta = zeros(0,0);
    vartype = vars["vartype"]; nb = length(vartype);
    bin_ind = collect(n-nb+1:n); #range(nb+1,n); #[nb+1:n];
  else
    #________________________________
    #Contrived QP example
    H = [0.97 0.19 0.15;0.19 0.98 0.05; 0.15 0.05 0.99];
    f = zeros(3); 
    f_theta = [11.3 -44.3; -3.66 -11.9; -32.6 7.81];
    H_theta = zeros(0,0);
    A = [2.5*0.15 2.5*0.88 2.5*0.17;0.49 0.57 0.22; 0.77 0.46 0.41];
    #b = zeros(3,1); ##for LP matrix!
    b = zeros(3); #,1 
    b[1:3] = [4.1;3.7;4.3];
    W = [0.19 -0.89; 0.64 -1.54; -0.59 -1.01];
    n = size(f,1); m = size(A,1); nth = size(W,2);
    bin_ind = [n-2;n-1;n];
    nb = length(bin_ind);
  end
  #________________________________________
  bounds_table = []; #collect(1:m);
  senses = zeros(Cint,m+2*nb); #for MI case #zeros(Cint,m); #for LP or QP
  #-----------------------------
  # Create mpQP
  mpMIQP = MPMIQP(H,f,f_theta,H_theta,A,b,W,bin_ind,bounds_table,senses)
  # Create region of interest and normalize primal and dual mpQP
  P_theta = (A = zeros(nth,0), b=zeros(0), ub=bnd*ones(nth),lb=-bnd*ones(nth)) 
  return mpMIQP,P_theta
end
# ______________________________________________________
#______________________________________________________
## Check containment in partition 
function param_location(th::Vector{Float64}, partition::Vector{BNBRegion};eps_gap=0.0)
  contained_in= Int64[]
  for (ind,region) in enumerate(partition)
    violation = minimum(region.bth-region.Ath'*th)
    if(violation>=-eps_gap)
      push!(contained_in,ind)
    end
  end
  return contained_in
end
#____________________________________________________
## Sampling from parameter space
function param_sampling(N_sam,P_theta, sample_way)
  nth = size(P_theta.lb,1);
  if (sample_way == "rand")
    th_sam = zeros(nth,N_sam);
    # random sampling in paramters space
    delta = P_theta.ub - P_theta.lb;
    for i = 1: N_sam
      th_sam[:,i] = P_theta.lb + rand(nth,nth)* (delta./2);
    end
    return th_sam
    #______________________________
    #_________________________________
  else
    # deterministic sampling 
    th_sam = zeros(nth,N_sam);   
    N = Int(floor(N_sam^(1/nth))); #otherwise N^nth will be generated 
    vect = range(P_theta.lb[1],P_theta.ub[1],length=N) #vect = range(P_theta.lb,P_theta.ub,length=N)
    if nth == 2 #2D parameter space
      th_all = collect(Iterators.product(vect, vect));
      for i = 1: length(th_all)
        th_sam[:,i] = [th_all[i][1];th_all[i][2]]
      end
      return th_sam
    elseif nth == 3 #3D
      th_all = collect(Iterators.product(vect, vect, vect));
      for i = 1: length(th_all)
        th_sam[:,i] = [th_all[i][1];th_all[i][2];th_all[i][end]]
      end
      return th_sam
    elseif nth == 4 #4D
      th_all = collect(Iterators.product(vect, vect, vect, vect));
      for i = 1: length(th_all)
        th_sam[:,i] = [th_all[i][1];th_all[i][2];th_all[i][3];th_all[i][end]]
      end
      return th_sam
    end 
  end
  
end

#___________________________________________
function fix_warm_AS(node,AS,lam,fix_id) 
  dlam = (node.A[AS,:]')\(node.A[fix_id,:]);
  block_ids = findall(dlam .< 0);
  if(isempty(block_ids)) 
    return Int64[]; # Dual unbounded => primal infeasible 
  end
  alphas = -lam[block_ids]./dlam[block_ids]; 
  alpha,rm_id = findmin(alphas);
  rm_id = block_ids[rm_id]; # Map id in block_ids to corresponding id in AS;
  inds_update = setdiff([1:length(AS);], rm_id) 
  return AS[inds_update] #,exitexec
end
#_________________________________________________
#_________________________________________________
## Find avergae of iteration number for final regions
function find_avg(parts)
  N = length(parts); iter_cert = zeros(N); tree_cert = zeros(N); iter_hu_tot_cert = zeros(N); 
  for j in eachindex(parts)
    iter_cert[j] = parts[j].iter
    tree_cert[j] = parts[j].kappa.node
    iter_hu_tot_cert[j] = parts[j].kappa.iter_hu_tot
    #iter_hu_cert[j] = parts[j].kappa.iter_hu
  end
  iter_avg = mean(iter_cert); #sum(iter_cert)/N
  tree_avg = mean(tree_cert); #sum(tree_cert)/N
  #iter_hu_avg = mean(iter_hu_cert); #sum(iter_hu_cert)/N
  iter_med = median(iter_cert); #sum(iter_cert)/N
  tree_med = median(tree_cert); #sum(tree_cert)/N
  iter_avg_hu_tot = mean(iter_hu_tot_cert); #sum(iter_cert)/N
  return iter_avg, tree_avg, iter_med, tree_med, iter_avg_hu_tot
end
#___________________________________________

#_________________________________________________
## Find the number of lock downs and up for diving heuristic
function count_locks(A)
  #using StatsBase #for geaometrical mean
  (nr,nc) = size(A)
  locks_down = zeros(1,nc); locks_up = zeros(1,nc); triv_round_down = false; triv_round_up = false;
  #___________________________
  #Find lock/down and update_sol
  for j in 1:nc
    locks_down[j] = count(x_0->x_0<=0,A[:,j]) #eps_zero
    locks_up[j] = count(x_0->x_0>=0,A[:,j])
  end
  (all(locks_down .== 0)) && (triv_round_down = true)
  (all(locks_up .== 0)) && (triv_round_up = true)
  return locks_down,locks_up, triv_round_down, triv_round_up
end
#___________________________________________
#_________________________________________________
## Find arithmatic and geometrical mean of Matrix A for each column
function mean_rel_matrix(A)
  #using StatsBase #for geaometrical mean
  (nr,nc) = size(A)
  #___________________________
  #Find relative value w-r-t no heuristic case for each experiment first, then get the arithmatic an geometric mean of it
  rel_diff_NH = zeros(nr,nc); perc= ones(nr,nc);  rel_wrt_NH = ones(nr,nc); 
  for j in 2:nc
    rel_diff_NH[:,j] = abs.(A[:,j].-A[:,1])./A[:,1]
    rel_wrt_NH[:,j] =  (A[:,j].+1)./(A[:,1].+1)
    perc[:,j] =  (A[:,j])./(A[:,1])
  end
  #___________________________
  arith_mean = zeros(1,nc); geom_mean = zeros(1,nc); arith_mean_rel = zeros(1,nc); geom_mean_rel = zeros(1,nc); arith_mean_perc = zeros(1,nc); geom_mean_perc = zeros(1,nc); arith_mean_rel_diff = zeros(1,nc); geom_mean_rel_diff = zeros(1,nc); 
  for j in 1:nc
    arith_mean[j] = mean(A[:,j])
    geom_mean[j] = geomean(A[:,j])
    arith_mean_rel_diff[j] = mean(rel_diff_NH[:,j])
    geom_mean_rel_diff[j] = geomean(rel_diff_NH[:,j])
    arith_mean_rel[j] = mean(rel_wrt_NH[:,j])
    geom_mean_rel[j] = geomean(rel_wrt_NH[:,j])  
    arith_mean_perc[j] = mean(perc[:,j])
    geom_mean_perc[j] = geomean(perc[:,j])  
  end
  return [round.(arith_mean, digits = 2); round.(geom_mean, digits = 2);round.(arith_mean_rel_diff, digits = 3); round.(geom_mean_rel_diff, digits = 3);round.(arith_mean_rel, digits = 3); round.(geom_mean_rel, digits = 3);round.(arith_mean_perc, digits = 3); round.(geom_mean_perc, digits = 3)]
end
#___________________________________________

