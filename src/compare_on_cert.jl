function compare_on_cert(final_parts, prob, ths, strategies, opts, points_to_compare)      
    #final_parts = deepcopy([parts]); #for debug
    (isempty(final_parts[1])) && (final_parts = final_parts[2:end])
    N = size(ths,2); n_test = length(final_parts); 
    max_diff =zeros(n_test); min_diff = zeros(n_test); mean_diff = zeros(n_test); Nr_diff_iter= zeros(n_test);   
    info_diff = zeros(N,n_test*6);  #opts.online_gurobi ? info_diff = zeros(N,n_test*6) : info_diff = zeros(N,n_test*5); #save data of comparision with gurobi
    A_eps = opts.tol_domcut_abs_param[:,1:end-1]; B_eps = opts.tol_domcut_abs_param[:,end]; 
    opts.chosen_solver = "Gurobi"; th_chebyC =[]; max_iter_hu_on = zeros(n_test); max_iter_tot_on = zeros(n_test);  avg_time = 0;
    ## ___________________________________________________________________________________
    if points_to_compare[1]
        # compare online and offline iteration resultsN = size(ths,2);         
        for i in 1:n_test
            opts.sorting_method, opts.warm_start = strategies[i,:]; 
            ind_sam = zeros(Int64,N); iter_on_sam = zeros(Float64,N); iter_off_sam= zeros(Float64,N); gen_nodes = zeros(Int64,N); iter_hu_on = zeros(Int64,N);
            J_on_diff = zeros(Float64,N); x_on_diff = zeros(Float64,N); iter_hu_cert= zeros(Float64,N);  iter_tot_on = zeros(Int64,N); iter_hu_tot_on = zeros(Int64,N); t_onn=zeros(Float64,N); 
            part = deepcopy(final_parts[i]); #deepcopy(final_parts[i+1]); #first one is empty #pop!(final_part); #pop from end of array
            for j = 1: N #j=1; #j = j+1;
                #println(j);  #for debug
                th_sam = ths[:,j];
                #__________ offline: ponit location ____________       
                ind_sam1 = param_location(th_sam,part,eps_gap=opts.eps_gap) # 2* # find point in final partition
                if (length(ind_sam1)>0)
                    ind_sam[j] = ind_sam1[1];
                    iter_off_sam[j] = part[ind_sam1[1]].iter; 
                    iter_hu_cert[j] = part[ind_sam1[1]].kappa.iter_hu; 
                    #___________ online ____________        
                     #TODO: uncomment below for suboptimal case
                        #(opts.tol_abs_param) && (opts.tol_domcut_abs = th_sam'*A_eps*th_sam + B_eps'*th_sam + opts.tol_domcut_abs); #parametric absolute allowance
                        t3 = time(); 
                        x_on,J_on, iter_sam, iter_hu, iter_tot, tot_nodes, seq_nodes, tree_status, _, decomp_nodes, tot_gen_node, iter_hu_tot = bnb_on_hu(MIQP(prob, th_sam), opts);  #bnb_online_subopt_hu
                        t_onn[j] = time() - t3;  
                    iter_on_sam[j] = iter_sam; gen_nodes[j] = tot_gen_node;  iter_hu_on[j] = iter_hu;  iter_tot_on[j] = iter_tot; iter_hu_tot_on[j] = iter_hu_tot
                    #________________________________Online Gurobi_______________________________________                  
                    (opts.online_gurobi) && ((exitflag_g, xopt_g, Jopt_g, AS_g) = MIP_online(MIQP(prob, th_sam), opts)) #@time ()
                    (opts.online_gurobi) && (J_on_diff[j]=J_on-Jopt_g; x_on_diff[j]= norm(x_on-xopt_g,Inf))                    
                end #if
            end #for N
            #______________________
            # compare results
            diff_iters = iter_off_sam - iter_on_sam;  
            diff_ind = findall(diff_iters .!= 0); 
            if !isempty(diff_ind) #keep track of differences if there is any
                Nr_diff_iter[i] = length(diff_ind); 
                max_diff[i] = maximum(diff_iters); min_diff[i] = minimum(diff_iters); mean_diff[i] = (sum(abs.(diff_iters))/N)*100
            end
            diff_iters_hu = iter_hu_cert - iter_hu_on; 
            max_iter_hu_on[i] = maximum(iter_hu_on);
            max_iter_tot_on[i] = maximum(iter_tot_on); 
            #opts.arith_mean ? avg_time = mean(t_onn) : avg_time = geomean(t_onn)
            (opts.storage_level>=1 && !opts.online_gurobi) && (info_diff[:,6*i-5]= diff_iters;  info_diff[:,6*i-4]= ind_sam; info_diff[:,6*i-3]= iter_on_sam; info_diff[:,6*i-2]= t_onn; info_diff[:,6*i-1]= iter_hu_on; info_diff[:,6*i]= diff_iters_hu); #;); #gen_nodes); #(opts.storage_level>=1) && (pushfirst!(info_diff, [diff_iters ind_sam]));
            (opts.storage_level>=1 && opts.online_gurobi) &&  (info_diff[:,6*i-5]= diff_iters;  info_diff[:,6*i-4]= ind_sam; info_diff[:,6*i-3]= iter_on_sam; info_diff[:,6*i-2]= iter_hu_on; info_diff[:,6*i-1]= x_on_diff; info_diff[:,6*i]= J_on_diff;);
        end# for i test
    end
    ## ___________________________________________________________________________________
    ##____________________________________________________________________________________
    if points_to_compare[2]  
        #compute chebychev centers and apply online bnb
        max_diff_chC =zeros(n_test); min_diff_chC = zeros(n_test); mean_diff_chC = zeros(n_test); Nr_diff_iter_chC= zeros(n_test);   
        info_diff_chC = []; #
        N_reg = zeros(n_test); 
        for i in 1:n_test
            opts.sorting_method, opts.warm_start = strategies[i,:]; 
            part = deepcopy(final_parts[i]); 
            n_reg = length(part); N_reg[i] = n_reg; 
            iter_on_sam_chC = zeros(Float64,n_reg); iter_off_sam_chC= zeros(Float64,n_reg); gen_nodes_chC = zeros(Int64,N);  iter_hu_on = zeros(Int64,N);
            #___________compute chebychev centers_________
            th_chebyC, ~ =  compute_centers(part); 
            #______________________
            for j = 1: n_reg
                th_sam = th_chebyC[j];
                (isempty(th_sam)) && continue;
                #__________ offline: ponit location ____________       
                iter_off_sam_chC[j] = part[j].iter; 
                #___________ online ____________    
                _,_, iter_Ch, iter_huCh, tot_nodesCh, _, _, _, _, tot_gen_nodeCh = bnb_on_hu(MIQP(prob, th_sam), opts);  #bnb_online_subopt_hu
                iter_on_sam_chC[j] = iter_Ch; gen_nodes_chC[j] = tot_gen_nodeCh; iter_hu_on[j] = iter_huCh;                
            end #for N
            #______________________
            # compare results
            diff_iters_chC = iter_off_sam_chC - iter_on_sam_chC;  
            diff_ind_chC = findall(diff_iters_chC .!= 0); 
            if !isempty(diff_ind_chC) #keep track of differences if there is any
                Nr_diff_iter_chC[i] = length(diff_ind_chC); 
                max_diff_chC[i] = maximum(diff_iters_chC); min_diff_chC[i] = minimum(diff_iters_chC); mean_diff_chC[i] = (sum(abs.(diff_iters_chC))/N)*100
            end
            (opts.storage_level>=1) && (push!(info_diff_chC, diff_iters_chC)); #, gen_nodes_chC));
         end# for test
            ##   ______________ collect sampling and chebychev results ____________________
        compare_info = Array{Any}(undef, 2*n_test+1, 7)
        compare_info[1, :] = [" #diff_iter "  " max diff "  " min diff " " mean dif %" " #points " "----bnb strategy----" "CS:false, WS:true"]
        compare_info[2:end, :] =  [Nr_diff_iter max_diff min_diff mean_diff N*ones(n_test) strategies; Nr_diff_iter_chC max_diff_chC min_diff_chC mean_diff_chC N_reg strategies]
        return compare_info, info_diff, max_iter_hu_on, max_iter_tot_on, avg_time, th_chebyC, info_diff_chC # diff_iters_chC ,, th_chebyC_all
    else
            ##   ______________ collect sampling results  ____________________
        compare_info = Array{Any}(undef, n_test+1, 7)
        compare_info[1, :] = [" #diff_iter "  " max diff "  " min diff " " mean dif %" " #points " "----bnb strategy----" "CS:false, WS:true"]
        compare_info[2:end, :] =  [Nr_diff_iter max_diff min_diff mean_diff N*ones(n_test) strategies]
        return compare_info, info_diff, max_iter_hu_on, max_iter_tot_on, avg_time, [] , []
    end
end