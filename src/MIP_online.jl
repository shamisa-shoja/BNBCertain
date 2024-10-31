# Use Gurobi to solve MIP problem online to compare with online bnb code
function MIP_online(prob, opts) 
    n = length(prob.f)
    m = length(prob.b)
    nb = size(prob.bin_ind,1); nc = n-nb;
    Ac = prob.A[:,1:nc]; Ab =  prob.A[:,nc+1:n];
    b = prob.b; 
    fc = prob.f[1:nc]; fb =  prob.f[nc+1:n];
    Hc = prob.H[1:nc,1:nc]; Hb =  prob.H[nc+1:n,nc+1:n]; Hcb = prob.H[1:nc,nc+1:n];

    if opts.chosen_solver == "Gurobi"
        model = Model(Gurobi.Optimizer)
    elseif opts.chosen_solver == "GLPK"
        model = Model(GLPK.Optimizer)  #Define the model for LP
    else
        model = Model(Ipopt.Optimizer)   #Define the model for nonlinear LP and QP
    end

    set_silent(model); #disable printing output from the solver.
    #unset_silent(model) #enable ..

    @variable(model, xc[1:nc]) #Create the variables #>=0
    @variable(model, xb[1:nb], Bin) #Create binary variables #>=0 #x = [xc; xb];
  
    for i in 1:m
        @constraint(model, (sum(Ac[i,j]*xc[j] for j in 1:nc)+ sum(Ab[i,j]*xb[j] for j in 1:nb)) <= b[i])
    end
    
    if opts.problem_type_QP 
        # QP
         @objective(model,Min, sum(0.5*Hc[i,j]*xc[i]*xc[j] for j in 1:nc, i in 1:nc) + sum(0.5*Hb[i,j]*xb[i]*xb[j] for j in 1:nb, i in 1:nb) + sum(Hcb[i,j]*xc[i]*xb[j] for j in 1:nb, i in 1:nc) + sum(fc[i]*xc[i] for i in 1:nc)+ sum(fb[i]*xb[i] for i in 1:nb)) # objective FUNCTION
    else #LP
       @objective(model,Min, sum(fc[i]*xc[i] for i in 1:nc)+ sum(fb[i]*xb[i] for i in 1:nb)) # objective FUNCTION
    end
        #@show model  #summery
    if opts.chosen_solver == "Gurobi"
        set_optimizer_attribute(model, "NonConvex", 2) # to solve nonconvex QP comparision by Gurobi
    end
    
    optimize!(model) # Solve 
    exitflag = termination_status(model)  #status
    #display(typeof(MOI.OPTIMAL))  # status code

    x_sol = zeros(n); obj_val = 0; AS = falses(length(prob.b));
    if string(exitflag) == "OPTIMAL"
        xc_sol = value.(xc);
        xb_sol = value.(xb);
        x_sol = [xc_sol;xb_sol];
        #@show objective_value(model)
        obj_val = objective_value(model)  #objective function
        zero_tol = 1e-6;
        cons = abs.(prob.A*x_sol - prob.b)
        AS = cons.<=zero_tol;  #active set
    end

    return exitflag, x_sol, obj_val, AS
end