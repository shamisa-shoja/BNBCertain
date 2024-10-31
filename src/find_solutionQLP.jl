# Use to solve QP relaxation or do the QP comparision

function QLPSolve(prob, opts) 
    n = length(prob.f)
    m = length(prob.b)
    if opts.chosen_solver == "Gurobi"
        model = Model(Gurobi.Optimizer)
    elseif opts.chosen_solver == "GLPK"
        model = Model(GLPK.Optimizer)  #Define the model for LP
    else
        model = Model(Ipopt.Optimizer)   #Define the model for nonlinear LP and QP
    end

    set_silent(model); #disable printing output from the solver.
    #unset_silent(model) #enable ..
    
    @variable(model, x[1:n]) #Create the variables #>=0
    @constraint(model, prob.A*x  .<= prob.b) 
    if opts.problem_type_QP 
        # QP
        @objective(model,Min, sum(prob.H[i,j]*x[i]*x[j] for j in 1:n, i in 1:n) + sum(prob.f[i]*x[i] for i in 1:n)) # objective FUNCTION
    else #LP
        @objective(model,Min, sum(prob.f[i]*x[i] for i in 1:n)) # objective FUNCTION
    end
        #@show model  #summery
    if opts.chosen_solver == "Gurobi"
        set_optimizer_attribute(model, "NonConvex", 2) # to solve nonconvex QP comparision by Gurobi
    end
    
    optimize!(model) # Solve 
    exitflag = termination_status(model)  #status
    #display(typeof(MOI.OPTIMAL))  # status code

    x_var = zeros(n); obj_var = 0; AS = falses(length(prob.b));
    if string(exitflag) == "OPTIMAL"
        for i in 1:n 
            x_var[i] = value(x[i]) #solution result
        end
        #@show objective_value(model)
        obj_var = objective_value(model)  #objective function
        zero_tol = 1e-6;
        cons = abs.(prob.A*x_var - prob.b)
        AS = cons.<= 0; #zero_tol;  #active set
    end

    return exitflag, x_var, obj_var, AS
end