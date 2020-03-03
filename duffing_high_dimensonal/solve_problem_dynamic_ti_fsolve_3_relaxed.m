% - function that solves the control node selection and control problems
% using the mixed integer optimization approach
% - input parameters: 
%                   - time_steps        - discrete-time simulation time 
%                   - x0                - initial state
%                   - xd                - desired state
%                   - h                 - discretization constant
%                   - initial_solution  - initial solution for the control
%                                         problem
%                   - no_control_nodes     - number of control nodes
%                   - fcnHandle         - function handle that describes the system
%                                         dynamics
%                   - fcnHandleGradient - function handle that described
%                                         the gradient of the fcnHandle
% - output parameters:
%                   - solution          - solution vector- first N 
                        %                 variables are binary and they
                        %                 represent the control selections
                        %                 0 - node is not controlled and 1
                        %                 node is controlled
                        %                 "N+1" until "end" variables are control sequences (N*time_steps) 
%                   - fval              - final function value
%                   - exitflag          - solver's flag
%                   - info              - additional info
% Author: Aleksandar Haber
% December 2019 - February 2020
function [solution,fval,exitflag,info]=solve_problem_dynamic_ti_fsolve_3_relaxed(time_steps,x0,xd,h,initial_solution,no_control_nodes,fcnHandle,fcnHandleGradient)

[n,~]=size(x0);
I=speye(n,n);
N=n/2;

nB = N; %Number of Binary Variables
nC = N*time_steps; %Number of Continuous Variables

tol=1e-2; % this is a tolerance for relaxing the equality constraints
% not the equality constraints are posed as double inequality constraints
% this is a standard practice for relaxing the equality constraints
% without the tolerance the solver returns infeasible solutions...
xtype = [repmat('B',1,nB),repmat('C',1,nC)];
Aeq=[ones(1,nB), zeros(1,nC) ]; beq=[no_control_nodes+tol];

% Options for the mixed integer optimization problem
opts = optiset('solver','nomad','display','iter','maxfeval',200000,'maxtime',36000,'maxnodes',20000)
Opt = opti('fun',@cost_function_dynamic,'ineq', Aeq,beq,'xtype',xtype,'options',opts);

% this is a vector of lifted desired states that is necessary for defining
% the cost function
Xd=kron(ones(time_steps-1,1),xd);


%Solve the mixed integer optimization problem
[solution,fval,exitflag,info] = solve(Opt,initial_solution);

% this is the cost function
function cost_function_value_dynamic = cost_function_dynamic(z)
    
    % this variable contains the computed final states    
      STATE=zeros(n,time_steps-1);
    % optimization options for fsolve
    options_fsolve = optimoptions('fsolve','Algorithm', 'trust-region','Display','off','SpecifyObjectiveGradient',true,'UseParallel',true,'FunctionTolerance',1.0000e-8,'MaxIter',10000,'StepTolerance', 1.0000e-8);
    
    problem.options = options_fsolve;
    problem.objective = @objective_fun;
    problem.solver = 'fsolve';
             
    Bz=formBmatrix_3(n,z(1:N));

    for o=1:(time_steps-1)
        if o==1
           X=x0;
           tmp0=fcnHandle(X); 
           problem.x0 = X+h*(tmp0+Bz*z(o*N+1:(o+1)*N)) ;  % use the Forward Euler method to generate the initial guess
           X=fsolve(problem);
           STATE(:,o)=X;
        else
            tmp0=fcnHandle(X); 
            problem.x0 = X+h*(tmp0+Bz*z(o*N+1:(o+1)*N)) ;  % use the Forward Euler method to generate the initial guess
            X=fsolve(problem);
            STATE(:,o)=X;
        end
    end
    
    cost_function_value_dynamic=norm(Xd-STATE(:),2)^2;
    
function [f, jacobianz]=objective_fun(xk)
    tmp1=0.5*h*(fcnHandle(xk)+tmp0);% tmp0=fcnHandle(xk-1) is computed outside this function to speed up the computations (in the function above)
    f=xk-X-tmp1- 0.5*h*Bz*(z(o*N+1:(o+1)*N,1)+z((o+1)*N+1:(o+2)*N,1));
    jacobianz=I-0.5*h*fcnHandleGradient(xk)'; % here we need to transpose since the function "fcnHandleGradient(z)" returns the gradient, and the Jacobian is its transpose
end
end   
end