% function that computes the control solution for a given control matrix B
% - input parameters: 
%                   - time_steps         - discrete-time simulation time 
%                   - x0                 - initial state
%                   - xd                 - desired state
%                   - h                  - discretization constant
%                   - Bmatrix            - control matrix
%                   - initial_solution   - initial solution of the control
%                   problem
%                   - fcnHandle - function handle that describes the system
%                   dynamics
%                   - fcnHandleGradient - function handle that described
%                   the gradient of the system dynamics
% - output parameters: 
%                   - solution           - solution of the control problem
% Author: Aleksandar Haber 
% December 2019 - February 2020

function solution=solve_problem_dynamic_unconstrained_ti_3(time_steps,x0,xd,h,Bmatrix,initial_solution,fcnHandle,fcnHandleGradient)

options = optimoptions(@fminunc,'Algorithm','quasi-newton','Display','iter','MaxFunctionEvaluations',30000,'UseParallel',true)
[n,~]=size(x0);
I=speye(n,n);
[~,no_controlled_nodes]=size(Bmatrix);

% this is a vector of lifted desired states that is necessary for defining
% the cost function
Xd=kron(ones(time_steps-1,1),xd);

[solution,fval] = fminunc(@(z)cost_function_dynamic_unconstrained(z),initial_solution,options);

% this function computes the objective function for the control input
% sequence z
function cost_function_value_dynamic = cost_function_dynamic_unconstrained(z)
    
% this variable contains the computed final states    
STATE=zeros(n,time_steps-1);

% parameters for the fsolve function - the objective function is defined at
% the end of the file
    options_fsolve = optimoptions('fsolve','Algorithm', 'trust-region','Display','off','SpecifyObjectiveGradient',true,'UseParallel',true,'FunctionTolerance',1.0000e-8,'MaxIter',10000,'StepTolerance', 1.0000e-8);
    

    problem.options = options_fsolve;
    problem.objective = @objective_fun;
    problem.solver = 'fsolve';
 
    for o=1:(time_steps-1)
        if o==1
           X=x0;
           tmp0=fcnHandle(X); 
           problem.x0 = X+h*(tmp0+Bmatrix*z((o-1)*no_controlled_nodes+1:o*no_controlled_nodes,1)) ;  % use the Forward Euler method to generate the initial guess
           X=fsolve(problem);
           STATE(:,o)=X;
        else
           tmp0=fcnHandle(X); 
           problem.x0 = X+h*(tmp0+Bmatrix*z((o-1)*no_controlled_nodes+1:o*no_controlled_nodes,1)) ;  % use the Forward Euler method to generate the initial guess
           X=fsolve(problem);
           STATE(:,o)=X; 
        end
    end
    cost_function_value_dynamic=norm(Xd-STATE(:),2)^2;
function [f, jacobianz]=objective_fun(xk)
    tmp1=0.5*h*(fcnHandle(xk)+tmp0);% tmp0=fcnHandle(xk-1) is computed outside this function to speed up the computations (in the function above)
    f=xk-X-tmp1- 0.5*h*Bmatrix*(z((o-1)*no_controlled_nodes+1:(o)*no_controlled_nodes,1)+z((o)*no_controlled_nodes+1:(o+1)*no_controlled_nodes,1));
    jacobianz=I-0.5*h*fcnHandleGradient(xk)'; % here we need to transpose since the function "fcnHandleGradient(z)" returns the gradient, and the Jacobian is its transpose
end
end
end