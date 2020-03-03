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
% - output parameters: 
%                   - solution           - solution of the control problem
% Author: Aleksandar Haber 
% December 2019 - February 2020

function solution=solve_problem_dynamic_unconstrained(time_steps,x0,xd,h,Bmatrix,initial_solution,fcnHandle)

options = optimoptions(@fminunc,'Algorithm','quasi-newton','Display','iter','MaxFunctionEvaluations',30000,'UseParallel',true)

% this is a vector of lifted desired states that is necessary for defining
% the cost function
Xd=kron(ones(time_steps,1),xd);

[N,~]=size(x0);
[~,no_controlled_nodes]=size(Bmatrix);
[solution,fval] = fminunc(@(z)cost_function_dynamic_unconstrained(z),initial_solution,options);

function cost_function_value_dynamic = cost_function_dynamic_unconstrained(z)

% this variable contains the computed final states    
STATE=zeros(N,time_steps);
    for o=1:time_steps
        if o==1
           X=x0+h*(fcnHandle(x0)+Bmatrix*z((o-1)*no_controlled_nodes+1:o*no_controlled_nodes,1));
           STATE(:,o)=X;
        else 
           X=X+h*(fcnHandle(X)+Bmatrix*z((o-1)*no_controlled_nodes+1:(o)*no_controlled_nodes,1));
           STATE(:,o)=X;
        end
    end
    cost_function_value_dynamic=norm(Xd-STATE(:),2)^2;
end
end