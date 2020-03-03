% - function that solves the control node selection and control problems
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
function [solution,fval,exitflag,info]=solve_problem_dynamic(time_steps,x0,xd,h,initial_solution,no_control_nodes,fcnHandle)

[N,~]=size(x0);

nB = N; %Number of Binary Variables
nC = N*time_steps; %Number of Continuous Variables

tol=1e-2; % this is a tolerance for relaxing the equality constraints
% not the equality constraints are posed as double inequality constraints
% this is a standard practice for relaxing the equality constraints
% without the tolerance the solver returns infeasible solutions...
xtype = [repmat('B',1,nB),repmat('C',1,nC)];
Aeq=[ones(1,nB), zeros(1,nC) ; -ones(1,nB), zeros(1,nC) ]; beq=[no_control_nodes+tol; - no_control_nodes+tol];

% Options for the mixed integer optimization problem
opts = optiset('solver','nomad','display','iter','maxfeval',100000)
Opt = opti('fun',@cost_function_dynamic,'ineq', Aeq,beq,'xtype',xtype,'options',opts);

% this is a vector of lifted desired states that is necessary for defining
% the cost function
Xd=kron(ones(time_steps,1),xd);


%Solve the mixed integer optimization problem
[solution,fval,exitflag,info] = solve(Opt,initial_solution);


function cost_function_value_dynamic = cost_function_dynamic(z)
   % this variable contains the computed final states    
   STATE=zeros(N,time_steps);
    
    Bz=diag(z(1:N));
    for o=1:time_steps
        if o==1
           X=x0+h*(fcnHandle(x0)+Bz*z(o*N+1:(o+1)*N));
           STATE(:,o)=X;
        else 
            X=X+h*(fcnHandle(X)+Bz*z(o*N+1:(o+1)*N));
            STATE(:,o)=X;
        end
    end
   cost_function_value_dynamic=norm(Xd-STATE(:),2)^2;
end
end