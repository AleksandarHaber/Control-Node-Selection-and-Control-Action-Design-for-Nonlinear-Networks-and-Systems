% - function that simulates the controlled dynamics using the trapezoidal
% implicit (TI) method
% - input parameters: 
%                   - time_steps    - discrete-time simulation time 
%                   - x0            - initial state
%                   - h             - discretization constant
%                   - Bmatrix       - control input matrix
%                   - control       - control input vector
%                   - fcnHandle     - function handle that describes the
%                   system dynamics
%                   - fcnHandleGradient - function handle that describes
%                   the gradient of the system dynamics
% - output parameters:
%                   - STATE         - state trajectory
% - Author: Aleksandar Haber
% December 2019 - February 2020

function STATE=simulate_controlled_ti_fsolve_3(time_steps,x0,h,Bmatrix,control,fcnHandle,fcnHandleGradient)

[n,~]=size(x0);
I=speye(n,n);
STATE=zeros(n,time_steps);
[~,no_controlled_nodes]=size(Bmatrix);

% here we adjust the options of the MATLAB method for solving nonlinear
% system of equations
 options_fsolve = optimoptions('fsolve','Algorithm', 'trust-region','Display','off','SpecifyObjectiveGradient',true,'UseParallel',true,'FunctionTolerance',1.0000e-8,'MaxIter',10000,'StepTolerance', 1.0000e-8);
      
problem.options = options_fsolve;
problem.objective = @objective_fun;
problem.solver = 'fsolve';


for o=1:(time_steps-1)
        if o==1
           STATE(:,o)=x0; 
           tmp0=fcnHandle(STATE(:,o));
           problem.x0 = STATE(:,o)+h*(tmp0+Bmatrix*control((o-1)*no_controlled_nodes+1:o*no_controlled_nodes,1)) ;  % use the Forward Euler method to generate the initial guess
           STATE(:,o+1)=fsolve(problem);
        else
           tmp0=fcnHandle(STATE(:,o));
           problem.x0 = STATE(:,o)+h*(tmp0+Bmatrix*control((o-1)*no_controlled_nodes+1:o*no_controlled_nodes,1)) ;  % use the Forward Euler method to generate the initial guess
           STATE(:,o+1)=fsolve(problem);
        end
end

function [f, jacobianz]=objective_fun(xk)
    tmp1=0.5*h*(fcnHandle(xk)+tmp0);% tmp0=fcnHandle(xk-1) is computed outside this function to speed up the computations (in the function above)
    f=xk-STATE(:,o)-tmp1- 0.5*h*Bmatrix*(control((o-1)*no_controlled_nodes+1:(o)*no_controlled_nodes,1)+control((o)*no_controlled_nodes+1:(o+1)*no_controlled_nodes,1));
    jacobianz=I-0.5*h*fcnHandleGradient(xk)'; % here we need to transpose since the function "fcnHandleGradient(z)" returns the gradient, and the Jacobian is its transpose
end

end