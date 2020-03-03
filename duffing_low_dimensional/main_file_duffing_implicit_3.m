%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% - File that performs control node selection for the Duffing network whose
% continuous-time dynamics is:
% \dot{x}=f(x)+Bu
% 
% where x is the state, f(x) is the vector fuction describing an
% uncontrolled dynamics, B is the control matrix and u is the control input
% vector.
%
% - The discretization of the dynamics is performed using the Trapezoidal
% Implicit (TI) method
%
% - Before running this file, run  "generate_dynamics_duffing.m" file to
% generate two files:
%
%    1.) "duffing_network_dynamics.m" - function that computes f(x) for a given x
%    2.) "duffing_network_dynamics_gradient" - function that computes
%    \nabla f(x) (gradient) for a given x
%
%  IN THIS FILE WE GENERATE A SMALLER NETWORK AND PERFORM AN EXHAUSTIVE
%  SEARCH (ALSO KNOWN AS A BRUTE-FORCE SEARCH)
% - Author: Aleksandar Haber
% December 2019 - February 2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear, pack, clc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                   parameter selection
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% N is the number of subsystems, every subsystem is of the second order
% the precise value of N should be adjusted such that it matches the number
% of subsystems in "duffing_network_dynamics.m" and "duffing_network_dynamics_gradient"
N=10
% discretization constant 
h=0.01
% number of time steps for control and control node selection
time_steps=10; 
% number of control nodes
no_control_nodes=2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                   end of parameter selection
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rng('shuffle') 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   computation of initial and desired states and additional tests
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% simulate the uncontrolled dynamics to obtain the system equilibrium point
initial_point=0.5*rand(2*N,1);
time=0:h:300*h;
[time_tmp,eq_point] = simulate_uncontrolled_ode45(time,initial_point,@duffing_network_dynamics);

hold on
plot(eq_point(:,:),'m')

% compare the TI method with ode45 simulation
STATE_fE=simulate_uncontrolled_ti_fsolve_3(length(time),initial_point,h,@duffing_network_dynamics,@duffing_network_dynamics_gradient);
STATE_fE=STATE_fE';

for i=1:length(time)
   error_simulation(i)= norm(STATE_fE(i,:)-eq_point(i,:),2);%/norm(eq_point(i,:),2);
end

figure(1)
hold on
plot(error_simulation,'k');
figure(2)
plot(eq_point(:,2),'k')
hold on 
plot(STATE_fE(:,2),'m')

% initial state
x0=eq_point(end,:)';
% desired state
xd=initial_point;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% end of computation of initial and desired states and additional tests
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% computation of the initial control sequence when all the nodes are
% controlled
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% initial control solution
initial_solution=randn(N*time_steps,1);
% form the control input matrix when all the nodes are controlled
Bmatrix=formBmatrix_3(N*2,ones(1,N));
% this is the control solution - this vector contains the control inputs 
solution_initial(:,1)=solve_problem_dynamic_unconstrained_ti_3(time_steps,x0,xd,h,Bmatrix,initial_solution,@duffing_network_dynamics,@duffing_network_dynamics_gradient)
% compute the response for the calculated control inputs
STATE_controlled=simulate_controlled_ti_fsolve_3(time_steps,x0,h,Bmatrix,solution_initial(:,1),@duffing_network_dynamics,@duffing_network_dynamics_gradient);
% compute the control error
for i=1:time_steps
    error(i)=norm(STATE_controlled(:,i)-xd(:,1));
end
plot(error)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% end of the computation of the initial control sequence when all the nodes are
% controlled
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%        iterative solution of the control node selection problem
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:3
% generate the initial solution for optimization

% main variables

% - "solution_initial_selection"
% every column of "solution_initial_selection" contains two parts-
% the first part (N \times 1) is the vector of binary variables describing the nodes
% that are controlled (1 - controlled, 0 - not controlled)
% the second part (time_steps \times 1) is the vector of control sequences
% this vector is used as an initial solution for the mixed integer
% optimization problem

% - "solution_selection" is the solution vector of the mixed integer
% problem - its structure corresponds to the structure of the vector 
% "solution_initial_selection"

if i==1
    % no need to map back anything
    solution_initial_selection(:,i)=[ones(N,1); solution_initial(:,i)];
else
    % map back the reduced control input sequences to the larger vector and add
    % zeros for the nodes that are not controlled
    for idx1=1:time_steps
        solution_initial((idx1-1)*N+nodes_that_are_controlled,i)= solution_initial_reduced_tmp((idx1-1)*no_control_nodes+1:idx1*no_control_nodes,1);
        solution_initial((idx1-1)*N+nodes_that_are_not_controlled,i)= zeros(numel(nodes_that_are_not_controlled),1);
    end
    solution_initial_selection(:,i)=[nodes_controlled_and_not ; solution_initial(:,i)];
end
% here we solve the mixed integer optimization problem
[solution_selection(:,i),fval(i),exitflag(i),info{i}]=solve_problem_dynamic_ti_fsolve_3(time_steps,x0,xd,h,solution_initial_selection(:,i),no_control_nodes,@duffing_network_dynamics,@duffing_network_dynamics_gradient);


nodes_controlled_and_not=solution_selection(1:N,i);
control_input=solution_selection(N+1:end,i);

% reduce the Bmatrix and the control input to exclude the nodes that are
% not controlled - this is used to compute the control sequence for the
% next iteration of the optimization problem
nodes_that_are_controlled = find(nodes_controlled_and_not>0);
nodes_that_are_not_controlled=find(nodes_controlled_and_not==0);


Bmatrix_reduced=formBmatrix_3_reduced(2*N,nodes_controlled_and_not);

% this vector is the vector obtained from the solution of the mixed integer
% optimization problem by eliminating the control sequences corresponding
% to the nodes that are not controlled
control_input_reduced=[];
for i1=1:time_steps
    for i2=1:numel(nodes_that_are_controlled)
        control_input_reduced=[control_input_reduced; control_input(((i1-1)*N)+nodes_that_are_controlled(i2))];
    end
end

% here we simulate the dynamics for the Bmatrix_reduced and
% control_input_reduced to obtain the controlled state trajectory and to
% compute the final control error
clear STATE_controlled;
STATE_controlled=simulate_controlled_ti_fsolve_3(time_steps,x0,h,Bmatrix_reduced,control_input_reduced,@duffing_network_dynamics,@duffing_network_dynamics_gradient);
error_iterations(2*i-1)=norm(STATE_controlled(:,end)-xd(:,1));

% for computed Bmatrix_reduced, we compute a new control sequence that is
% used as an initial guess for the next mixed integer optimization problem
clear STATE_controlled;
solution_initial_reduced_tmp=solve_problem_dynamic_unconstrained_ti_3(time_steps,x0,xd,h,Bmatrix_reduced,control_input_reduced,@duffing_network_dynamics,@duffing_network_dynamics_gradient);
STATE_controlled=simulate_controlled_ti_fsolve_3(time_steps,x0,h,Bmatrix_reduced,solution_initial_reduced_tmp,@duffing_network_dynamics,@duffing_network_dynamics_gradient);
error_iterations(2*i)=norm(STATE_controlled(:,end)-xd(:,1));

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                   final computations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% compute the control error
for i=1:time_steps
    error(i)=norm(STATE_controlled(:,i)-xd(:,1));
end
plot(error)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%        verify the results by exhaustive search 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

v = 1:1:N;
% generate all the possible selections of "no_control_nodes" control nodes
% out of possible "N" control nodes
nodes_random_selection = nchoosek(v,no_control_nodes)


for i=1:numel(nodes_random_selection(:,1))
i
Brandom=zeros(2*N,no_control_nodes);

for j=1:no_control_nodes
    Brandom(2*nodes_random_selection(i,j),j)=1;
end
solution_random=solve_problem_dynamic_unconstrained_ti_3(time_steps,x0,xd,h,Brandom,randn(time_steps*no_control_nodes,1),@duffing_network_dynamics,@duffing_network_dynamics_gradient);
% compute the response for the calculated control inputs

STATE_controlled_random=simulate_controlled_ti_fsolve_3(time_steps,x0,h,Brandom,solution_random,@duffing_network_dynamics,@duffing_network_dynamics_gradient)
% compute the control error
for indx1=1:time_steps
    error_random(i,indx1)=norm(STATE_controlled_random(:,indx1)-xd(:,1));
end
end

hist(error_random(:,time_steps))
histogram(error_random(:,time_steps), 10, 'Normalization','probability' )
ytix = get(gca, 'YTick')
set(gca, 'YTick',ytix, 'YTickLabel',ytix*100)
hold on 
line([error(end) error(end)], [0 0.25],'color','r','linewidth',2);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  end verify the results by randomly selecting the control nodes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





