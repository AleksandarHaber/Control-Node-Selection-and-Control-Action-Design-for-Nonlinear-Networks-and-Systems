%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% file that performs control node selection for the memory network whose
% continuous-time dynamics is defined in the file:
% "memory_network_dynamics.m"
% Author: Aleksandar Haber
% December 2019 - February 2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear, pack, clc
% The size of the network is N=N1xN1 
N1=5; N=N1*N1;

% discretization constant for the forward Euler method
h=0.01

% number of time steps for control and control node selection
time_steps=10; 
rng('shuffle') 
% number of control nodes
no_control_nodes=10;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                   letter definitions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define the letters
L=[1 -1 -1 -1 -1;
   1 -1 -1 -1 -1; 
   1 -1 -1 -1 -1;
   1 -1 -1 -1 -1;
   1  1  1  1  1;];

T=[ 1  1  1   1  1;  
   -1 -1  1  -1  -1;
   -1 -1  1  -1  -1;
   -1 -1  1  -1  -1;
   -1 -1  1  -1  -1;];
    
H=[ 1 -1 -1 -1  1 ;
    1 -1 -1 -1  1 ;
    1  1  1  1  1 ;
    1 -1 -1 -1  1 ;
    1 -1 -1 -1  1 ;
    ]
E=[1   1   1   1   1;
   1  -1  -1   -1 -1;
   1   1   1   -1 -1;
   1  -1  -1   -1 -1;
   1   1   1    1  1;];

F=[1   1   1    1   1;
   1  -1  -1   -1  -1;
   1   1   1   -1  -1;
   1  -1  -1   -1  -1;
   1  -1  -1   -1  -1;];

S=[1  1  1  1  1;
   1 -1 -1 -1 -1;
   1  1  1  1  1;
  -1  -1 -1 -1 1;
  1    1  1  1 1;];

O=[1   1   1   1  1;
   1  -1  -1  -1  1;
   1  -1  -1  -1  1;
   1  -1  -1  -1  1;
   1   1   1   1  1;];

Y=[1 -1 -1 -1  1;
   -1 1 -1  1 -1;
   -1 -1 1 -1 -1;
   -1 -1 1 -1 -1;
   -1 -1 1 -1 -1;];

X=[1 -1 -1 -1 1;
   -1 1 -1  1 -1;
   -1 -1 1  -1 -1;
   -1 1 -1  1 -1;
   1 -1 -1 -1  1;];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                end of the letter definitions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  computation of initial and desired states and additional tests
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% the initial state is the equilibrium corresponding to H
% the desired state is the equilibrium corresponding to T

% simulate the dynamics to obtain the system equilibrium points
x_H=H(:)+0.5*randn(N,1);
time=0:h:1000*h;
[time_tmp,states_H] = simulate_uncontrolled_ode45(time,x_H,@memory_network_dynamics);

x_T=T(:)+0.5*randn(N,1);
time=0:h:1000*h;
[time_tmp,states_T]=simulate_uncontrolled_ode45(time,x_T,@memory_network_dynamics);

% - plot the uncontrolled network to verify the equilibrium points 
% - this piece of code is to verify that the equilibrium points are correct
% for i=1:size(time_tmp)
%     Y2=reshape(states_H(i,:)',5,5);
%     imagesc(Y2)
%     pause(0.1)
% end

% compare the forward Euler method simulation with ode45 simulation
STATE_fE=simulate_uncontrolled_forward_Euler(length(time),x_T,h,@memory_network_dynamics);

STATE_fE=STATE_fE';
for i=1:length(time)
   error_simulation(i)= norm(STATE_fE(i,:)-states_T(i,:),2)/norm(states_T(i,:),2);
end
plot(error_simulation,'k');


plot(states_T(:,1))
hold on 
plot(STATE_fE(:,1),'m')

plot(states_T(1:200,:),'r','LineWidth',1.5)


% initial state
x0=states_H(end,:)';
% desired state
xd=states_T(end,:)';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% end of computation of initial and desired states and additional tests
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% computation of the initial control sequence when all the nodes are
% controlled
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% initial control solution
initial_solution=randn(N*time_steps,1);
% control matrix
Bmatrix=eye(N,N);
% this is the control solution
solution_initial(:,1)=solve_problem_dynamic_unconstrained(time_steps,x0,xd,h,Bmatrix,initial_solution,@memory_network_dynamics)

% compute the response for the calculated control inputs
STATE_controlled=simulate_controlled_forward_Euler(time_steps,x0,h,Bmatrix,solution_initial(:,1),@memory_network_dynamics);

% compute the control error
for i=1:time_steps+1
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
if i==1
    solution_initial_selection(:,i)=[diag(Bmatrix); solution_initial(:,i)];
else
    % map back the reduced control inputs to the larger vector and add
    % zeros for the nodes that are not controlled
    for idx1=1:time_steps
        solution_initial((idx1-1)*N+nodes_that_are_controlled,i)= solution_initial_reduced_tmp((idx1-1)*no_control_nodes+1:idx1*no_control_nodes,1);
        solution_initial((idx1-1)*N+nodes_that_are_not_controlled,i)= zeros(numel(nodes_that_are_not_controlled),1);
    end
    solution_initial_selection(:,i)=[diag(Bmatrix); solution_initial(:,i) ];
end

% here we solve the main problem
[solution_selection(:,i),fval(i),exitflag(i),info{i}]=solve_problem_dynamic(time_steps,x0,xd,h,solution_initial_selection(:,i),no_control_nodes,@memory_network_dynamics);

nodes_controlled_and_not=solution_selection(1:N,i);
Bmatrix=diag(nodes_controlled_and_not);
control_input=solution_selection(N+1:end,i);

% reduce the Bmatrix and the control input to exclude the nodes that are
% not controlled - this is used to compute the control sequence for the
% next iteration of the optimization problem
nodes_that_are_controlled = find(nodes_controlled_and_not>0);
nodes_that_are_not_controlled=find(nodes_controlled_and_not==0);
Bmatrix_reduced=zeros(N,numel(nodes_that_are_controlled));

indx1=1;
for i1=1:numel(nodes_that_are_controlled)
    Bmatrix_reduced(nodes_that_are_controlled(i1),indx1)=1;
    indx1=indx1+1;    
end
control_input_reduced=[]
for i1=1:time_steps
    for i2=1:numel(nodes_that_are_controlled)
        control_input_reduced=[control_input_reduced; control_input(((i1-1)*N)+nodes_that_are_controlled(i2))]
    end;
end
clear STATE_controlled;
STATE_controlled=simulate_controlled_forward_Euler(time_steps,x0,h,Bmatrix_reduced,control_input_reduced,@memory_network_dynamics);
error_iterations(2*i-1)=norm(STATE_controlled(:,end)-xd(:,1));

clear STATE_controlled;
solution_initial_reduced_tmp=solve_problem_dynamic_unconstrained(time_steps,x0,xd,h,Bmatrix_reduced,control_input_reduced,@memory_network_dynamics);
STATE_controlled=simulate_controlled_forward_Euler(time_steps,x0,h,Bmatrix_reduced,solution_initial_reduced_tmp,@memory_network_dynamics);
error_iterations(2*i)=norm(STATE_controlled(:,end)-xd(:,1));
end

% visualize the nodes that are controlled
controlled_nodes=diag(Bmatrix);
controlled_nodes_visualize=reshape(controlled_nodes,N1,N1);
imagesc(controlled_nodes_visualize)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    end of iterative solution of the control node selection problem
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%







%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                   final computations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% compute the response for the final calculated control inputs


STATE_controlled=simulate_controlled_forward_Euler(time_steps,x0,h,Bmatrix_reduced,solution_initial_reduced_tmp,@memory_network_dynamics)
% compute the control error
for i=1:time_steps+1
    error(i)=norm(STATE_controlled(:,i)-xd(:,1));
end
plot(error)

% open loop simulation from the final state
STATE_fE=simulate_uncontrolled_forward_Euler(length(time),STATE_controlled(:,time_steps+1),h,@memory_network_dynamics);
STATE_fE=STATE_fE';

for i=1:length(time)
   error_simulation(i)= norm(STATE_fE(i,:)-states_T(i,:),2);%/norm(states_T(i,:),2);
end
plot(error_simulation);

% visualize the open-loop simulation
for i=1:size(time_tmp)
    Y2=reshape(STATE_fE(i,:)',5,5);
    imagesc(Y2)
    pause(0.5)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         end of final computations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%        verify the results by randomly selecting the control nodes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i=1:1000

ss = sort(randsample(25,no_control_nodes));
Brandom=zeros(25,no_control_nodes);
for j=1:no_control_nodes
    Brandom(ss(j),j)=1;
end
solution_random=solve_problem_dynamic_unconstrained(time_steps,x0,xd,h,Brandom,randn(time_steps*no_control_nodes,1),@memory_network_dynamics);

% compute the response for the calculated control inputs
STATE_controlled_random=simulate_controlled_forward_Euler(time_steps,x0,h,Brandom,solution_random,@memory_network_dynamics)
% compute the control error
for indx1=1:time_steps+1
    error_random(i,indx1)=norm(STATE_controlled_random(:,indx1)-xd(:,1));
end
end



hist(error_random(:,time_steps))
histogram(error_random(:,time_steps), 10, 'Normalization','probability' )
ytix = get(gca, 'YTick')
set(gca, 'YTick',ytix, 'YTickLabel',ytix*100)
hold on 
line([error(end) error(end)], [0 0.25],'color','r','linewidth',2);









hist(error_random(:,9))
hold on 
line([error_iterations(end) error_iterations(end)], [0 20],'color','r','linewidth',2);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  end verify the results by randomly selecting the control nodes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%