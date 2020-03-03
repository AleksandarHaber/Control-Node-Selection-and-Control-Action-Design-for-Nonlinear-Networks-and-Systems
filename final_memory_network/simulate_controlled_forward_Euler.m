% - function that simulates the uncontrolled dynamics using the feedforward
%   Euler method 
% - input parameters: 
%                   - time_steps    - discrete-time simulation time 
%                   - x0            - initial state
%                   - h             - discretization constant
%                   - B             - control input matrix
%                   - control       - control sequence
%                   - fcnHandle     - function handle that describes the
%                   system dynamics
%                   
% - output parameters:
%                   - STATE         - state trajectory
% Author: Aleksandar Haber
% December 2019 - February 2020

function STATE=simulate_controlled_forward_Euler(time_steps,x0,h,B,control,fcnHandle)   
[N,~]=size(x0);
[~,no_controlled_nodes]=size(B);
STATE=zeros(N,time_steps+1);
for o=1:time_steps
        if o==1
           STATE(:,o)=x0; 
           STATE(:,o+1)=x0+h*(fcnHandle(x0)+B*control((o-1)*no_controlled_nodes+1:o*no_controlled_nodes,1));
        else 
           STATE(:,o+1)=STATE(:,o)+h*(fcnHandle(STATE(:,o))+B*control((o-1)*no_controlled_nodes+1:(o)*no_controlled_nodes,1));
        end
end
end