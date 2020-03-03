% - function that simulates the uncontrolled dynamics using the ode45 or ode23s MATLAB
%   functions
% - input parameters: 
%                   - time      - simulation time 
%                   - x0        - initial state
%                   - fcnHandle - function handle that describes the system
%                   dynamics
% - output parameters:
%                   - ts - time for the state trajectory
%                   - xs - state trajectory
% Author: Aleksandar Haber 
% December 2019 - February 2020
function [ts,xs]=simulate_uncontrolled_ode45(time,x0,fcnHandle)

%[ts,xs] = ode45(@simulate_dynamics_open_loop,time,x0);
[ts,xs] = ode23s(@simulate_dynamics_open_loop,time,x0);

function dx = simulate_dynamics_open_loop(t,x)
    dx=fcnHandle(x);
end

end