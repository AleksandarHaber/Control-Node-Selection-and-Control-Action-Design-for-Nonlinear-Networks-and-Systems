%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% - this code creates and stores MATLAB functions that describe the
% dynamics and gradient of the dynamics of the Duffing network
% 
% - the dynamics is stored in "duffing_network_dynamics.m"
% - the gradient is stored in "duffing_network_dynamics_gradient.m"
% - the network adjacency matrix is a random geometric graph 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear

% the largest size of the network, before the largest connected component
% has been taken
nsize=12;
% is the final network size
N=0;
% this is the minimal acceptable network size
N_min=9;
% this code generates the NxN adjacency matrix denoted by "Atmp"
while (N<N_min)
Atmp=geo(nsize);
    % ensure that it is connected
    [ci sizs] = components(Atmp);
    k = find(ci==1);
    Atmp = Atmp(k,k);
    Atmp=Atmp+speye(size(Atmp));
    % Cuthill-Mckee permutation
    r=symrcm(Atmp);
    I=speye(size(Atmp));
    P=I(r,:);
    Apattern=P*Atmp*P';
    [N,N]=size(Atmp);
end

x_variables=sym('x',[2*N,1]);

% k1 and k2 are the constants of the Duffing equations.
% spring forces: k1*x-k2*x^3
% k3 is the damping constant

k1=(10)*ones(N,N)+(10)*rand(N,N); k1=(k1+k1')/2;
k2=(1)*ones(N,N)+(1)*rand(N,N); k2=(k2+k2')/2;
k3=(1)*ones(N,N)+(1)*rand(N,N); k3=(k3+k3')/2;


% this is used to verify the code accuracy
% N=3
% Atmp=[1 1 0; 1 1 1; 0 1 1]
% x_variables=sym('x',[2*N,1]);
% k1=ones(N,N); k2=ones(N,N); k3=ones(N,N);


 for i=1:N
     dx((i-1)*2+1,1)= x_variables((i-1)*2+2);
     dx((i-1)*2+2,1)= 0;
     tmp=find(Atmp(i,:)>0);
     
     for s=1:numel(tmp)
         j=tmp(s);
         if i==j
             dx((i-1)*2+2,1)=dx((i-1)*2+2,1)+(-k1(j,j)*x_variables((j-1)*2+1)+k2(j,j)*x_variables((j-1)*2+1)^3-k3(j,j)*x_variables((j-1)*2+2));
         else
             dx((i-1)*2+2,1)=dx((i-1)*2+2,1)+(-k1(i,j)*(x_variables((i-1)*2+1)-x_variables((j-1)*2+1))+k2(i,j)*(x_variables((i-1)*2+1)-x_variables((j-1)*2+1))^3-k3(i,j)*(x_variables((i-1)*2+2)-x_variables((j-1)*2+2)));
         end
     end
    
     dx((i-1)*2+2,1)=simplify(dx((i-1)*2+2,1));
 
 end
 
gradient_dx=jacobian(dx)';

% here we store the variables
GG=matlabFunction(dx,'file','duffing_network_dynamics','vars',{x_variables});
GG2=matlabFunction(gradient_dx,'file','duffing_network_dynamics_gradient','vars',{x_variables});