%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% this code creates and stores MATLAB functions that describe the
% dynamics and gradient of the memory network
% the dynamics is stored in "memory_network_dynamics.m"
% the gradient is stored in "memory_network_dynamics_gradient"
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% The size of the network is N=N1xN1 
N1=5; N=N1*N1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                   model and letter definitions
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


% the parameter p defines how many letters will be stored, 
% we store ksi{1},ksi{2},ksi{3},\ldots, ksi{p}
p=3;
ksi{1}=H(:);ksi{2}=T(:);ksi{3}=L(:); ksi{4}=E(:); ksi{5}=F(:); ksi{6}=S(:); ksi{7}=O(:); ksi{8}=Y(:); ksi{9}=X(:);

% This is the interconnection matrix that stores the learned patterns
A=zeros(N,N);
for i=1:N
    for j=1:N
       sum=0;
        for s=1:p
          sum=sum+ksi{s}(i)*ksi{s}(j);
        end
       sum=sum/N;
       A(i,j)=sum; 
    end
end

% epsilon parameter
epsi=0.8;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                   end of the model definition
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


x_variables=sym('x',[N,1]);

 for i=1:N
         sum1=0;
         sum2=0;
             for l=1:N
                sum1=sum1+A(i,l)*sin(x_variables(l)-x_variables(i));
                sum2=sum2+sin(2*(x_variables(l)-x_variables(i)));
             end
                dx(i,1)=sum1+(epsi/N)*sum2;
end

GG=matlabFunction(dx,'file','memory_network_dynamics','vars',{x_variables});
gradient_dx=jacobian(dx)';
GG2=matlabFunction(gradient_dx,'file','memory_network_dynamics_gradient','vars',{x_variables});