% function that generates the control input matrix B 
% - input parameters: 
%       - n - state order
%       - controlled_nodes - vector whose i-th entry is 1 (0) is the node i
%       is (is not) controlled
% - output parameters:
%       - Bmatrix_reduced - generated control input matrix with eliminated
%       control that are not active
%       
% - Author: Aleksandar Haber
% December 2019 - February 2020
% Note: This function only works with dynamics whose local  subsystem  states are 
% [\dot{x}_{i} ; x_{i} ]
% that is, if the i-th subsystem is actuated, then i-th block diagonal
% element of the matrix B is [0 ; 1], otherwise it is [0 ; 0]

function Bmatrix_reduced=formBmatrix_3_reduced(n,controlled_nodes)
only_controlled_nodes=find(controlled_nodes>0);
Bmatrix_reduced=zeros(n,numel(only_controlled_nodes));
for i = 1:numel(only_controlled_nodes)
    Bmatrix_reduced(2*only_controlled_nodes(i),i)=1; 
end
end