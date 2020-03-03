% function that generates the control input matrix B 
% - input parameters: 
%       - n - state order
%       - controlled_nodes - vector whose i-th entry is 1 (0) is the node i
%       is (is not) controlled
% - output parameters:
%       - Bmatrix - generated control input matrix
% - Author: Aleksandar Haber
% December 2019 - February 2020
% Note: This function only works with dynamics whose local  subsystem  states are 
% [\dot{x}_{i} ; x_{i} ]
% that is, if the i-th subsystem is actuated, then i-th block diagonal
% element of the matrix B is [0 ; 1], otherwise it is [0 ; 0]

function Bmatrix=formBmatrix_3(n,controlled_nodes)
Bmatrix=zeros(n,numel(controlled_nodes));
for i = 1:numel(controlled_nodes)
    Bmatrix(2*i,i)=controlled_nodes(i); 
end
end