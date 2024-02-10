%%************************************************************************************
%% Compute the operator
%   Input:  gamma - m*1 vector, Phi - m*m matrix, Q - m*m matrix
%           z - m*1 vector, reg1 - 1*1 scalar, reg2 - 1*1 scalar
%   Output: g - m*1 vector
%%************************************************************************************
function g = gradient_AGD(gamma, Phi, Q, z, reg1, reg2)

T = -Phi*diag(gamma)*Phi'; 
[V, D] = eig(T); 
T_new = V*max(D, 0)*V';
H = Phi'*(T-T_new)*Phi / reg1; 
g = (Q * gamma - z) / (2 * reg2) - diag(H); 

end