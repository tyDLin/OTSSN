%%************************************************************************************
%% Compute the operator
%   Input:  gamma - m*1 vector, Phi - m*m matrix, Q - m*m matrix
%           z - m*1 vector, reg1 - 1*1 scalar, reg2 - 1*1 scalar
%   Output: grad - m*1 vector, hess - m*m matrix
%%************************************************************************************
function [g_gamma, g_X] = gradient(gamma, X, Phi, Q, z, reg1, reg2)

m = length(z);
H = Phi'*X*Phi; 
g_gamma = (Q * gamma - z) / (2 * reg2) - diag(H); 
g_X = Phi * diag(gamma) * Phi' + reg1*eye(m); 

end