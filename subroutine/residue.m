%%************************************************************************************
%% Compute the residue
%   Input:  gamma - m*1 vector, Phi - m*m matrix, Q - m*m matrix
%           z - m*1 vector, reg1 - 1*1 scalar, reg2 - 1*1 scalar
%   Output: r_gamma - m*1 vector, r_X - m*m matrix
%%************************************************************************************
function [r_gamma, r_X] = residue(gamma, X, Phi, Q, z, reg1, reg2)

[r_gamma, g_X] = gradient(gamma, X, Phi, Q, z, reg1, reg2); 
X_new = X - g_X;
[V, D] = eig(X_new); 
X_new = V*max(D, 0)*V';
r_X = X - X_new; 

end