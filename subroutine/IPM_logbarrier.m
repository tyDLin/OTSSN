%%************************************************************************************
%% Compute the gradient and the Hessian of the barrier function
%   Input:  gamma - m*1 vector, Phi - m*m matrix, reg - 1*1 scalar
%   Output: grad - m*1 vector, hess - m*m matrix
%%************************************************************************************
function [grad, hess] = IPM_logbarrier(gamma, Phi, reg)

m = length(gamma); 
H = Phi' * inv(Phi * diag(gamma) * Phi' + reg * eye(m)) * Phi; 
H = 0.5*(H + H'); 
grad = diag(H); 
hess = -H.^2; 

end