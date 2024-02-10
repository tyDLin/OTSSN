%%************************************************************************************
%% Compute the kernel SoS estimator for OT cost
%   Input:  gamma - m*1 vector, reg - 1*1 scalar
%           KX2 - n*m matrix, KY2 - n*m matrix
%           KX3 - n*n matrix, KY3 - n*n matrix
%   Output: c - 1*1 scalar
%%************************************************************************************
function c = kernel_cost(gamma, data, reg)

%% input data
KX2 = data.KX2; 
KY2 = data.KY2; 
KX3 = data.KX3; 
KY3 = data.KY3; 

tmp1 = mean(KX3(:)) + mean(KY3(:)); 
tmp2 = (mean(KX2) + mean(KY2)) * gamma;
c = (tmp1 - tmp2) / (2 * reg); 

end