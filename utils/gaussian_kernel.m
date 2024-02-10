%%*************************************************************************
%% Compute the gaussian kernel matrix between X and Y
%   Input:  X - n*k matrix, Y - m*k matrix
%   Output: K - n*m matrix
%%*************************************************************************
function K = gaussian_kernel(X, Y, sigma)

if isscalar(X) && isscalar(Y)
    arg = (X - Y)^2 / sigma; 
else
    [n, ~] = size(X);
    [m, ~] = size(Y); 
    tmpX = repmat(sum(X.^2, 2), 1, m);  
    tmpY = repmat(sum(Y.^2, 2), 1, n); 
    arg = (tmpX - 2*X*Y' + tmpY') / sigma; 
end

K = exp(-arg); 

end