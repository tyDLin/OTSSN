%%************************************************************************************
%% Compute all the gaussian kernel matrices between X, Y, X_f and Y_f
%   Input:  X - n*k matrix, Y - n*k matrix, X_f - m*k matrix, Y_f - m*k matrix
%   Output: Phi - m*m matrix, M - m*1 vector
%           KX1 - m*m matrix, KY1 - m*m matrix
%           KX2 - n*m matrix, KY2 - n*m matrix
%           KX3 - n*n matrix, KY3 - n*n matrix
%%************************************************************************************
function data = kernel_make(X, Y, X_f, Y_f, sigma)

XY_f = [X_f Y_f]; 

%% compute the distance matrix
M = 0.5*sum((X_f - Y_f).^2, 2); 

%% compute the kernel matrices

KX1 = gaussian_kernel(X_f, X_f, sigma); 
KY1 = gaussian_kernel(Y_f, Y_f, sigma); 
KX2 = gaussian_kernel(X, X_f, sigma); 
KY2 = gaussian_kernel(Y, Y_f, sigma); 
KX3 = gaussian_kernel(X, X, sigma); 
KY3 = gaussian_kernel(Y, Y, sigma);
K   = gaussian_kernel(XY_f, XY_f, sigma); 

%% compute the Cholesky factorization
Phi = chol(K); 

%% output
data.Phi = Phi; 
data.M = M; 
data.KX1 = KX1; 
data.KY1 = KY1; 
data.KX2 = KX2; 
data.KY2 = KY2; 
data.KX3 = KX3; 
data.KY3 = KY3; 

end