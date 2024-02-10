%%************************************************************************************
%% Generate the sample k-dimensional data X and Y
%   Input:  n - the number of samples/filling points
%   Output: X - n*k matrix, Y - n*k matrix, X_f - n*k matrix, Y_f - n*k matrix
%%************************************************************************************
function [X, Y, X_fill, Y_fill] = generate_data(nsamples, k)

%% initialization
mu1 = rand(2, k); 
mu2 = rand(3, k); 

t1 = [0.4, 0.6]; 
t2 = [0.2, 0.2, 0.6];

%% generate data samples
u1 = rand(nsamples, 1); 
u2 = rand(nsamples, 1); 

X = randn(nsamples, k)*0.09; 
Y = randn(nsamples, k)*0.075; 

ind1 = find(u1 < t1(1));
ind2 = find(u1 >= t1(1));
X(ind1, :) = X(ind1, :) + repmat(mu1(1, :), length(ind1), 1); 
X(ind2, :) = X(ind2, :) + repmat(mu1(2, :), length(ind2), 1); 

ind1 = find(u2 < t2(1));
ind2 = find(u2 >= t2(1) & u2 < t2(1) + t2(2));
ind3 = find(u2 >= t2(1) + t2(2));
Y(ind1, :) = Y(ind1, :) + repmat(mu2(1, :), length(ind1), 1); 
Y(ind2, :) = Y(ind2, :) + repmat(mu2(2, :), length(ind2), 1); 
Y(ind3, :) = Y(ind3, :) + repmat(mu2(3, :), length(ind3), 1); 

%% generate filling points
p = sobolset(2*k, 'Skip', 3000); 
sob = net(p, nsamples); 
X_fill = sob(:, 1:k); 
Y_fill = sob(:, k+1:end); 

end