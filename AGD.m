%%************************************************************************************
%% Solve the conic optimization model using an extragradient method
%   Input:  data - 1*1 structure, reg1 - 1*1 scalar, reg2 - 1*1 scalar
%   Output: gamma - m*1 vector, t - 1*1 scalar
%%************************************************************************************
function [gamma, c, t, res_time, res_norm] = AGD(data, reg1, reg2, verbose)

%% input data
M = data.M; 
Phi = data.Phi;  
KX1 = data.KX1; 
KY1 = data.KY1; 
KX2 = data.KX2; 
KY2 = data.KY2; 

%% initialization
m = length(M); 
Q = KX1 + KY1; 
z = mean(KX2)' + mean(KY2)' - 2 * reg2 * M; 
nIter = 300; 

gamma = ones(m, 1) / m; 
gamma_hat = gamma; 
eta = 0.01; 
lambda = 0;
g = gradient_AGD(gamma, Phi, Q, z, reg1, reg2);
mu = norm(g); 
res_time = [0]; 
res_norm = [mu]; 

if verbose
    fprintf('\n-------------- AGD ---------------\n');
    fprintf('iter |  cost  |  residue  |  time\n');
end

tstart = clock;

%% main loop
for iter = 1:nIter

    % compute the step
    g_hat = gradient_AGD(gamma_hat, Phi, Q, z, reg1, reg2);
    gamma_next = gamma_hat - eta*g_hat; 

    % add the momentum
    lambda_next = (1+sqrt(1+4*lambda*lambda))/2; 
    beta = (1-lambda)/lambda_next; 
    gamma_hat = (1-beta)*gamma_next + beta*gamma;
    
    % update the next iterate
    gamma = gamma_next; 

    % compute the residue function
    g = gradient_AGD(gamma, Phi, Q, z, reg1, reg2);
    mu = norm(g); 

    if mu < 1e-8
        c = kernel_cost(gamma, data, reg2);
        t = etime(clock, tstart); 
        res_time = [res_time; t]; 
        res_norm = [res_norm; mu];
        break; 
    end

    if mod(iter, 30) == 0
        c = kernel_cost(gamma, data, reg2);
        t = etime(clock, tstart);
        res_time = [res_time; t]; 
        res_norm = [res_norm; mu];
        if verbose
            fprintf('%5.0f|%3.2e|%3.2e|%3.2e\n', iter, c, mu, t);
        end
    end

    c = kernel_cost(gamma, data, reg2);
    t = etime(clock, tstart); 
end

end