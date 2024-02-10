%%************************************************************************************
%% Solve the conic optimization model using an extragradient method
%   Input:  data - 1*1 structure, reg1 - 1*1 scalar, reg2 - 1*1 scalar
%   Output: gamma - m*1 vector, t - 1*1 scalar
%%************************************************************************************
function [gamma, c, t, res_time, res_norm] = EG(data, reg1, reg2, verbose)

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
nIter = 3000; 

gamma = ones(m, 1) / m; 
X = ones(m) / m; 
eta = 0.05; 
[r_gamma, r_X] = residue(gamma, X, Phi, Q, z, reg1, reg2);
mu = norm(r_gamma) + norm(r_X, 'fro'); 
res_time = [0]; 
res_norm = [mu]; 

if verbose
    fprintf('\n-------------- EG ---------------\n');
    fprintf('iter |  cost  |  residue  |  time\n');
end

tstart = clock;

%% main loop
for iter = 1:nIter
    
    % compute the residue function
    [r_gamma, r_X] = residue(gamma, X, Phi, Q, z, reg1, reg2);
    mu = norm(r_gamma) + norm(r_X, 'fro'); 

    % compute the step
    [g_gamma, g_X] = gradient(gamma, X, Phi, Q, z, reg1, reg2); 
    gamma_hat = gamma - eta*g_gamma; 
    X_hat = X - eta*g_X;

    % compute the extra step
    [g_gamma_hat, g_X_hat] = gradient(gamma_hat, X_hat, Phi, Q, z, reg1, reg2); 
    gamma = gamma - eta*g_gamma_hat; 
    X = X - eta*g_X_hat;

    if mu < 1e-8
        c = kernel_cost(gamma, data, reg2);
        t = etime(clock, tstart); 
        res_time = [res_time; t]; 
        res_norm = [res_norm; mu];
        break; 
    end

    if mod(iter, 300) == 0
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