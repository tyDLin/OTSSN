%%************************************************************************************
%% Solve the conic optimization model using a semi-smooth Newton method
%   Input:  data - 1*1 structure, reg1 - 1*1 scalar, reg2 - 1*1 scalar
%   Output: gamma - m*1 vector, t - 1*1 scalar
%%************************************************************************************
function [gamma, c, t, res_time, res_norm] = SSN(data, reg1, reg2, verbose)

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
X = ones(m) / (m*m); 
kappa = 1.0; 
[r_gamma, r_X] = residue(gamma, X, Phi, Q, z, reg1, reg2);
mu = norm(r_gamma) + norm(r_X, 'fro'); 
res_time = [0]; 
res_norm = [mu]; 

if verbose
    fprintf('\n-------------- SSNEG ---------------\n');
    fprintf('iter |  cost  |  residue  |  time\n');
end

tstart = clock;

%% main loop
for iter = 1:nIter
    
    % compute the residue function
    mu = norm(r_gamma) + norm(r_X, 'fro'); 

    % compute SSN step
    [d_gamma, d_X] = SSN_main(r_gamma, r_X, gamma, X, (m/5)*kappa*mu, Q, Phi, reg1, reg2); 
    
    % compute the next iterate
    gamma = gamma + d_gamma; 
    X = X + d_X;

    % update the parameter kappa. 
    [r_gamma, r_X] = residue(gamma, X, Phi, Q, z, reg1, reg2);
    rho = -(r_gamma'*d_gamma + trace(r_X'*d_X))/(norm(d_gamma)^2 + norm(d_X, 'fro')^2); 
    if rho >= 1
        kappa = max(0.5*kappa, 1e-16);
    elseif rho >= 1e-6
        kappa = 1.2*kappa;
    else
        kappa = 25*kappa; 
    end

    if mu < 1e-8 % 5e-3
        c = kernel_cost(gamma, data, reg2);
        t = etime(clock, tstart);
        res_time = [res_time; t]; 
        res_norm = [res_norm; mu];
        if verbose
            fprintf('%5.0f|%3.2e|%3.2e|%3.2e\n', iter, c, mu, t);
        end
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