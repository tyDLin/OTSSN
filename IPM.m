%%************************************************************************************
%% Solve the conic optimization model using an interior-point method
%   Input:  data - 1*1 structure, reg1 - 1*1 scalar, reg2 - 1*1 scalar
%   Output: gamma - m*1 vector, t - 1*1 scalar
%%************************************************************************************
function [gamma, c, t] = IPM(data, reg1, reg2, verbose)

%% input data
M = data.M; 
Phi = data.Phi;  
KX1 = data.KX1; 
KY1 = data.KY1; 
KX2 = data.KX2; 
KY2 = data.KY2; 

%% initialization
m = length(M); 
gamma = ones(m, 1) / m; 
Q = KX1 + KY1; 
z = mean(KX2)' + mean(KY2)' - 2 * reg2 * M; 
eta = m; 
nIter = 5000; 

if verbose
    fprintf('\n-------------- IPM ---------------\n');
    fprintf('iter |  cost  |  eta  |  time\n');
end

tstart = clock;

%% main loop
for iter = 1:nIter

    % compute the gradient and Hessian of log-barrier functions
    [g_barrier, H_barrier] = IPM_logbarrier(gamma, Phi, reg1); 

    % form the Hessian
    hess = - eta * H_barrier + Q / (2 * reg2); 
    hess = 0.5 * (hess + hess'); 
    inv_hess = inv(hess+1e-12); 
    
    % symmetrize inverse Hessian (numerical stability)
    inv_hess = 0.5 * (inv_hess + inv_hess');  
    
    % compute the gradient
    g = (Q * gamma - z) / (2 * reg2) - eta * g_barrier; 
    
    % compute the next iterate
    d = inv_hess * g; 
    decr = sqrt(g'*d); 
    gamma_new = gamma - 0.8*d / (1 + decr / sqrt(eta));  

    % check that nothing has blown up
    if isnan(decr)
        fprintf('Issue: NaN decrement is detected\n');
    else
        gamma = gamma_new; 
    end
    
    % check decrement and decrease eta if small than certain threshold
    if decr < 1e-8
        eta = 0.8 * eta; 
    end

    if eta < 1e-8
        c = kernel_cost(gamma, data, reg2);
        t = etime(clock, tstart); 
        break; 
    end

    if verbose && mod(iter, 100) == 1
        c = kernel_cost(gamma, data, reg2);
        t = etime(clock, tstart); 
        fprintf('%5.0f|%3.2e|%3.2e|%3.2e\n', iter, c, eta, t);
    end  
end

if iter == nIter
    c = kernel_cost(gamma, data, reg2);
    t = etime(clock, tstart); 
end

end