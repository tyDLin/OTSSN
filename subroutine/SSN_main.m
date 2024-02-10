%%************************************************************************************
%% Compute the SSN step
%   Input:  r_gamma - m*1 vector, r_X - m*m matrix
%           gamma - m*1 vector, X - m*m matrix, mu - 1*1 scalar
%           Q - m*m matrix, z - m*1 vector
%           Phi - m*m matrix, reg1 - 1*1 scalar, reg2 - 1*1 scalar         
%   Output: d_gamma - m*1 vector, d_X - m*m matrix
%%************************************************************************************
function [d_gamma, d_X] = SSN_main(r_gamma, r_X, gamma, X, mu, Q, Phi, reg1, reg2)

m = length(gamma);

%% the first step.
Z = X - (Phi * diag(gamma) * Phi' + reg1*eye(m));
[P, Sigma] = eig(Z);
alpha = find(diag(Sigma) > 0); 
beta = find(diag(Sigma) <= 0);
Omega = zeros(m); 
Omega(alpha, alpha) = ones(length(alpha)); 
sigma = diag(Sigma); 
eta = 1 - repmat(sigma(beta)', length(alpha), 1)./repmat(sigma(alpha), 1, length(beta)); 
eta = 1./eta; 
Omega(alpha, beta) = eta;  
Omega(beta, alpha) = eta'; 
L = Omega./(mu+1-Omega); 

T = r_X + P*(L.*(P'*r_X*P))*P'; 
H = Phi'*T*Phi; 
d_gamma = - r_gamma - diag(H)/(1+mu);  
d_X = -r_X; 

%% the second step (CG).
y = d_gamma; 
K = P'*Phi; 
H = K'*(L.*(K*diag(y)*K'))*K; 
r = d_gamma - ((0.5/reg2)*Q*y + mu*y + diag(H)); 
p = r;
rr = r'*r; 
for i = 1:min(m/5, 50)
    H = K'*(L.*(K * diag(p) * K'))*K;
    Ap = (0.5/reg2)*Q*p + mu*p + diag(H); 
    ss1 = rr/(p'*Ap); 
    y = y + ss1*p; 
    r = r - ss1*Ap;
    if norm(r) < 1e-6
        break;
    end
    ss2 = r'*r/rr;
    p = r + ss2*p;
    H = K'*(L.*(K * diag(y) * K'))*K;
    r = d_gamma - ((0.5/reg2)*Q*y + mu*y + diag(H)); 
    rr = r'*r; 
end
d_gamma = y; 
d_X = (d_X + P*(L.*(P'*d_X*P))*P')/(1+mu); 

%% the third step.
d_X = d_X - (P*(L.*(K * diag(d_gamma) * K'))*P'); 

end