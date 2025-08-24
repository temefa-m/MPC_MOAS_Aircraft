function [G,F,J,W,c,E_total,b_total] = MPC_Matrices(A,B,C,Q,R,N,P,u_min,u_max,du_min,du_max,y_min,y_max)

% -------------------------------------------------------------------------
% Function to build MPC cost function matrices and constraint matrices
% Inputs:
%   A, B, C    : State-space model matrices (discrete-time)
%   Q, R       : State/output weighting and input weighting matrices
%   N          : Prediction horizon
%   P          : Terminal cost matrix
%   u_min/u_max: Input bounds
%   du_min/du_max: Input rate-of-change bounds
%   y_min/y_max: Output bounds
%
% Outputs:
%   G, F       : Cost function matrices for quadprog (1/2 * U'GU + x0'FU)
%   J, W, c    : Constraint matrices (J*U <= W*x0 + c)
%   E_total,
%   b_total    : Input and delta-input inequality constraints
% -------------------------------------------------------------------------



n = size(A,1);      % number of states
m = size(B,2);      % number of inputs

%% --- Build Gamma (input response matrix) and Phi (state response matrix)

Gamma = [];
for i = 1:N
    Rw = [];
    for j = 1:i
        % Contribution of input at step j to state at step i
        Rw = [(A^(j-1))*B, Rw];
    end
    % Pad with zeros so that Rw has size n x (N*m)
    Rw = [Rw , zeros(n, N*m - size(Rw,2))];
    Gamma = [Gamma; Rw];
end

% Phi collects the free evolution of the states
Phi = [];
for i = 1:N
    Phi = [Phi; A^i];
end

%% --- Cost function matrices
Psi = kron(eye(N), R);
Omega = kron(eye(N-1), Q);
Omega = blkdiag(Omega, P);

% Quadratic cost matrix
G = 2 * (Psi + Gamma' * Omega * Gamma);
G = (G+G')/2;  % Enforce symmetry for numerical stability

% Linear cost matrix
F = 2 * Gamma' * Omega * Phi;

%% --- Constraints setup
% Basic bounds (inputs and outputs)
bi = [-u_min; u_max; -y_min; y_max];
bN = [-y_min; y_max]; % terminal step output bounds
c  = [kron(ones(N,1), bi); bN]; % Collect constraints into c

% Build stage constraint matrices (Mi for intermediate steps, MN for terminal step)
Mi = [ -eye(m), zeros(m, n-m);   % -u <= u_min
        eye(m), zeros(m, n-m);   %  u <= u_max
       -C;                       % -y <= y_min
        C ];                     %  y <= y_max

MN = [-C; C];
% Block-diagonal constraint structure
Dd = [Mi; zeros(size(c,1)-size(Mi,1), n)];
M  = [blkdiag(kron(eye(N-1), Mi), MN)];
M  = [zeros(size(c,1)-size(M,1), size(M,2)); M];
%% --- Delta input (rate of change) constraints
Ed = zeros((N-1)*m, N*m);
for i = 1:N-1
    % Enforce Δu = u(k+1)-u(k)
    Ed(i, (i-1)*m+1:i*m) = -1;
    Ed(i, i*m+1:(i+1)*m) =  1;
end
Ed = [Ed; -Ed]; % Add symmetric bounds
bd = repmat([du_max; -du_min], N-1, 1);
%% --- Input bounds
Eu = kron(eye(N), [eye(m); -eye(m)]);
bu = repmat([u_max; -u_min], N, 1);
%% --- Combine all constraints
E_total = [Eu; Ed];
b_total = [bu; bd];

% Final condensed inequality constraint form: J*U <= W*x0 + c
eps = [E_total; zeros(size(c,1)-size(E_total,1), m*N)];
J   = M*Gamma + eps;
W   = -(Dd + M*Phi);

end
