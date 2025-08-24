clc; clear all; close all;

%% ------------------ Continuous-time aircraft model ------------------

% System matrices (linearized longitudinal dynamics)
Ac = [-1.2822 0 0.98 0;
      0 0 1 0;
     -5.4293 0 -1.8366 0;
    -128.2 128.2 0 0];

% Input matrix (elevator deflection)
Bc = [-0.3; 0; -17; 0];

% Output matrix (outputs: pitch angle, altitude)
Cc = [0 1 0 0;
      0 0 0 1];

% Feedthrough matrix
D = [0; 0];

% Sampling time for discretization
Ts = 0.25;  

% Continuous-time state-space system
sys_c = ss(Ac, Bc, Cc, D);

% Discretized system
sys_d = c2d(sys_c, Ts);

A = sys_d.A;  % Discrete-time state matrix
B = sys_d.B;  % Discrete-time input matrix
C = sys_d.C;  % Discrete-time output matrix

%% ---------------Simulation parameters-------------------

T = 20;                  
n = size(A,1);
m = size(B,2);
q = size(C,1);
Ns = T/Ts;               
N = 10;             % Prediction horizon

%% -------------------Constraints--------------------------

u_max = 0.262; u_min = -0.262;              % Input (elevator angle) limits (rad)

du_max = 0.524*Ts; du_min = -0.524*Ts;      % Input rate limits (rad/sec scaled by Ts)

y1_min = -0.349;  y1_max = 0.349;           % Output constraints

y2_min = 0;       y2_max = 15545;


y_min = [y1_min; y2_min];
y_max = [y1_max; y2_max];

%% --------------Cost function weights-------------

Q = eye(n);
R = 10 * eye(m);

% Solve discrete-time LQR problem to obtain terminal weight P
[Klqr, P, ~] = dlqr(A, B, Q, R);  

% Initial state (small pitch, altitude 10 m)
x0 = [0;0;0;10];

%% ------------------ MPC matrices ------------------

% Function builds all MPC matrices (Hessian, gradient, constraints)
% Outputs:
%   G - Hessian of QP
%   F - Gradient matrix
%   J,W,c - Matrices for state/output constraints
%   E_total, b_total - Input and input rate constraints

[G,F,J,W,c,E_total,b_total] = MPC_Matrices(A,B,C,Q,R,N,P,u_min,u_max,du_min,du_max,y_min,y_max);

%% ------------------ MPC simulation (with constraints) ------------------

Xhh = [];
Yhh = [];
Uhh = [];

for k = 1:Ns
    % Solve quadratic program at each time step
    [Uu, ~, exitflag] = quadprog(G, F*x0, [J; E_total], [c + W*x0; b_total]);
    if exitflag ~= 1
        warning('QP doesnt solve in %d ', k);
        break;
    end
    % Apply first optimal control input
    uff = Uu(1:m);  
    Uhh = [Uhh, uff];

    % State update
    x1 = A * x0 + B * uff;
    Xhh = [Xhh, x1];

    % Output update
    yy = C * x0 + D * uff;
    Yhh = [Yhh, yy];

    % Update state for next iteration
    x0 = x1;
end

% Plot closed-loop trajectories
figure(1);
subplot(3,1,1);
plot(Yhh(1,:), 'LineWidth',1.5);
ylabel('y_1 (pitch angle)');

subplot(3,1,2);
plot(Yhh(2,:), 'LineWidth',1.5);
ylabel('y_2 (altitude)');

subplot(3,1,3);
stairs(Uhh, 'LineWidth',1.5);
ylabel('u (elevator)');
xlabel('Time step');

%% ------------------ MOAS computation ------------------

% Define ranges for state-space grid search

x2range = y1_min:0.01:y1_max;
x4range = y2_min:100:y2_max;

% Compute Maximal Output Admissible Set (MOAS)
% Inputs:
%   A,C - system matrices
%   y_min,y_max - output bounds
%   x2range,x4range - grid ranges for state variables
%   10 - time horizon for feasibility check

moas = moasApprox(A,C,y_min,y_max,x2range,x4range,10);


% Plot MOAS approximation
figure(2);
plot(moas(:,1), moas(:,2), 'b.');
xlabel('pitch angle (rad)');
ylabel('altitude (m)');
title('Approximation of MOAS');
grid on;

