function moas = moasApprox(A,C,y_min,y_max,x2range,x4range,tmax)
% moasApprox  Approximation of the Maximal Output Admissible Set (MOAS)
%
%   moas = moasApprox(A,C,y_min,y_max,x2range,x4range,tmax)
%
%   This function computes an approximation of the MOAS for a linear
%   discrete-time system:
%       x(k+1) = A x(k)
%       y(k)   = C x(k)
%
%   The method checks output constraint satisfaction over a finite horizon.
%
%   Inputs:
%       A        - State transition matrix (n x n)
%       C        - Output matrix (q x n)
%       y_min    - Lower bounds on outputs (q x 1)
%       y_max    - Upper bounds on outputs (q x 1)
%       x2range  - Range of pitch angle values (for grid search)
%       x4range  - Range of altitude values (for grid search)
%       tmax     - Maximum horizon for constraint checking
%
%   Output:
%       moas     - Set of (x2, x4) pairs that satisfy constraints
%
%   The algorithm loops over candidate states defined by x2 (pitch angle)
%   and x4 (altitude), forms inequalities for outputs y = C*A^t*x0,
%   and checks whether all constraints are satisfied up to time tmax.
%
%   Example:
%       moas = moasApprox(A,C,y_min,y_max,-0.35:0.01:0.35,0:100:15545,10);
moas = [];
for x20 = x2range
    for x40 = x4range
        % Candidate initial state: vary only x2 and x4
        x0 = [0;x20;0;x40];

        % Initialize constraint matrices
        O = [];
        bx = [];
        % Build constraints for each time step up to tmax
        for t = 1:tmax
            O  = [O;  C*(A^t)];     % y <= y_max
            O  = [O; -C*(A^t)];     % -y <= -y_min
            bx = [bx; y_max];
            bx = [bx; -y_min];
        end

        % Check feasibility of candidate state
        if all(O*x0 <= bx)
            moas = [moas; x20, x40]; % Add feasible state to MOAS
        end
    end
end

figure;
plot(moas(:,1), moas(:,2), 'b.');
xlabel('pitch angle (rad)');
ylabel('altitude (m)');
title('Approximation of MOAS');
grid on;

end
