clc
clear

A1=[1 0 2 0; 0 1 0 2; 0 0 3 0; 0 0 0 3];
B1=[2 0;0 2;3 0;0 3];
x01=[-10;10;-1;1];

A2=[1 0 3 0; 0 1 0 3; 0 0 7 0; 0 0 0 7];
B2=[3 0; 0 3; 7 0; 0 7];
x02=[10;10;1;1];

A3=[1 0 1 0; 0 1 0 1; 0 0 1.1 0; 0 0 0 1.1];
B3=[1 0; 0 1; 1.1 0; 0 1.1];
x03=[10;-10;1;-1];

A4=[1 0 6 0; 0 1 0 6; 0 0 20 0; 0 0 0 20];
B4=[6 0;0 6;20 0; 0 20];
x04=[-10;-10;-1;-1];

%%
% Define variables and parameters
% Tfinal = 5;   % Horizon length
umax = 60;    % Control limit

% Define the parameters and variables
N = 4;                % Number of aircraft
Tfinal = 5;           % Time horizon
T = Tfinal + 1;       % Number of time steps

% Define system matrices A and B for each aircraft
A = {A1, A2, A3, A4};
B = {B1, B2, B3, B4};

% Set up the optimization variables
X = sdpvar(4, Tfinal+1);     % Position variables
U = sdpvar(4, Tfinal);       % Input variables

%%
% Define the objective function
objective = sum(sum(X.^2)) + sum(sum(U.^2));

% Initialize the constraints
constraints = [];

% Build the constraints for each aircraft
for i = 1:N
    equalityConstraint = [];
    inequalityConstraint = [];

    for t = 1:Tfinal-1 % Adjusted loop index

        AX = A{i} * X(:, t);
        BU = B{i} * U(:, t);
        equalityConstraint = [equalityConstraint, X(i, t+1) == AX + BU];
        
        inequalityConstraint = [inequalityConstraint, abs(U(i, t)) <= umax / Tfinal];
    end

    % Add final state constraint
    equalityConstraint = [equalityConstraint, X(i, Tfinal) == x_f(i)];

    constraints = [constraints, equalityConstraint, inequalityConstraint];
end

% Set up the optimization problem
ops = sdpsettings('solver', 'mosek', 'verbose', 2);

% Solve the optimization problem
optimize(constraints, objective, ops);

% Extract the optimal solutions
X_opt = value(X);
U_opt = value(U);

%%
% Plot the aircraft state trajectories
figure;
hold on;
colors = {'r', 'g', 'b', 'm'};  % Color for each aircraft

for i = 1:N
    plot(0:Tfinal, X_opt(i, :), 'LineWidth', 1.5, 'Color', colors{i});
end

xlabel('Time step');
ylabel('Position');
legend('Aircraft 1', 'Aircraft 2', 'Aircraft 3', 'Aircraft 4');
title('Aircraft State Trajectories');
grid on;
hold off;
