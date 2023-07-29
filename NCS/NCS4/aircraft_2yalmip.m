% Add path to YALMIP and a solver (e.g., SDPT3)
addpath('/path/to/yalmip')
addpath('/path/to/solver')

% Define variables and parameters
Tfinal = 5;   % Horizon length
umax = 60;    % Control limit

% Define system matrices and initial states for each agent
A1 = [1 0 2 0; 0 1 0 2; 0 0 3 0; 0 0 0 3];
B1 = [2 0; 0 2; 3 0; 0 3];
x01 = [-10; 10; -1; 1];

A2 = [1 0 3 0; 0 1 0 3; 0 0 7 0; 0 0 0 7];
B2 = [3 0; 0 3; 7 0; 0 7];
x02 = [10; 10; 1; 1];

A3 = [1 0 1 0; 0 1 0 1; 0 0 1.1 0; 0 0 0 1.1];
B3 = [1 0; 0 1; 1.1 0; 0 1.1];
x03 = [10; -10; 1; -1];

A4 = [1 0 6 0; 0 1 0 6; 0 0 20 0; 0 0 0 20];
B4 = [6 0; 0 6; 20 0; 0 20];
x04 = [-10; -10; -1; -1];

% Create optimization variables
x1 = sdpvar(4, Tfinal, 'full');
u1 = sdpvar(2, Tfinal-1, 'full');
x2 = sdpvar(4, Tfinal, 'full');
u2 = sdpvar(2, Tfinal-1, 'full');
x3 = sdpvar(4, Tfinal, 'full');
u3 = sdpvar(2, Tfinal-1, 'full');
x4 = sdpvar(4, Tfinal, 'full');
u4 = sdpvar(2, Tfinal-1, 'full');

% Create objective function
obj = 0;
for t = 1:Tfinal-1
    obj = obj + x1(:,t)' * x1(:,t) + u1(:,t)' * u1(:,t) ...
                + x2(:,t)' * x2(:,t) + u2(:,t)' * u2(:,t) ...
                + x3(:,t)' * x3(:,t) + u3(:,t)' * u3(:,t) ...
                + x4(:,t)' * x4(:,t) + u4(:,t)' * u4(:,t);
end

% Create constraints
constr = [];
for t = 1:Tfinal-1
    constr = [constr, x1(:,t+1) == A1 * x1(:,t) + B1 * u1(:,t)];
    constr = [constr, x2(:,t+1) == A2 * x2(:,t) + B2 * u2(:,t)];
    constr = [constr, x3(:,t+1) == A3 * x3(:,t) + B3 * u3(:,t)];
    constr = [constr, x4(:,t+1) == A4 * x4(:,t) + B4 * u4(:,t)];
    
    constr = [constr, abs(u1(:,t)) <= umax / Tfinal];
    constr = [constr, abs(u2(:,t)) <= umax / Tfinal];
    constr = [constr, abs(u3(:,t)) <= umax / Tfinal];
    constr = [constr, abs(u4(:,t)) <= umax / Tfinal];
end

% Set initial states
constr = [constr, x1(:,1) == x01];
constr = [constr, x2(:,1) == x02];
constr = [constr, x3(:,1) == x03];
constr = [constr, x4(:,1) == x04];

% Initialize Lagrange multipliers
theta1 = zeros(2, Tfinal-1);
theta2 = zeros(2, Tfinal-1);
theta3 = zeros(2, Tfinal-1);
theta4 = zeros(2, Tfinal-1);
alpha = 0.01;   % Step size

% Run projected subgradient method for dual decomposition
for iter = 1:100
    % Solve the optimization problem with fixed Lagrange multipliers
    ops = sdpsettings('solver', 'sdpt3');
    optimize(constr + [abs(B1' * x1(:,2:Tfinal)) <= umax / Tfinal + theta1, ...
                       abs(B2' * x2(:,2:Tfinal)) <= umax / Tfinal + theta2, ...
                       abs(B3' * x3(:,2:Tfinal)) <= umax / Tfinal + theta3, ...
                       abs(B4' * x4(:,2:Tfinal)) <= umax / Tfinal + theta4], obj, ops);
    
    % Retrieve the optimal values
    x1_opt = value(x1);
    u1_opt = value(u1);
    x2_opt = value(x2);
    u2_opt = value(u2);
    x3_opt = value(x3);
    u3_opt = value(u3);
    x4_opt = value(x4);
    u4_opt = value(u4);
    
    % Update Lagrange multipliers for control input constraints
    theta1 = max(0, theta1 - alpha * (abs(B1' * x1_opt(:,2:Tfinal)) - umax / Tfinal));
    theta2 = max(0, theta2 - alpha * (abs(B2' * x2_opt(:,2:Tfinal)) - umax / Tfinal));
    theta3 = max(0, theta3 - alpha * (abs(B3' * x3_opt(:,2:Tfinal)) - umax / Tfinal));
    theta4 = max(0, theta4 - alpha * (abs(B4' * x4_opt(:,2:Tfinal)) - umax / Tfinal));
end

% Print final states
disp('Final states:');
disp(['Agent 1: ', mat2str(x1_opt(:,Tfinal))]);
disp(['Agent 2: ', mat2str(x2_opt(:,Tfinal))]);
disp(['Agent 3: ', mat2str(x3_opt(:,Tfinal))]);
disp(['Agent 4: ', mat2str(x4_opt(:,Tfinal))]);
