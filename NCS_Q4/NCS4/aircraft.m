% state vector is (x,y,xdot,ydot)

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

Tfinal=5;   %horizon length
umax=60;    %control limit
%% 
% Define variables and parameters
Tfinal = 5;   % Horizon length
umax = 60;    % Control limit
epsilon = 1e-6; % Convergence threshold

% Define system matrices and initial conditions
A = {A1, A2, A3, A4};
B = {B1, B2, B3, B4};
x0 = {x01, x02, x03, x04};

% Create CVX variables
cvx_begin
    variable u(2, Tfinal, 4) % u(i,t,k) represents the control input of aircraft k at time t
    variable x(4, Tfinal+1, 4) % x(:,t,k) represents the state of aircraft k at time t
    variable xf(4) % Common terminal state
    
    % Define objective function
    obj = 0;
    for k = 1:4
        obj = obj + sum(sum(square(x(:,1:Tfinal,k)))) + sum(sum(square(u(:,:,k))));
    end
    
    % Define constraints
    constr = [];
    for k = 1:4
        % Initial states
        constr = [constr, x(:,1,k) == x0{k}];
        
        % Dynamics constraints for each time step
        for t = 1:Tfinal
            constr = [constr, x(:,t+1,k) == A{k}*x(:,t,k) + B{k}*u(:,t,k)];
            
            % Control input constraints
            constr = [constr, norm(u(:,t,k), 2) <= umax/Tfinal];
        end
        
        % Terminal state constraint
        constr = [constr, x(:,Tfinal+1,k) == xf];
    end
    
    % Dual decomposition with Lagrangian function
    lambda = sdpvar(4, Tfinal);
    lagrangian = 0;
    for k = 1:4
        for t = 1:Tfinal
            lagrangian = lagrangian...
            + (x(:,t+1,k) - A{k}*x(:,t,k) - B{k}*u(:,t,k) - xf)'*lambda(:,t);
        end
    end
    obj = obj + lagrangian;
    
    % Set CVX options
    cvx_solver('sedumi');
    
    minimize(obj)
    subject to
        constr
cvx_end

% Initialize variables
x_hist = zeros(4, Tfinal+1, 4);
error_seq = zeros(1, Tfinal+1);

% Store the optimal values and compute the error
optimal_u = u;
optimal_x = x;
optimal_xf = xf;

% Compute the error sequence
for iteration = 0:Tfinal
    x_hist(:, :, iteration+1) = optimal_x(:, :, 1);
    error_seq(iteration+1) = norm(optimal_x(:, :, 1) - optimal_xf) / norm(optimal_xf);
end

% Plot the convergence
figure;
plot(0:Tfinal, error_seq);
xlabel('Subgradient Iteration');
ylabel('Error');
title('Convergence of Dual Decomposition');
