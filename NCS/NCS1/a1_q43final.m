% Define the system matrices A and B based on the given system
a = 5;
b = 3;
c = 1;
A = [a-b, 0.5-c; 0, 1];
B = [0; 1];

% Initialize variables to store the maximum h retaining stability
max_h_stable = 0;
max_h_range = 0;

% Initialize arrays to store eigenvalues and h values
eigenvalues = [];
h_values = [];

for h = linspace(0, 5, 100)
    % Define the time-varying inter-sample times hk and delays tauk
    hk = 0.5 * h;     % Inter-sample time
    tau = 0.5 * hk;   % Time-varying delays

    % Compute Fhk and Ghk using matrices_q4 function
    [Fhk, Ghk] = matrices_q4_ag(A, B, h, tau);

    x = sdpvar(4, 4); % x = inv(P)
    y = sdpvar(1, 4); % y = K*inv(P), dim(K) = dim(y)

    % Check stability using LMIs (for discrete-time case)
    con1 = [[x, Fhk*x-Ghk*y; x*(Fhk)'-(Ghk*y)', x] >= 0];
    con2 = [x >= 0];
    cons = [con1, con2];
    options = sdpsettings('solver', 'sedumi');
    sol = optimize(cons, [], options);

    if sol.problem == 0
        disp('The sampled-data controller is stable.');
        X = value(x);
        Y = value(y);
        K = Y * inv(X);
        % Update the maximum h retaining stability
        max_h_stable = h;
        max_h_range = max_h_range + 1;
    else
        disp('The sampled-data controller is unstable.');
        break; % Exit the loop if the controller becomes unstable
    end
    disp(sol.info);

    % Compute the absolute eigenvalues and store them
    Ke1 = Fhk - Ghk*K;
    eigvals = max(abs(eig(Ke1)));
    eigenvalues = [eigenvalues; eigvals.'];
    h_values = [h_values; h];
end

% Display the maximum h retaining stability
disp(['Maximum h retaining stability: ' num2str(max_h_stable)]);
disp(['Number of stable h values: ' num2str(max_h_range)]);


% Plot the absolute eigenvalues with respect to h
figure;
plot(h_values, eigenvalues, 'o-');
xlabel('h');
ylabel('Absolute Eigenvalue Magnitude');
legend('maximum absolute eigenvalue magnitude of the controller');
title('Absolute Eigenvalue Magnitudes of the Controller');
