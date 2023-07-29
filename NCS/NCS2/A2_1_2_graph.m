% Question 1-2
%% To zero (is better than to hold)
clc
clear 

a = 5;
b = 3;
c = 1;
A = [a-b, 0.5-c; 0, 1];
B = [0; 1];
p1 = [-1+2j, -1-2j];
K = place(A, B, p1);

numSamples = 100; % number of sampling times to test
maxTime = 0.45; % maximum sampling time to test
hVals = linspace(0.001, maxTime, numSamples); % generate sampling times to test
maxEigVals = zeros(numSamples, 1); % initialize matrix to store maximum eigenvalues

figure;
hold on;

for i = 1:numSamples
    h = hVals(i);
    F = expm(A*h);
    G = (expm(A*h) - eye(size(A))) / A * B;
    Fz1 = [F, zeros(size(A));  eye(size(A)), zeros(size(A))];
    Fz0 = [F-G*K, zeros(size(A));  eye(size(A)), zeros(size(A))];
    Fzf = Fz1 * Fz0 *Fz0; % x_{k+1} = Fhf * x_{k}

    eigh = eig(Fzf);
    maxEigVals(i) = max(abs(eigh));

end

% Find the stable range
stableIdx = find(maxEigVals < 1);
unstableIdx = find(maxEigVals >= 1);

% Plot stable range
scatter(hVals(stableIdx), maxEigVals(stableIdx), 20, 'g', 'filled');

% Plot unstable range
scatter(hVals(unstableIdx), maxEigVals(unstableIdx), 20, 'r', 'filled');

% Customize plot
xlabel('Sampling Time (h)');
ylabel('Maximum Eigenvalue');
title('Stability Analysis: To-Zero Approach');
legend('Stable Range', 'Unstable Range');

% Display results
disp("To-Zero Approach:");
disp("Stable Range:");
disp(hVals(stableIdx)');
disp("Best Stable Range:");
disp(max(hVals(stableIdx)));

%% To hold
clc
clear 

a = 5;
b = 3;
c = 1;
A = [a-b, 0.5-c; 0, 1];
B = [0; 1];
p1 = [-1+2j, -1-2j];
K = place(A, B, p1);

numSamples = 100; % number of sampling times to test
maxTime = 0.45; % maximum sampling time to test
hVals = linspace(0.001, maxTime, numSamples); % generate sampling times to test
maxEigVals = zeros(numSamples, 1); % initialize matrix to store maximum eigenvalues

figure;
hold on;

% To hold 
for i = 1:numSamples
    h = hVals(i);
    F = expm(A*h);
    G = (expm(A*h) - eye(size(A))) / A * B;
    Fh1 = [F, -G*K;  eye(size(A)), zeros(size(A))];
    Fh0 = [F-G*K, zeros(size(A));  eye(size(A)), zeros(size(A))];
    Fhf = Fh1 * Fh0 *Fh0; % x_{k+1} = Fhf * x_{k}

    eigh = eig(Fhf);
    maxEigVals(i) = max(abs(eigh));


end

% Find the stable range
stableIdx = find(maxEigVals < 1);
unstableIdx = find(maxEigVals >= 1);

% Plot stable range
scatter(hVals(stableIdx), maxEigVals(stableIdx), 20, 'g', 'filled');

% Plot unstable range
scatter(hVals(unstableIdx), maxEigVals(unstableIdx), 20, 'r', 'filled');

% Customize plot
xlabel('Sampling Time (h)');
ylabel('Maximum Eigenvalue');
title('Stability Analysis: To-Hold Approach');
legend('Stable Range', 'Unstable Range');

% Display results
disp("To-Hold Approach:");
disp("Stable Range:");
disp(hVals(stableIdx)');
disp("Best Stable Range:");
disp(max(hVals(stableIdx)));

