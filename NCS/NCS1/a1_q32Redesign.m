% Define the system matrices A and B based on the given system
a = 5;
b = 3;
c = 1;
A = [a-b, 0.5-c; 0, 1];
B = [0; 1];

% Define the selected sampling interval(tau = h)
h = 0.05;
tau = linspace(0, 1.5*h, 15);


x = sdpvar(4,4); % x = inv(P)
y = sdpvar(1,4); % y = K*inv(P), dim(K) = dim(y)

con1 =[];
con2 = [x >= eye(4,4)^(-10)];
for i =1:length(tau)
%when 0=<tau<h
if tau(i) < h
        %new augumented system to match dimension
        Fx = expm(A*h);
        Fu = (expm(A*h) - expm(A*(h-tau(i))))/A * B;
        G1 = (expm(A*(h-tau(i))) - eye(size(A)))/A * B;
        Fe = [Fx, Fu, zeros(2,1); zeros(1,4); zeros(1,2), 1, 0];
        Ge = [G1; 1; 0];

%when h<tau<=1.5h
else 
        Fxh = expm(A*h);
        Fuh = (expm(A*h) - expm(A*(2*h-tau(i))))/A * B;
        Gh = (expm(A*(2*h-tau(i))) - eye(size(A)))/A * B;
        Fe = [Fxh, Fuh, Gh; zeros(size(Fxh)),[0;1],zeros(2,1)]; 
        Ge = [zeros(size(Fxh(:,1))); 1; 0];
end
% Check stability using LMIs(for discrete-time case)
con1 = [con1; [x, Fe*x-Ge*y; x*(Fe)'-(Ge*y)', x] >= eye(8,8)^(-10)];
end
cons = [con1, con2];
options = sdpsettings('solver', 'sedumi');
sol = optimize(cons, [], options);
disp(sol.info);

if sol.problem == 0
    disp('The sampled-data controller is stable.');
    X = value(x)
    Y = value(y)
    K = Y * inv(X)
else
    disp('The sampled-data controller is unstable.');
end
Ke1 = Fe - Ge*K;
eig(Ke1); %stable

% Plot
result = zeros(size(tau));
for i = 1:length(tau)
    % When 0 <= tau < h
    if tau(i) < h
        Fx = expm(A*h);
        Fu = (expm(A*h) - expm(A*(h-tau(i))))/A * B;
        G1 = (expm(A*(h-tau(i))) - eye(size(A)))/A * B;
        Fe = [Fx, Fu, zeros(2,1); zeros(1,4); zeros(1,2), 1, 0];
        Ge = [G1; 1; 0];
    % When h < tau <= 1.5h
    elseif tau(i) <= 1.5*h
        Fxh = expm(A*h);
        Fuh = (expm(A*h) - expm(A*(2*h-tau(i))))/A * B;
        Gh = (expm(A*(2*h-tau(i))) - eye(size(A)))/A * B;
        Fe = [Fxh, Fuh, Gh; zeros(size(Fxh)),[0;1],zeros(2,1)]; 
        Ge = [zeros(size(Fxh(:,1))); 1; 0];
    end
    
    Ke1 = Fe - Ge*K;
    result(i) = max(abs(eig(Ke1))) ;
    %result(i) = max(abs(eig(Ke1))) < 1; %it should be all '1'
end

% Plotting
plot(tau, result, 'o-');
xlabel('tau');
ylabel('Stability Measure');
legend('maximum absolute eigenvalue magnitude of the controller');