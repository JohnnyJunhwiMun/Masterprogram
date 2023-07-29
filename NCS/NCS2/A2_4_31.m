clc
clear

syms h tau
a = 5;
bl = 3;
c = 1;
A = [a-bl, 0.5-c; 0, 1];
B = [0; 1];
p1 = [-1+2j, -1-2j];
K = place(A, B, p1);


Fx = expm(A*h);
Fu = (expm(A*h) - expm(A*(h-tau)))/A * B;
Gu = (expm(A*(h-tau)) - eye(size(A)))/A * B;
F = [Fx, Fu; zeros(1,2) zeros(1)];
Fs = simplify(F);
G = [Gu; 1];
Gs = simplify(G);


 a1=exp(2*(h-tau));
 a2=exp(h-tau);

F1 = [0,0,1/4;0,0,0;0,0,0];
F2 = [0,0,-1/2; 0,0,-1; 0,0,0];
F0 = Fs-a1*F1 -a2*F2;

G1 = [-1/4; 0; 0];
G2 = [1/2; 1; 0];
G0 = [-1/4; -1; 1];
Gsk = G0 +a1*G1+ a2*G2;


% Define the grid of h and τ values
h_min = 0;
h_max = 0.3;
num_h = 30;  % Increase the number of h values for a smoother plot

h_vals = linspace(h_min, h_max, num_h);  % Values of h

% Initialize a matrix to store stability results
stability_results = zeros(1, num_h);

% Iterate over the values of h
for i = 1:num_h
    h_val = h_vals(i);

    a1=exp(2*(h_val-tau));
    a2=exp(h_val-tau);

    tau_max = h_val-0.01;
    tau_min = 0;

    % Evaluate a1 at specific values of tau
    tau_vals = linspace(tau_min, tau_max, 100);  % Values of tau
    a1_vals = double(subs(a1, tau, tau_vals));   % Evaluate a1 at tau_vals
    max_a1 = max(a1_vals);  % Maximum value of a1
    min_a1 = min(a1_vals);  % Minimum value of a1

    a2_vals = double(subs(a2, tau, tau_vals));   % Evaluate a1 at tau_vals
    max_a2 = max(a2_vals);  % Maximum value of a2
    min_a2 = min(a2_vals);  % Minimum value of a2

%F
    F_over_approx{1}=double(subs(F0+max_a1*F1+max_a2*F2,h,h_val));
    Fk=double(subs(F0+min_a1*F1+max_a2*F2,h,h_val));
    F_over_approx{2}=double(subs(F0+min_a1*F1+min_a2*F2,h,h_val));

    % Calculate the additional vertices on the lines
    Fk3 = F_over_approx{2} + (F_over_approx{2} + Fk) / 40;
    Fk4 = F_over_approx{1} + (Fk - F_over_approx{1}) / 40;

    F_over_approx{3} =Fk3 + (F_over_approx{2}- Fk3)/2;
    F_over_approx{6} =F_over_approx{1} + (-F_over_approx{1}+ Fk4)/2;

    F_over_approx{5} =Fk4 + (Fk3- Fk4)/3;
    F_over_approx{4} =Fk3 + (-Fk3+ Fk4)/3;

%G
    G_over_approx{1}=double(subs(G0+max_a1*G1+max_a2*G2,h,h_val));
    Gk=double(subs(G0+min_a1*G1+max_a2*G2,h,h_val));
    G_over_approx{2}=double(subs(G0+min_a1*G1+min_a2*G2,h,h_val));

    % Calculate the additional vertices on the lines
    Gk3 = G_over_approx{2} + (G_over_approx{2} + Gk) / 40;
    Gk4 = G_over_approx{1} + (Gk - G_over_approx{1}) / 40;

    G_over_approx{3} =Gk3 + (G_over_approx{2}- Gk3)/2;
    G_over_approx{6} =G_over_approx{1} + (-G_over_approx{1}+ Gk4)/2;

    G_over_approx{5} =Gk4 + (Gk3- Gk4)/3;
    G_over_approx{4} =Gk3 + (-Gk3+ Gk4)/3;

    hold on
    verticesFF = zeros(length(F_over_approx), 2);
    for j = 1:length(F_over_approx)
        verticesFF(j, :) = [F_over_approx{j}(1, 3), F_over_approx{j}(2, 3)];
    end
    
    hold on
    verticesGG = zeros(length(G_over_approx), 2);
    for j = 1:length(G_over_approx)
        verticesGG(j, :) = [G_over_approx{j}(1), G_over_approx{j}(2)];
    end
    
    con = [];
    
    % Check stability using LMIs
    P = sdpvar(2, 2);
    con1 = P >= eye(2) * 1e-8;
    for q = 1:6
        % GES
        con = [con; (verticesFF(q, :)' - verticesGG(q, :)'*K)'*P*(verticesFF(q, :)' - verticesGG(q, :)'*K) - P <= eye(2) * 1e-8];
    end
    cons = [con, con1];
    options = sdpsettings('solver', 'sedumi');
    sol = optimize(cons, [], options);
    
    % Update stability results
    if sol.problem == 0
        stability_results(i) = 1;  % System is stable
    else
        stability_results(i) = 0;  % System is unstable
    end
end
disp(sol.info);

% Plot the stable and unstable regions
figure()
plot(h_vals(stability_results == 1), zeros(1, sum(stability_results)), 'go', 'MarkerFaceColor', 'green');
hold on %0.197
plot(h_vals(stability_results == 0), zeros(1, sum(~stability_results)), 'ro', 'MarkerFaceColor', 'red');
xlabel('h');
title('Stability Analysis of h');
legend('Stable', 'Unstable');