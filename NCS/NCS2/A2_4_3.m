
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

b = 1;
h = 0.2;
tau_max = h-0.01;
tau_min = 0;

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

F13 = F(1,3);
F23 = F(2,3);
G11 = G(1);
G21 = G(2);

% Evaluate a1 at specific values of tau
tau_vals = linspace(tau_min, tau_max, 100);  % Values of tau
a1_vals = double(subs(a1, tau, tau_vals));   % Evaluate a1 at tau_vals
max_a1 = max(a1_vals);  % Maximum value of a1
min_a1 = min(a1_vals);  % Minimum value of a1

a2_vals = double(subs(a2, tau, tau_vals));   % Evaluate a1 at tau_vals
max_a2 = max(a2_vals);  % Maximum value of a1
min_a2 = min(a2_vals);  % Minimum value of a1

% Display the maximum and minimum values of a1,a2
fprintf('Maximum value of a1: %f\n', max_a1);
fprintf('Minimum value of a1: %f\n', min_a1);
fprintf('Maximum value of a2: %f\n', max_a2);
fprintf('Minimum value of a2: %f\n', min_a2);

%% F Plot
% Convert F13 and F23 to numeric arrays
F13_vals = double(subs(F13, tau, tau_vals));
F23_vals = double(subs(F23, tau, tau_vals));

% Plot the relationship between F13 and F23
figure;
plot(F13_vals, F23_vals);
xlabel('F13');
ylabel('F23');
title('Plot of F13 against F23');

% F Vertices 
h_val = 0.2;

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


hold on 
for i=1:1:length(F_over_approx)
scatter(F_over_approx{i}(1,3),F_over_approx{i}(2,3),'green','filled');hold on
verticesF(i,:)=[F_over_approx{i}(1,3),F_over_approx{i}(2,3)];
end
p=patch(verticesF(:,1),verticesF(:,2),'cyan','FaceAlpha',0.1);
uistack(p, 'bottom');
hold off

%% G Plot 
% Convert G13 and G23 to numeric arrays
G11_vals = double(subs(G11, tau, tau_vals));
G21_vals = double(subs(G21, tau, tau_vals));

% Plot the relationship between G11 and G21
figure;
plot(G11_vals, G21_vals);
xlabel('G11');
ylabel('G21');
title('Plot of G11 against G21');

% G Vertices 
h_val = 0.2;

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
for i=1:1:length(G_over_approx)
scatter(G_over_approx{i}(1),G_over_approx{i}(2),'g','filled');hold on
verticesG(i,:)=[G_over_approx{i}(1),G_over_approx{i}(2)];
end
p=patch(verticesG(:,1),verticesG(:,2),'yellow','FaceAlpha',0.1);
uistack(p, 'bottom');



