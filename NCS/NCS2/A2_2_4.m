%% To zero
clc
clear 

a = 5;
b = 3;
c = 1;
A = [a-b, 0.5-c; 0, 1];
B = [0; 1];
p1 = [-1+2j, -1-2j];
K = place(A, B, p1);

h = linspace(0, 5, 50);

x = sdpvar(4,4); % x = inv(P)
y = sdpvar(1,4);

con1 =[];
con2 =[];
con3 = x >= eye(4)^(-10);
for i = 1:length(h)
    hl1 = h(i);
    hl2 = 2 * h(i);
    
    % Case hl = h
    F1 = expm(A * hl1);
    G1 = (expm(A * hl1) - eye(size(A))) /A * B * K;
    l1 = eye(size(A));
    Fz1 = [F1-l1 * G1, zeros(size(A)); zeros(size(A)), zeros(size(A))];

    % Check stability using LMIs for hl = h
    con1 = [con1; [x, Fz1 * x; (Fz1 * x)', x] >= eye(8)^(-10)];
    
    % Case hl = 2h
    F2 = expm(A * hl2);
    G2 = (expm(A * hl2) - eye(size(A))) /A * B * K;
    l2 = expm(A * (hl2 - h(i)));
    Fz2 = [F2-l2 * G2, zeros(size(A)); zeros(size(A)), zeros(size(A))];

    % Check stability using LMIs for hl = 2h
    con2 = [con2; [x, Fz2 * x; (Fz2 * x)', x] >= eye(8)^(-10)];
end
cons = [con1, con2, con3]; 
options = sdpsettings('solver', 'sedumi');
sol = optimize(cons, [], options);
disp(sol.info);
