clc
clear
clear('yalmip')

a = 5;
b = 3;
c = 1;
A = [a-b, 0.5-c; 0, 1];
B = [0; 1];
C = eye(2);
D = zeros(2,1);

p1 = [-1+2j, -1-2j];
K = place(A, B, p1);
k1 = [-26, 5];

% Declare the optimization variables
P = sdpvar(2, 2, 'symmetric');
Q = eye(2);

% Define the LMI constraints
cons = [P >= eye(2)*1e-3, Q >= eye(2)*1e-3,... 
        Q >= -((A-B*k1)'*P+P*(A-B*k1)),...
        Q <= -((A-B*k1)'*P+P*(A-B*k1))];

options = sdpsettings('solver', 'sedumi');
sol = optimize(cons, [], options);
disp(sol.info);

P=double(P);

%% Q5-3
havg =  3*10/(23+ 23+ 23) %largest average inter sample time
sys = ss(A,B,eye(2),0);
sysd = c2d(sys, havg);
max(abs(eig(sysd.A-sysd.B*k1)))
% Q5-4
havgop = havg/2