clc
clear('yalmip')

%To zero
h=0.001;
h_max=h;

a = 5;
b = 3;
c = 1;
A = [a-b, 0.5-c; 0, 1];
B = [0; 1];
p1 = [-1+2j, -1-2j];
K = place(A, B, p1);

while 1


P1=sdpvar(2,2);
P2=sdpvar(2,2);
P3=sdpvar(2,2);
P4=sdpvar(2,2);
P5=sdpvar(2,2);
P6=sdpvar(2,2);

obj=0;
F=expm(A*h);
G=(expm(A*h)-eye(2))/A*B;

F_m0=F-G*K;
F_m1=F;

p0=0.70;%no packet loss, when non collision
p1=0.49;%no packet loss, when collision

%no packet loss
A1 = F_m0;
A3 = F_m0;
A5 = F_m0;
%packet loss
A2 = F_m1;
A4 = F_m1;
A6 = F_m1;


con7=[];
con8=[];
con9=[];
con10=[];
con11=[];
con12=[];


con1 = P1 >=eye(2)*1e-8;
con2 = P2 >=eye(2)*1e-8;
con3 = P3 >=eye(2)*1e-8;
con4 = P4 >=eye(2)*1e-8; 
con5 = P5 >=eye(2)*1e-8;
con6 = P6 >=eye(2)*1e-8;


con7 = [con7; P1-p0*A3'*P3*A3-(1-p0)*A4'*P4*A4 >= eye(2)*1e-8];
con8 = [con8; P2-p0*A3'*P3*A3-(1-p0)*A4'*P4*A4 >= eye(2)*1e-8];

con9 = [con9; P3-p1*A5'*P5*A5-(1-p1)*A6'*P6*A6 >= eye(2)*1e-8];
con10= [con10; P4-p1*A5'*P5*A5-(1-p1)*A6'*P6*A6 >= eye(2)*1e-8];

con11= [con11; P5-p0*A1'*P1*A1-(1-p0)*A2'*P2*A2 >= eye(2)*1e-8];
con12= [con12; P6-p0*A1'*P1*A1-(1-p0)*A2'*P2*A2 >= eye(2)*1e-8];


cons = [con1, con2, con3, con4, con5, con6, con7, con8, con9, con10, con11, con12];
options = sdpsettings('solver', 'sedumi');
sol = optimize(cons, [], options);
if sol.problem==0
h_max=h;
h=h+0.01;
else
disp(['The maximal possible h is ', num2str(h_max)] );
break;
end
end
%h = 0.291
% 0.031 when P0 = 0.70
disp(sol.info);


%% To hold 

clc
clear
clear('yalmip')

h=0.001;
h_max=h;

a = 5;
b = 3;
c = 1;
A = [a-b, 0.5-c; 0, 1];
B = [0; 1];
p1 = [-1+2j, -1-2j];
K = place(A, B, p1);

while 1


P1=sdpvar(4);
P2=sdpvar(4);
P3=sdpvar(4);
P4=sdpvar(4);
P5=sdpvar(4);
P6=sdpvar(4);

obj=0;
F=expm(A*h);
G=(expm(A*h)-eye(2))/A*B;


F_m1=[F, -G*K;  eye(size(A)), zeros(size(A))];
F_m0=[F-G*K, zeros(size(A));  eye(size(A)), zeros(size(A))];

p0=0.99;%no packet loss, when non collision
p1=0.49;%no packet loss, when collision

%no packet loss
A1 = F_m0;
A3 = F_m0;
A5 = F_m0;
%packet loss
A2 = F_m1;
A4 = F_m1;
A6 = F_m1;


con7=[];
con8=[];
con9=[];
con10=[];
con11=[];
con12=[];


con1 = P1 >=eye(4)*1e-8;
con2 = P2 >=eye(4)*1e-8;
con3 = P3 >=eye(4)*1e-8;
con4 = P4 >=eye(4)*1e-8; 
con5 = P5 >=eye(4)*1e-8;
con6 = P6 >=eye(4)*1e-8;


con7 = [con7; P1-p0*A3'*P3*A3-(1-p0)*A4'*P4*A4 >= eye(4)*1e-8];
con8 = [con8; P2-p0*A3'*P3*A3-(1-p0)*A4'*P4*A4 >= eye(4)*1e-8];

con9 = [con9; P3-p1*A5'*P5*A5-(1-p1)*A6'*P6*A6 >= eye(4)*1e-8];
con10= [con10; P4-p1*A5'*P5*A5-(1-p1)*A6'*P6*A6 >= eye(4)*1e-8];

con11= [con11; P5-p0*A1'*P1*A1-(1-p0)*A2'*P2*A2 >= eye(4)*1e-8];
con12= [con12; P6-p0*A1'*P1*A1-(1-p0)*A2'*P2*A2 >= eye(4)*1e-8];


cons = [con1, con2, con3, con4, con5, con6, con7, con8, con9, con10, con11, con12];
options = sdpsettings('solver', 'sedumi');
sol = optimize(cons, [], options);
if sol.problem==0
h_max=h;
h=h+0.01;
else
disp(['The maximal possible h is ', num2str(h_max)] );
break;
end
end
%0.251
% 0.171 when P0=0.7
disp(sol.info);