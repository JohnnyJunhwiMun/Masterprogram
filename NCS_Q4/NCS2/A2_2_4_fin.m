%% to zero
clc
clear
clear('yalmip')

a = 5;
b = 3;
c = 1;
A = [a-b, 0.5-c; 0, 1];
B = [0; 1];
p1 = [-1+2j, -1-2j];
K = place(A, B, p1);

h=0.001;
h_max=h;
while 1
X=sdpvar(4);
obj=0;
F=expm(A*h);
G=(expm(A*h)-eye(2))/A*B;

%To zero

F_m1 = [F, zeros(size(A));  eye(size(A)), zeros(size(A))];
F_m0 = [F-G*K, zeros(size(A));  eye(size(A)), zeros(size(A))];
F_cl_1=F_m1*F_m0*F_m0;
F_cl_2=F_m0*F_m0*F_m0;

cons =X>=eye(4)*1e-8;
cons=[cons;[X F_cl_1*X;X*F_cl_1' X]>=eye(8)*1e-3];
cons=[cons;[X F_cl_2*X;X*F_cl_2' X]>=eye(8)*1e-3];

result=optimize(cons,obj);
x=value(X);
P=eye(4)/x;
if result.problem==0
h_max=h;
h=h+0.01;
else
disp(['The maximal possible h is ', num2str(h_max)] );
break;
end
end
%0.271
disp(result.info);

%% To hold 
clc
clear
clear('yalmip')

a = 5;
b = 3;
c = 1;
A = [a-b, 0.5-c; 0, 1];
B = [0; 1];
p1 = [-1+2j, -1-2j];
K = place(A, B, p1);

h=0.001;
h_max=h;
while 1
X=sdpvar(4);
obj=0;
F=expm(A*h);
G=(expm(A*h)-eye(2))/A*B;

%To hold
F_m1=[F, -G*K;  eye(size(A)), zeros(size(A))];
F_m0=[F-G*K, zeros(size(A));  eye(size(A)), zeros(size(A))];

F_cl_1=F_m1*F_m0*F_m0;
F_cl_2=F_m0*F_m0*F_m0;

cons =X>=eye(4)*1e-8;
cons=[cons;[X F_cl_1*X;X*F_cl_1' X]>=eye(8)*1e-3];
cons=[cons;[X F_cl_2*X;X*F_cl_2' X]>=eye(8)*1e-3];

result=optimize(cons,obj);
x=value(X);
P=eye(4)/x;
if result.problem==0
h_max=h;
h=h+0.01;
else
disp(['The maximal possible h is ', num2str(h_max)] );
break;
end
end
%0.231
disp(result.info);
