%% Changing Q

yalmip('clear')
clear all

% Model data
Ra= 1;                              %Armature electric resistance[Ohm]
La= 0.5;                            %Electric inductance[H]
Ki= 0.023;                          %Motor torque constant= Electromotive force constant[Nm/A]
Jm= 0.01;                           %Moment of inertia of the rotor[kgm^2]
Kb= 0.023;
Bm= 0.00003;

A=[-Ra/La -Kb/La; Ki/Jm -Bm/Jm];
B=[1/La; 0];
C=[0 1];
D=0;
sys=ss(A,B,C,D);           %Denominator and Nominator of Transfer function
Ts=0.2;
sys_d=c2d(sys,Ts);              %Transfer function

nx = 2; % Number of states
nu = 1; % Number of inputs

%% Q0
Q = diag([100,10]); % State weighting matrix
R = 1; % Input weighting matrix
N=7;
K=10.^10;%soften the constraint
[P,gain,L] = idare(sys_d.A,sys_d.B,Q,R); gain = -gain; %for terminal cost(V_f)

% Control invariant set
system = LTISystem('A', sys_d.A, 'B', sys_d.B);
system.x.min = [-3; -12];
system.x.max = [3; 12];
system.u.min = -5;
system.u.max = 5;
system.x.penalty = QuadFunction(Q);
system.u.penalty = QuadFunction(R);

system.x.with('terminalSet');
system.x.terminalSet = system.LQRSet();
system.x.with('terminalPenalty');
system.x.terminalPenalty = system.LQRPenalty();

InvSet = system.LQRSet();

% Set Optimizer
u_opt = sdpvar(repmat(nu,1,N),repmat(1,1,N));
x_opt = sdpvar(repmat(nx,1,N+1),repmat(1,1,N+1));
z_opt=sdpvar(1,1); %for constraints
q_opt=sdpvar(1,1); %for constraints
q1_opt=sdpvar(1,1); %for constraints
constraints = [];
objective = 0;

% Simulation
for k = 1:N
    % x{k+1} = sys_d.A*x{k} + sys_d.B*u{k};
    objective = objective + transpose(x_opt{k})*Q*(x_opt{k}) + transpose(u_opt{k})*R*(u_opt{k}) +K*z_opt+ K*q_opt + K*q1_opt;
    constraints = [constraints, x_opt{k+1} == sys_d.A*x_opt{k} + sys_d.B*u_opt{k}];
    constraints = [constraints, u_opt{k}<= 5+z_opt, -5-z_opt<= u_opt{k}, z_opt>=0];           %voltage
    constraints = [constraints, x_opt{k}(1)<= 3+q_opt, -3-q_opt<= x_opt{k}(1), q_opt>=0];     %current
    constraints = [constraints, x_opt{k}(2)<= 12+q1_opt, -12-q1_opt<= x_opt{k}(2), q1_opt>=0];%angular velocity
end
constraints=[constraints,InvSet.A*x_opt{N+1}<=InvSet.b]; %Terminal constraint
objective = objective + transpose(x_opt{N+1})*P*(x_opt{N+1}); %adding the terminal cost at the end
controller = optimizer(constraints, objective,[],x_opt{1},[u_opt{:}]);

x = [6;15];
hold on
implementedU = [];
X=[];
for i = 1:50
    X=[X,x];
    U = controller{x};
    stairs(i:i+length(U)-1,U,'r')
    x = sys_d.A*x + sys_d.B*U(1);
    implementedU = [implementedU;U(1)];
end
%% Q1
Q1 = diag([100,50]); % State weighting matrix
[P1,gain1,L1] = idare(sys_d.A,sys_d.B,Q1,R); gain1 = -gain1; %for terminal cost(V_f)

% Control invariant set
system = LTISystem('A', sys_d.A, 'B', sys_d.B);
system.x.min = [-3; -12];
system.x.max = [3; 12];
system.u.min = -5;
system.u.max = 5;
system.x.penalty = QuadFunction(Q1);
system.u.penalty = QuadFunction(R);

system.x.with('terminalSet');
system.x.terminalSet = system.LQRSet();
system.x.with('terminalPenalty');
system.x.terminalPenalty = system.LQRPenalty();

InvSet = system.LQRSet();

% Set Optimizer
u_opt = sdpvar(repmat(nu,1,N),repmat(1,1,N));
x_opt= sdpvar(repmat(nx,1,N+1),repmat(1,1,N+1));
z_opt=sdpvar(1,1); %for constraints
q_opt=sdpvar(1,1); %for constraints
q1_opt=sdpvar(1,1); %for constraints
constraints = [];
objective = 0;

% Simulation
for k = 1:N
    objective = objective + transpose(x_opt{k})*Q1*(x_opt{k}) + transpose(u_opt{k})*R*(u_opt{k}) +K*z_opt+ K*q_opt + K*q1_opt;
    constraints = [constraints, x_opt{k+1} == sys_d.A*x_opt{k} + sys_d.B*u_opt{k}];
    constraints = [constraints, u_opt{k}<= 5+z_opt, -5-z_opt<= u_opt{k}, z_opt>=0];           %voltage
    constraints = [constraints, x_opt{k}(1)<= 3+q_opt, -3-q_opt<= x_opt{k}(1), q_opt>=0];     %current
    constraints = [constraints, x_opt{k}(2)<= 12+q1_opt, -12-q1_opt<= x_opt{k}(2), q1_opt>=0];%angular velocity
end
constraints=[constraints,InvSet.A*x_opt{N+1}<=InvSet.b]; %Terminal constraint
objective = objective + transpose(x_opt{N+1})*P1*(x_opt{N+1}); %adding the terminal cost at the end
controller1 = optimizer(constraints, objective,[],x_opt{1},[u_opt{:}]);

x = [6;15];
hold on
implementedU1 = [];
X1=[];
for i = 1:50
    X1=[X1,x];
    U1 = controller1{x};
    stairs(i:i+length(U1)-1,U1,'r')
    x = sys_d.A*x + sys_d.B*U1(1);
    implementedU1 = [implementedU1;U1(1)];
end

%% Q2
Q2 = diag([100,100]); % State weighting matrix
[P2,gain2,L2] = idare(sys_d.A,sys_d.B,Q2,R); gain2 = -gain2; %for terminal cost(V_f)

% Control invariant set
system = LTISystem('A', sys_d.A, 'B', sys_d.B);
system.x.min = [-3; -12];
system.x.max = [3; 12];
system.u.min = -5;
system.u.max = 5;
system.x.penalty = QuadFunction(Q2);
system.u.penalty = QuadFunction(R);

system.x.with('terminalSet');
system.x.terminalSet = system.LQRSet();
system.x.with('terminalPenalty');
system.x.terminalPenalty = system.LQRPenalty();

InvSet = system.LQRSet();

% Set Optimizer
u_opt = sdpvar(repmat(nu,1,N),repmat(1,1,N));
x_opt= sdpvar(repmat(nx,1,N+1),repmat(1,1,N+1));
z_opt=sdpvar(1,1); %for constraints
q_opt=sdpvar(1,1); %for constraints
q1_opt=sdpvar(1,1); %for constraints
constraints = [];
objective = 0;

% Simulation
for k = 1:N
    % x{k+1} = sys_d.A*x{k} + sys_d.B*u{k};
    objective = objective + transpose(x_opt{k})*Q2*(x_opt{k}) + transpose(u_opt{k})*R*(u_opt{k}) +K*z_opt+ K*q_opt + K*q1_opt;
    constraints = [constraints, x_opt{k+1} == sys_d.A*x_opt{k} + sys_d.B*u_opt{k}];
    constraints = [constraints, u_opt{k}<= 5+z_opt, -5-z_opt<= u_opt{k}, z_opt>=0];           %voltage
    constraints = [constraints, x_opt{k}(1)<= 3+q_opt, -3-q_opt<= x_opt{k}(1), q_opt>=0];     %current
    constraints = [constraints, x_opt{k}(2)<= 12+q1_opt, -12-q1_opt<= x_opt{k}(2), q1_opt>=0];%angular velocity
end
constraints=[constraints,InvSet.A*x_opt{N+1}<=InvSet.b]; %Terminal constraint
objective = objective + transpose(x_opt{N+1})*P2*(x_opt{N+1}); %adding the terminal cost at the end
controller2 = optimizer(constraints, objective,[],x_opt{1},[u_opt{:}]);

x = [6;15];
hold on
implementedU2 = [];
X2=[];
for i = 1:50
    X2=[X2,x];
    U2 = controller2{x};
    stairs(i:i+length(U2)-1,U2,'r')
    x = sys_d.A*x + sys_d.B*U2(1);
    implementedU2 = [implementedU2;U2(1)];
end

%% Q3
Q3 = diag([100,1000]); % State weighting matrix
[P3,gain3,L3] = idare(sys_d.A,sys_d.B,Q3,R); gain3 = -gain3; %for terminal cost(V_f)

% Control invariant set
system = LTISystem('A', sys_d.A, 'B', sys_d.B);
system.x.min = [-3; -12];
system.x.max = [3; 12];
system.u.min = -5;
system.u.max = 5;
system.x.penalty = QuadFunction(Q3);
system.u.penalty = QuadFunction(R);

system.x.with('terminalSet');
system.x.terminalSet = system.LQRSet();
system.x.with('terminalPenalty');
system.x.terminalPenalty = system.LQRPenalty();

InvSet = system.LQRSet();

% Set Optimizer
u_opt = sdpvar(repmat(nu,1,N),repmat(1,1,N));
x_opt= sdpvar(repmat(nx,1,N+1),repmat(1,1,N+1));
z_opt=sdpvar(1,1); %for constraints
q_opt=sdpvar(1,1); %for constraints
q1_opt=sdpvar(1,1); %for constraints
constraints = [];
objective = 0;

for k = 1:N
    % x{k+1} = sys_d.A*x{k} + sys_d.B*u{k};
    objective = objective + transpose(x_opt{k})*Q3*(x_opt{k}) + transpose(u_opt{k})*R*(u_opt{k}) +K*z_opt+ K*q_opt + K*q1_opt;
    constraints = [constraints, x_opt{k+1} == sys_d.A*x_opt{k} + sys_d.B*u_opt{k}];
    constraints = [constraints, u_opt{k}<= 5+z_opt, -5-z_opt<= u_opt{k}, z_opt>=0];           %voltage
    constraints = [constraints, x_opt{k}(1)<= 3+q_opt, -3-q_opt<= x_opt{k}(1), q_opt>=0];     %current
    constraints = [constraints, x_opt{k}(2)<= 12+q1_opt, -12-q1_opt<= x_opt{k}(2), q1_opt>=0];%angular velocity
end
constraints=[constraints,InvSet.A*x_opt{N+1}<=InvSet.b]; %Terminal constraint
objective = objective + transpose(x_opt{N+1})*P3*(x_opt{N+1}); %adding the terminal cost at the end
controller3 = optimizer(constraints, objective,[],x_opt{1},[u_opt{:}]);

x = [6;15];
hold on
implementedU3 = [];
X3=[];
for i = 1:50
    X3=[X3,x];
    U3 = controller3{x};
    stairs(i:i+length(U3)-1,U3,'r')
    x = sys_d.A*x + sys_d.B*U3(1);
    implementedU3 = [implementedU3;U3(1)];
end


%% Simulation
im_u=stairs(implementedU, 'g');
hold on
im_u1=stairs(implementedU1, 'r');
hold on
im_u2=stairs(implementedU2, 'b');
hold on
im_u3=stairs(implementedU3, 'k');
legend([im_u, im_u1, im_u2, im_u3], 'q=10', 'q=50', 'q=100', 'q=1000')
xlabel('Steps')
grid



figure()
x21=stairs(X(2,:),'g');
hold on
x22=stairs(X1(2,:),'r');
hold on
x23=stairs(X2(2,:),'b');
hold on
x24=stairs(X3(2,:),'k');
xlabel('Steps')
ylabel('RPM')
legend([x21, x22, x23, x24],'q=10', 'q=50', 'q=100', 'q=1000')
grid