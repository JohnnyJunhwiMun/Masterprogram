%% Standard MPC with no reference and disturbance
yalmip('clear')
clear all

%% Model data
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

C=ctrb(sys_d.A,sys_d.B); %Check controllability
rank(C)
%% Calculate invariant set
system = LTISystem('A', sys_d.A, 'B', sys_d.B);
system.x.min = [-3; -10];
system.x.max = [3; 10];
system.u.min = -3;
system.u.max = 3;
system.x.penalty = QuadFunction(diag([100,10]));
system.u.penalty = QuadFunction(1);

system.x.with('terminalSet');
system.x.terminalSet = system.LQRSet();
system.x.with('terminalPenalty');
system.x.terminalPenalty = system.LQRPenalty();

InvSet = system.LQRSet();

%% Simulation
nx = 2; % Number of states
nu = 1; % Number of inputs


Q = diag([100,10]); % State weighting matrix
R = 1; % Input weighting matrix
N=7;
K=10.^10;%soften the constraint
[P,gain,L] = idare(sys_d.A,sys_d.B,Q,R); gain = -gain; %for terminal cost(V_f)


% Set Optimizer
u_opt = sdpvar(repmat(nu,1,N),repmat(1,1,N));
x_opt = sdpvar(repmat(nx,1,N+1),repmat(1,1,N+1));
z_opt=sdpvar(1,1); %for constraints
q_opt=sdpvar(1,1); %for constraints
q1_opt=sdpvar(1,1); %for constraints
constraints = [];
objective = 0;

for k = 1:N
    objective = objective + transpose(x_opt{k})*Q*(x_opt{k}) + transpose(u_opt{k})*R*(u_opt{k}) +K*z_opt+ K*q_opt + K*q1_opt;
    constraints = [constraints, x_opt{k+1} == sys_d.A*x_opt{k} + sys_d.B*u_opt{k}];
    constraints = [constraints, u_opt{k}<= 3+z_opt, -3-z_opt<= u_opt{k}, z_opt>=0];           %voltage
    constraints = [constraints, x_opt{k}(1)<= 3+q_opt, -3-q_opt<= x_opt{k}(1), q_opt>=0];     %current
    constraints = [constraints, x_opt{k}(2)<= 10+q1_opt, -10-q1_opt<= x_opt{k}(2), q1_opt>=0];%angular velocity
end
constraints=[constraints,InvSet.A*x_opt{N+1}<=InvSet.b]; %Terminal constraint
objective = objective + transpose(x_opt{N+1})*P*(x_opt{N+1}); %adding the terminal cost at the end
controller = optimizer(constraints, objective,[],x_opt{1},[u_opt{:}]);

x = [6;12];
hold on
implementedU = [];
X=[];
for i = 1:50
    X=[X,x];
    U = controller{x};
    stairs(i:i+length(U)-1,U,'r')
    x = sys_d.A*x + sys_d.B*U(1);
    pause(0.1)
    stairs(i:i+length(U)-1,U,'k')
    pre_u=stairs(i:i+length(U)-1,U,'k');
    implementedU = [implementedU;U(1)];
end

%% Plot
im_u=stairs(implementedU,'b');
xlabel('Steps')
legend([pre_u,im_u],'Predctied inputs','Implemented inputs');
grid


figure()
stairs(X(1,:))
xlabel('Steps')
ylabel('system states')
hold on;
stairs(X(2,:))
legend('Current','Rotor speed')
grid

figure()
Lyap1=zeros(1,size(X,2)-1);
Lyap2=zeros(1,size(X,2)-1);
for i=2:1:size(X,2)
    Lyap1(i-1)=transpose(X(:,i))*P*X(:,i);
    Lyap2(i-1)=transpose(X(:,i-1))*P*X(:,i-1)-(transpose(X(:,i-1))*Q*(X(:,i-1)) + transpose(implementedU(i-1))*R*implementedU(i-1));
end
stairs(Lyap1);hold on
stairs(Lyap2);
legend('Vf(Ax+Bu)','Vf(x)-l(x,u)');
grid

