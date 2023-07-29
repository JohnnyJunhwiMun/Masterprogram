%% Reference tracking with disturance and states observer
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
Ts=0.5;                    %Sampling period
sys_d=c2d(sys,Ts);         %Discretized system


%%Setting

%Augmented system
Bd=[1;0];
Cd=1;
A_aug=[sys_d.A Bd;zeros(1,2) 1];
C_aug=[sys_d.C Cd];

%Test controllability
C=ctrb(sys_d.A,sys_d.B); %Test controllability
C1=rank(C);

%Test observability
ob=obsv(sys_d.A,sys_d.C); %A,C observable
O1=rank(ob);
O2=rank([eye(2)-sys_d.A -Bd;sys_d.C Cd]); %Should be n+nd

%Observer design
pob=[0.2 0.3 0.15];%desired rate
L=place(A_aug',C_aug',pob)';%pole placement
L1=L(1:2);
L2=L(3);
rank(A_aug-L*C_aug) %A-LC stable

%Parameters
nx = 2; % Number of states
nu = 1; % Number of inputs

T=40;   % Simulation time
Q = diag([100,10]); % State weighting matrix
R = 1; % Input weighting matrix
N=7;   % Prediction horizon
K=10.^10;% Soften the constraint

y_ref=4; %Output reference
x=[0;0]; %Real states-Initial value [current, rotor speed]
x_ob=[0;0];%Initial observerd states
d=[ones(1,20),1.5*ones(1,T-20)];     %Real disturbance
d_ob=0;  %Observed disturbance-Initial value

[P,gain,L] = idare(sys_d.A,sys_d.B,Q,R); gain = -gain; %for terminal cost(V_f)

implementedU = [];
X=[];
X_ob=[];
D=[];
Y=[];

%% Calculate invariant set

system = LTISystem('A', sys_d.A, 'B', sys_d.B);
system.x.min = [-3;-10];
system.x.max = [3;10];
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

for i = 1:T

    X=[X,x];
    X_ob=[X_ob,x_ob];
    D=[D,[d_ob;d(i)]];

    % OTS-online
    x_ref=sdpvar(2,1);
    u_ref=sdpvar(1,1);
    objective_ref= transpose(x_ref)*Q*x_ref + transpose(u_ref)*R*u_ref;
    constraints_ref= [x_ref==sys_d.A*x_ref + sys_d.B*u_ref+Bd*d_ob,sys_d.C*x_ref==y_ref-Cd*d_ob];
    optimize(constraints_ref,objective_ref);
    x_ref=value(x_ref);
    u_ref=value(u_ref);



    % Set Optimizer
    u_opt = sdpvar(repmat(nu,1,N),repmat(1,1,N));
    x_opt = sdpvar(repmat(nx,1,N+1),repmat(1,1,N+1));
    z_opt=sdpvar(1,1); %for constraints
    q_opt=sdpvar(1,1); %for constraints
    q1_opt=sdpvar(1,1); %for constraints
    constraints = [];
    objective = 0;
    for k = 1:N
        %J=stagecost+finalcost
        objective = objective + transpose(x_opt{k}-x_ref)*Q*(x_opt{k}-x_ref) + transpose(u_opt{k}-u_ref)*R*(u_opt{k}-u_ref) +K*z_opt+ K*q_opt + K*q1_opt;
        constraints = [constraints, x_opt{k+1} == sys_d.A*x_opt{k}+sys_d.B*u_opt{k}+Bd*d_ob];
        constraints = [constraints, u_opt{k}<=3+z_opt,-3-z_opt<= u_opt{k}, z_opt>=0];           %voltage
        constraints = [constraints, x_opt{k}(1)<= 3+q_opt, -3-q_opt<= x_opt{k}(1), q_opt>=0];     %current
        constraints = [constraints, x_opt{k}(2)<= 10+q1_opt, -10-q1_opt<= x_opt{k}(2), q1_opt>=0];%angular velocity
    end
    constraints=[constraints,InvSet.A*x_opt{N+1}<=InvSet.b]; %Terminal constraint
    objective = objective + transpose(x_opt{N+1}-x_ref)*P*(x_opt{N+1}-x_ref); %adding the terminal cost at the end
    controller = optimizer(constraints, objective,[],x_opt{1},[u_opt{:}]);

    %Simulate with first input
    U = controller{x_ob};
    stairs(i:i+length(U)-1,U,'r')

    y=sys_d.C*x+d(i);     %Output update
    y_ob=sys.C*x_ob+d_ob; %Observed output update
    Y=[Y,[y;y_ref]];
    x = sys_d.A*x + sys_d.B*U(1)+Bd*d(i); %State update
    x_ob=sys_d.A*x_ob+sys_d.B*U(1)+Bd*d_ob+L1*(y-y_ob);%Observed state update
    d_ob=d_ob+L2*(y-y_ob); %Disturbance observer update
    pause(0.05)
    pre_u=stairs(i:i+length(U)-1,U,'k');hold on
    implementedU = [implementedU;U(1)];%Record all implemented inputs

end
im_u=stairs(implementedU,'b');
xlabel('Steps')
legend([pre_u,im_u],'Predctied inputs','Implemented inputs');
grid

figure()
stairs(X(1,:));hold on
stairs(X_ob(1,:));
xlabel('Steps')
ylabel('amp')
legend('Current','Observed current')
grid

figure()
stairs(X(2,:));hold on
stairs(X_ob(2,:));
xlabel('Steps')
ylabel('RPM')
legend('Rotor speed','Observed Rotor speed')
grid

figure()
stairs(D(1,:));hold on
stairs(D(2,:));
xlabel('Steps')
legend('Observed disturbance','Real disturbance')
grid

figure()
stairs(Y(1,:));hold on
stairs(Y(2,:));
xlabel('Steps')
ylabel('System output')
legend('Output','Reference output')
grid