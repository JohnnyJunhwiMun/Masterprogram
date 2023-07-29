%% Reference tracking with disturbance
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

C=ctrb(sys_d.A,sys_d.B);   %Test controllability
C1=rank(C);

%% Setting

nx = 2; % Number of states
nu = 1; % Number of inputs
T=40;

Q = diag([100,10]); % State weighting matrix
R = 1; % Input weighting matrix
N=7;   %Prediction horizon
K=10.^10;%soften the constraint

y_ref=4; %Reference 
x=[0;0]; %Initial states
d=[ones(1,20),1.5*ones(1,T-20)];     %Real disturbance
d_ob=0;  %Observed disturbance-Initial value
L=0.5;   %Observer gain

[P,gain,L1] = idare(sys_d.A,sys_d.B,Q,R); gain = -gain;  %For terminal cost

implementedU = [];
X=[];
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
  y=sys_d.C*x+d_ob;
  Y=[Y,y];
  X=[X,x];
  D=[D,[d_ob;d(i)]];


  % OTS-Online
  x_ref=sdpvar(2,1);
  u_ref=sdpvar(1,1);
  objective_ref=transpose(x_ref)*Q*x_ref + transpose(u_ref)*R*u_ref;
  constraints_ref= [x_ref==sys_d.A*x_ref + sys_d.B*u_ref,sys_d.C*x_ref==y_ref-d_ob];%d_ob will iterate
  optimize(constraints_ref,objective_ref);
  x_ref=value(x_ref);
  u_ref=value(u_ref);

  % Set Opimizer
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
        constraints = [constraints, x_opt{k+1} == sys_d.A*x_opt{k}+sys_d.B*u_opt{k}];
        constraints = [constraints, u_opt{k}<=3+z_opt,-3-z_opt<= u_opt{k}, z_opt>=0];           %voltage
        constraints = [constraints, x_opt{k}(1)<= 3+q_opt, -3-q_opt<= x_opt{k}(1), q_opt>=0];     %current
        constraints = [constraints, x_opt{k}(2)<= 10+q1_opt, -10-q1_opt<= x_opt{k}(2), q1_opt>=0];%angular velocity
    end
    constraints=[constraints,InvSet.A*x_opt{N+1}<=InvSet.b]; %Terminal constraint
    objective = objective + transpose(x_opt{N+1}-x_ref)*P*(x_opt{N+1}-x_ref); %adding the terminal cost at the end
    controller = optimizer(constraints, objective,[],x_opt{1},[u_opt{:}]);
  
  %Simulate with first input
  U = controller{x};  
  stairs(i:i+length(U)-1,U,'r')
  x = sys_d.A*x + sys_d.B*U(1);         %State update
  d_ob=(1-L)*d_ob+L*d(i);                  %Disturbance observer update
  pause(0.05)
  pre_u=stairs(i:i+length(U)-1,U,'k');hold on
  implementedU = [implementedU;U(1)];   %Record all implemented inputs
end

im_u=stairs(implementedU,'b');
legend([pre_u,im_u],'Predicted inputs','Implemented inputs');
grid
figure()
stairs(X(1,:));hold on
stairs(X(2,:))
legend('Current','Rotor speed','Location','northwest')
grid
figure()
stairs(D(1,:));hold on;
stairs(D(2,:));
legend('Observed disturbance','Real disturbance','Location','northwest')
grid
figure()
stairs(Y);
legend('System output','Location','northwest')
grid