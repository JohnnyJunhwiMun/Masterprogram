%% Comparison MPC and LQR

%% State space model of a DC motor

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
Ts=0.5;                    %Sampling period
Ds=c2d(sys,Ts);         %Discretized system     

%% LQR
Q   = diag([100,10]);
R   = 1;
K   = dlqr(Ds.A, Ds.B, Q, R);

Gcl = ss(Ds.A-Ds.B*K, Ds.B, Ds.C, Ds.D, Ts);  % closed-loop system
Ndc= 1/dcgain(Gcl);
[T,D] = eig(Ds.A-Ds.B*K);

%define the desired state 
xd=[0; 0];
% set initial condition
x0=[0,0];
%final simulation time 
tFinal=10;
time_total=0:Ts:tFinal;

% input_ud=ud*ones(size(time_total));
input_ud=ones(size(time_total));
%open-loop step response
[output_ud_only,time_ud_only,state_ud_only] = lsim(Ds,input_ud,time_total,x0);

closed_loop_input= ones(size(time_total));
% adding the reference
[output_closed_loop,time_closed_loop,state_closed_loop] = lsim(Ndc*Gcl,closed_loop_input,time_total,x0);




%% MPC

nx = 2; % Number of states
nu = 1; % Number of inputs

    
Q = diag([100, 10]); % State weighting matrix
R = 1; % Input weighting matrix
N=7;  % Prediction horizon


% OTS
y_ref=1;
x_ref=sdpvar(2,1);
u_ref=sdpvar(1,1);
objective_ref=transpose(x_ref)*Q*x_ref + transpose(u_ref)*R*u_ref;
constraints_ref= [x_ref==Ds.A*x_ref + Ds.B*u_ref,Ds.C*x_ref==y_ref];

optimize(constraints_ref,objective_ref);
x_ref=value(x_ref);
u_ref=value(u_ref);

% Calculate invariant set
system = LTISystem('A', sys_d.A, 'B', sys_d.B);
system.x.min = [-3; -12];
system.x.max = [3; 12];
system.u.min = -0.1;
system.u.max = 0.1;
system.x.penalty = QuadFunction(Q);
system.u.penalty = QuadFunction(R);

system.x.with('terminalSet');
system.x.terminalSet = system.LQRSet();
system.x.with('terminalPenalty');
system.x.terminalPenalty = system.LQRPenalty();

InvSet = system.LQRSet();

% Optimizer setup
u = sdpvar(repmat(nu,1,N),repmat(1,1,N));
x = sdpvar(repmat(nx,1,N+1),repmat(1,1,N+1));
[P,gain,L] = idare(Ds.A,Ds.B,Q,R); gain = -gain;
constraints = [];
objective = 0;
for k = 1:N
 objective = objective + transpose(x{k}-x_ref)*Q*(x{k}-x_ref) + transpose(u{k}-u_ref)*R*(u{k}-u_ref);
 constraints = [constraints, x{k+1} == Ds.A*x{k} + Ds.B*u{k}];
 constraints = [constraints, -0.1<= u{k}<= 0.1, -3<=x{k+1}(1)<=3, -12<=x{k+1}(2)<=12];
end
constraints=[constraints,InvSet.A*x_opt{N+1}<=InvSet.b]; %Terminal constraint
objective = objective + transpose(x{N+1}-x_ref)*P*(x{N+1}-x_ref); %adding the terminal cost at the end

controller = optimizer(constraints, objective,[],x{1},[u{:}]);

%Simulation
x = [0;0];
implementedU = [];
X=[];
for i = 1:20
  X=[X,x];
  U = controller{x};  
  stairs(i:i+length(U)-1,U,'r')
  x = Ds.A*x + Ds.B*U(1);
  pause(0.05)
  stairs(i:i+length(U)-1,U,'k');hold on
  implementedU = [implementedU;U(1)];
end
%Visualization
stairs(implementedU,'b');
figure()
stairs(X(2,:))
hold on
stairs(output_closed_loop)
xlabel('Steps')
ylabel('System output')
legend('MPC','LQR-closed loop')
grid