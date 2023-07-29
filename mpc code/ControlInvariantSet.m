%% Plot control invariant set

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

Q = diag([100,10]); % State weighting matrix
R = 1; % Input weighting matrix

sys1=ss(A,B,C,D);           %Denominator and Nominator of Transfer function
Ts=0.5;                    %Sampling period
Ds=c2d(sys1,Ts);         %Discretized system     

system = LTISystem('A', Ds.A, 'B', Ds.B);
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

InvSet = system.LQRSet()

InvSet.plot()
