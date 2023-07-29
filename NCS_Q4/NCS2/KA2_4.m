clc
clear

syms h tau

a = 5;
bl = 3;
c = 1;
A = [a-bl, 0.5-c; 0, 1];
B = [0; 1];
p1 = [-1, -2];
K = place(A, B, p1);


F_x=expm(A*h);
F_u=(expm(A*h)-expm(A*(h-tau)))/A*B;
G1=(expm(A*(h-tau))-eye(2))/A*B;
F=[F_x F_u;zeros(1,2) 0];
F=simplify(F)
G=[G1;1];
G=simplify(G)


a1=exp(-2*(h-tau));
a2=exp(h-tau);

F_1=zeros(3,3);
F_2=zeros(3,3);
F_1(1,3)=-1/12;
F_2(1,3)=-1/6;
F_2(2,3)=-1;
F_0=simplify(F-a1*F_1-a2*F_2)

G_1=[1/12;0;0];
G_2=[1/6;1;0];
G_0=simplify(G-a1*G_1-a2*G_2)

h_val=0.2;
val=[];
for tau_val=0:0.01:h_val
F_val=double(subs(F, [h, tau], [h_val,tau_val]));
val=[val;[F_val(1,3),F_val(2,3)]];
end
figure()
plot(val(:,1),val(:,2));hold on
%%


a1_min = double(subs(a1,[h,tau],[h_val,0]))
a1_max = double(subs(a1,[h,tau],[h_val,h_val]))
a2_min = double(subs(a2,[h,tau],[h_val,h_val]))
a2_max = double(subs(a2,[h,tau],[h_val,0]))


F_over_approx{1}=double(subs(F_0+a1_max*F_1+a2_max*F_2,h,h_val));
F_over_approx{2}=double(subs(F_0+a1_min*F_1+a2_max*F_2,h,h_val));
F_over_approx{3}=double(subs(F_0+a1_min*F_1+a2_min*F_2,h,h_val));
F_over_approx{4}=double(subs(F_0+a1_max*F_1+a2_min*F_2,h,h_val));

G_over_approx{1}=double(subs(G_0+a1_max*G_1+a2_max*G_2,h,h_val));
G_over_approx{2}=double(subs(G_0+a1_min*G_1+a2_max*G_2,h,h_val));
G_over_approx{3}=double(subs(G_0+a1_min*G_1+a2_min*G_2,h,h_val));
G_over_approx{4}=double(subs(G_0+a1_max*G_1+a2_min*G_2,h,h_val));

for i=1:1:length(F_over_approx)
scatter(F_over_approx{i}(1,3),F_over_approx{i}(2,3),'green','filled');hold on
vertices_val(i,:)=[F_over_approx{i}(1,3),F_over_approx{i}(2,3)];
end
p=patch(vertices_val(:,1),vertices_val(:,2),'m');
uistack(p, 'bottom');