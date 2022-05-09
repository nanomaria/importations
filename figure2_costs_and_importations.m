%------------------------------------------------------------------------
% AUTHORSHIP AND CONTACT
% Original code by Maria M. Martignoni (Memorial University of Newfoundland)
% Contact: mmartignonim@mun.ca
%-------------------------------------------------



clear all
close all

global m gamma betamax L h 

% model parameters

for z = 1:4
    
choice = z;
m = [0.2;1];
L = [20;500];
if choice == 1  %L low m low
m = m(1);
L = L(1);
elseif choice == 2 % L high m low
m = m(2);
L = L(1);
elseif choice == 3 % L low m high
m = m(1);
L = L(2);
elseif choice == 4 % L high m high
m = m(2);
L = L(2);
end

h = 1;
c=0.05;
gamma = 1/14;
betamax = 2.4*gamma;

i0 = 0;
lambdaN = 0;


Nsteps2 = 500;  
Tf = 180;


%Step size
hs = 0.001;
time = 0:hs:(Nsteps2);
Nsteps = Nsteps2/hs;

%Initial guess for beta_opt
beta_opt = ones(1,Nsteps+1)*1*gamma;
beta_opt1= ones(1,Nsteps+1)*2.5*gamma;

I = zeros(1,Nsteps+1);
lambda = zeros(1,Nsteps+1);

%Number of iteraction to find the optimum (at some point I will put a while
%loop, stopping when convergence is observed)
Msteps = 100;

%Iteractions to convergence
%for z = 1:Msteps

 while ((beta_opt1 - beta_opt)/beta_opt1) > 0.001

%%RK to solve I
for i = 1:Nsteps
    h2 = hs/2;
I(1) = i0;
k1 = (beta_opt(i)- gamma)*I(i) + m;
k2 = ((beta_opt(i)+beta_opt(i+1))/2- gamma)*(I(i)+h2*k1) + m;
k3 = ((beta_opt(i)+beta_opt(i+1))/2- gamma)*(I(i)+h2*k2) + m;
k4 = ((beta_opt(i+1))- gamma)*(I(i)+hs*k3) + m;
I(i+1) = I(i) + hs/6*(k1+2*k2+2*k3+k4) ;
end



%RK backward in time to solve for lambda
for i = 1:(Nsteps)
j = Nsteps + 2 - i;
h2 = hs/2;
k1 = (beta_opt(j))*(1 - (lambda(j)))+lambda(j)*gamma;
k2 = (beta_opt(j)+beta_opt(j-1))/2*(1-(lambda(j)-h2*k1))+(lambda(j)-h2*k1)*gamma;
k3 = (beta_opt(j)+beta_opt(j-1))/2*(1-(lambda(j)-h2*k2))+(lambda(j)-h2*k2)*gamma;
k4 = (beta_opt(j-1))*(1 - (lambda(j)-hs*k3))+(lambda(j)-hs*k3)*gamma;
lambda(j-1) = lambda(j) - (hs/6)*(k1+2*k2+2*k3+k4);
end


%save old beta vector
beta_opt_old = beta_opt;

%update beta_opt vector
for i = 1:(Nsteps+1)
beta_opt1(i) = min(betamax,max(0,betamax*(1-I(i)*betamax*(h/(2*L))*(1-lambda(i)))));
end


%write new beta vector
beta_opt  = c*(beta_opt1) + (1-c)*beta_opt_old;

%plots (to see what the code is doing at each iteraction)
%figure(2)
%subplot(3,1,1)
%plot(time,beta_opt/gamma,'r.')
%hold on
%subplot(3,1,2)
%plot(time,I,'k.')
%hold on
%subplot(3,1,3)
%plot(time,lambda,'b.')
%hold on
end


%%If this number is larger than 0.001, then no convergence is observed
disp((beta_opt1 - beta_opt)/beta_opt1)
conv = (beta_opt1 - beta_opt)/beta_opt1;


%if convergence is observed, we plot the optimum
if conv < 0.001    
figure(1)
subplot(3,1,1)
plot(time,beta_opt/gamma,'r.')
ylabel('\beta^*/\gamma')
subplot(3,1,2)
plot(time,I,'k.')
ylabel('I')
subplot(3,1,3)
plot(time,lambda,'b.')
ylabel('\lambda')
xlabel('Time')

end



figure(2)

Imax = max(I(1:Tf/hs));

subplot(1,2,1)
if choice == 1
plot(time,beta_opt/gamma,'r:','Linewidth',2)
elseif choice == 2
plot(time,beta_opt/gamma,'r:','Linewidth',1)
elseif choice == 3 
plot(time,beta_opt/gamma,'r--','Linewidth',2)
elseif choice == 4 
plot(time,beta_opt/gamma,'r--','Linewidth',1)
end
hold on
xlabel('time [days]')
ylabel('R_0 = \beta^*/\gamma')
axis([ 0 Tf 0.65 2.5])
set(gca,'fontsize',12)


subplot(1,2,2)
if choice == 1
plot(time,I,'k:','Linewidth',2)
elseif choice == 2 
plot(time,I,'k:','Linewidth',1)
elseif choice == 3 
plot(time,I,'k--','Linewidth',2)
elseif choice == 4
plot(time,I,'k--','Linewidth',1)
end
axis([ 0 Tf 0 Imax+0.1*Imax])
set(gca,'fontsize',12)
hold on
ylabel('I(t)')
xlabel('time [days]')



%% Compute cost of optimal strategy
%cost at each point in time
u_p = L/h*(1-beta_opt/betamax).^2;

%cost cumulative
u = zeros(1,(Tf/hs));
for i = 1:Tf/hs
    u(1) = u_p(1) ; 
    u(i+1) = u(i)+ u_p(i)*hs;
end
u = u(1:Tf/hs);


%% costs of infections (calculate number of cumulative cases)
C = zeros(1,(Tf/hs));
for i = 1:Tf/hs
C(1) = I(1);
C(i+1) = C(i) + (beta_opt(i)*I(i) + m)*hs;
end
C = C(1:Tf/hs);
%close all
%figure(4); plot(time(1:Tf/hs),I(1:Tf/hs));hold on;plot(time(1:Tf/hs),C(1:Tf/hs));

GDPu = zeros(1,Tf/hs);
GDPi = zeros(1,Tf/hs);
for i = 1:(Tf/hs)
GDPu(i) = 10*u(i)/(L/h*(i*hs));
GDPi(i) = 10*C(i)/(L/h*(i*hs));
end

if choice ==1
GDPu1 = GDPu;
GDPi1 = GDPi;
elseif choice ==2
GDPu2 = GDPu;
GDPi2 = GDPi;
elseif choice ==3
GDPu3 = GDPu;
GDPi3 = GDPi;
elseif choice ==4
GDPu4 = GDPu;
GDPi4 = GDPi;
end


end


%% Figures
figure(3)
fig2 = figure(3);


x2 = categorical({'L=20, m=0.2','L=20, m=1','L=500, m=0.2','L=500, m=1'});
x2 = reordercats(x2,string(x2));


GDPl = [GDPu1(Tf/hs) GDPi1(Tf/hs) GDPu1(Tf/hs)+GDPi1(Tf/hs) ; GDPu2(Tf/hs) GDPi2(Tf/hs) GDPu2(Tf/hs)+GDPi2(Tf);GDPu3(Tf/hs) GDPi3(Tf/hs) GDPu3(Tf/hs)+GDPi3(Tf/hs);GDPu4(Tf/hs) GDPi4(Tf/hs) GDPu4(Tf/hs)+GDPi4(Tf/hs)];
set(gca,'fontsize',12)

b3 = bar(x2,GDPl,0.7,'FaceColor','flat');
b3(1).CData(1,:) = [1 0 0]; % group 1 1st bar
b3(2).CData(1,:) = [0 0 0]; % group 1 2nd bar
b3(3).CData(1,:) = [0 0 1]; % group 1 3rd bar
b3(1).CData(2,:) = [1 0 0]; % group 2 1st bar
b3(2).CData(2,:) = [0 0 0]; 
b3(3).CData(2,:) = [0 0 1]; 
b3(1).CData(3,:) = [1 0 0]; 
b3(2).CData(3,:) = [0 0 0]; 
b3(3).CData(3,:) = [0 0 1]; 
b3(1).CData(4,:) = [1 0 0]; 
b3(2).CData(4,:) = [0 0 0]; 
b3(3).CData(4,:) = [0 0 1]; 
ylabel('Cost [% GDP loss]')

A = [GDPu1(Tf/hs)+GDPi1(Tf/hs),GDPu2(Tf/hs)+GDPi2(Tf),GDPu3(Tf/hs)+GDPi3(Tf/hs),GDPu4(Tf/hs)+GDPi4(Tf/hs)];
ylimt = max(A)+0.1*max(A);
ylim([ 0 ylimt]);
set(gca,'fontsize',12)


