%------------------------------------------------------------------------
% AUTHORSHIP AND CONTACT
% Original code by Maria M. Martignoni (Memorial University of Newfoundland)
% Contact: mmartignonim@mun.ca
%-------------------------------------------------


clear all
close all

global m gamma betamax L h 

Lv = [20;100;500;1000];
mv = 0.02:0.02:1.2;

GDPu = zeros(length(mv),length(Lv));
GDPi = zeros(length(mv),length(Lv));
GDPt = zeros(length(mv),length(Lv));
 Iss = zeros(length(mv),length(Lv));
betass = zeros(length(mv),length(Lv));

for z = 1:length(Lv)
    
L = Lv(z);

for w = 1:length(mv)

m = mv(w);

h = 1;
c=0.05;
gamma = 1/14;
betamax = 2.4*gamma;

i0 = 0;
lambdaN = 0;


Nsteps2 = 500;  


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

end



Tf = Nsteps2/2;

Iss(w,z) = I(Tf/hs);
betass(w,z) = beta_opt(Tf/hs);

u_p = L/h*(1-beta_opt/betamax).^2;
%cost cumulative
u = zeros(1,(Tf/hs));
for i = 1:Tf/hs
    u(1) = u_p(1) ; 
    u(i+1) = u(i)+ u_p(i)*hs;
end
u = u(Tf/hs);


%% costs of infections (calculate number of cumulative cases)
C = zeros(1,(Tf/hs));
for i = 1:Tf/hs
C(1) = I(1);
C(i+1) = C(i) + (beta_opt(i)*I(i) + m)*hs;
end
C = C(Tf/hs);

GDPu(w,z) = 10*u/(L/h*(Tf));
GDPi(w,z) = 10*C/(L/h*(Tf));
GDPt(w,z) = GDPu(w,z)+GDPi(w,z);

end
end

%% Figures
figure(1)
subplot(1,2,1)
plot(mv,betass/gamma,'r-','linewidth',1)
xlabel('Importation rate m')
%ylabel('Optimal transmission rate (at equilibrium) [ R_0 = \beta^*/\gamma]')
set(gca,'FontSize',11)
axis([0 max(mv) min(min(betass))/gamma*(0.9) max(max(betass))/gamma*(1.1)])

subplot(1,2,2)
plot(mv,Iss,'k-','linewidth',1)
xlabel('Importation rate m')
%ylabel('Number of active cases (at equilibrium) [I(t)]')
set(gca,'FontSize',11)
axis([0 max(mv) 0 max(max(Iss))*1.1])

figure(2)
subplot(1,3,1)
plot(mv,GDPu,'r-','linewidth',1)
xlabel('Importation rate m')
%ylabel('Daily cost (at equilibrium) [% GDP loss]')
set(gca,'FontSize',11)

subplot(1,3,2)
plot(mv,GDPi,'k-','linewidth',1)
xlabel('Importation rate m')
%ylabel('Daily cost (at equilibrium) [% GDP loss]')
set(gca,'FontSize',11)

subplot(1,3,3)
plot(mv,GDPt,'b-','linewidth',1)
xlabel('Importation rate m')
%ylabel('Daily cost (at equilibrium) [% GDP loss]')
set(gca,'FontSize',11)


