%------------------------------------------------------------------------
% AUTHORSHIP AND CONTACT
% Original code by Maria M. Martignoni (Memorial University of Newfoundland)
% Contact: mmartignonim@mun.ca
%-------------------------------------------------



close all
clear all
global betal gamma Id I0 betamax betaeq

gamma = 1/14;
betal = 0.3*gamma;
betamax = 2.4*gamma;
betaeq = 0.9*gamma;

Idv = [1,10,20,50];




for z = 1:length(Idv)

 Id = Idv(z);
 I0v = Id:1:500;
T = zeros(1,length(I0v));
  
for i = 1:length(I0v)

I0 = I0v(i);
T(i) = (log(Id/I0)*(betamax-betal)^2)/((betal-gamma) * (betamax-betaeq)^2);


end


figure(1)
if Id == 1
plot(I0v,T,'k-','Linewidth',1.5)
else
    plot(I0v,T,'k--','Linewidth',0.5)
end
xlabel('Number of active cases at lockdown time [I_0]')
ylabel('Days with no cases escaping self-isolation [T = 1/m]')
%axis([0 1 0 100])
set(gca,'fontsize',12)
hold on



end




