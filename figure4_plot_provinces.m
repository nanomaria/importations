%------------------------------------------------------------------------
% AUTHORSHIP AND CONTACT
% Original code by Maria M. Martignoni (Memorial University of Newfoundland)
% Contact: mmartignonim@mun.ca
%-------------------------------------------------



close all
clear all

Places = ['St. Johns','Moncton','Halifax','Chalottentown','Whitehorse','Montréal-Pierre-Elliot-Trudeau, Quebec ', ...
     'Toronto International Airport T1','John Diefenbaker Airport, Saskatoon', 'Winnipeg','Vancouver airport', ...
     'Calgary','Iqaluit'];
 Labels = {'NL','NB','NS','PE','YT','QC','ON','SK','MB','BC','AB','NU'};
 Labels_c = {'QC','ON','MB','BC','AB'};
 Labels_c_s = {'NL','NB','NS','PE','YT','SK','NU'};
 
 
 m_pre = [328,75,1610,9,992,16938,24603,300,961,19000,7595,18];
 m_pre_c_s = [328,75,1610,9,992,300,18];
 m_pre_c = [16938,24603,961,19000,7595];

 m_post = [32,3, 67,  1,10,  1170,1801, 1,  21, 1016, 323, 2];
 
 L = [122, 61,333,30,38,1552,3004,323,693,666,1797,11];
 L_c_s = [122, 61,333,30,38,323,11];
 L_c = [1552,3004,693,666,1797];

 
 figure(1)
 subplot(1,2,1)
 plot(L_c,m_pre_c,'s','MarkerSize',8,...
    'MarkerEdgeColor','red',...
    'MarkerFaceColor',[1 .6 .6])
 text(L_c,m_pre_c,Labels_c,'VerticalAlignment','bottom','HorizontalAlignment','left')
 xlabel('L/h ratio')
 ylabel('International arrivals (pre-pandemic)')
 hold on
 plot(L_c_s,m_pre_c_s,'s','MarkerSize',6,...
    'MarkerEdgeColor','red',...
    'MarkerFaceColor',[1 .6 .6])
 axis([0 3500 0 27000])
set(gca,'FontSize',11)
subplot(1,2,2)
 plot(L_c_s,m_pre_c_s,'s','MarkerSize',8,...
    'MarkerEdgeColor','red',...
    'MarkerFaceColor',[1 .6 .6])
xlabel('L/h ratio')
 ylabel('International arrivals (pre-pandemic)')
set(gca,'FontSize',11)
axis([0 400 0 1800])
  text(L_c_s,m_pre_c_s,Labels_c_s,'VerticalAlignment','bottom','HorizontalAlignment','left')
 