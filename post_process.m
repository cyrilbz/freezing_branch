% % This file is part of "Freezing Branch", a code that simulates 
% % freeze induced cell dehydration, pressure and diameter changes 
% % in a tree stem.
% % 
% % Author: Cyril Bozonnet (cyril.bozonnet@inrae.fr; github: cyrilbz) 
% %         INRAE, PIAF, Clermont-Ferrand
% %         
% %         The code structure has been inspired from an existing code 
% %         written by Isabell Graf (Konrad) and John M. Stockie
% %         Department of Mathematics
% %         Simon Fraser University
% %          
% % Developped using Matlab version R2018a.
% % 
% % The code is distributed under the CeCILL-B free software 
% % license agreement.
% % (https://cecill.info/licences/Licence_CeCILL-B_V1-en.html)

clear all
close all

%% Change directory if needed
olddr= pwd ;
mydr = strcat(pwd,'\base_case_2') ;
cd(mydr)
%% plot program
 p = parameters ; % load parameters and functions
 
load ('t.mat') % load time vector
load ('H.mat') % load H matrix
load ('rf.mat') % load rf matrix
load ('rv.mat') % load rv matrix
load ('Uroot.mat') % load Root matrix
load ('Ufv.mat') % load fiber-vessel exchange matrix
load ('Upv.mat') % load parenchyma-vessel exchange matrix
load ('Upp.mat') % load parenchyma-parenchyma exchange matrix
load ('Vp.mat') % load parenchyma volume matrix
load ('pp.mat') % load parenchyma pressure matrix
load ('nsv.mat')
load ('nslc.mat')
load ('ngv.mat')
load('Vbark_cell.mat','Vbark_cell') % load bark volume matrix
load('pbark.mat','pbark') % load pressure bark matrix
load('nsb.mat')


%% Reconstruct vessel water pressure signal and other state variables

nt=length(t) ; % stored times
pwv=zeros(size(rf)) ; % Initialize water pressure matrix

for i=1:nt
    % Get temperature and related informations
        % get microscale sugar concentration
        Cs_v(i,:) = nsv(i,:)./(pi*(p.Rv^2-rv(i,:).^2)*p.Lv) ; % This is on the vessel grid !
        % Upscale it (interpolate on the fine grid + fill end points using nearest value
        Cs_v_nc = interp1(p.ru,Cs_v(i,:)',p.r,'linear') ;
        Cs_v_nc = fillmissing(Cs_v_nc,'linear','EndValues', 'nearest') ;
        [Tm, Tirv, Twrv, Hirv, Hwrv, Tirf, Twrf, Hirf, Hwrf] = p.RegParam(Cs_v_nc) ;
        Temp = p.HtoT(H(i,:), Tirv, Twrv, Hirv, Hwrv, Tirf, Twrf, Hirf, Hwrf) ; % converts enthalpy to temperature
        Temp_nv = Temp(1:p.MR:end) ; % Temperature on the coarser grid
        Temp_save (i,:) = Temp ; % Temperature to be saved on the vessel grid 
        Text = p.Tempout(t(i)) ; % Compute external temperature
        H_nv = H(i,1:p.MR:end) ;
        Hwrv_nv = Hwrv(1:p.MR:end) ; % Vessel water enthalpy on the coarse grid
        Tm_nv = Tm(1:p.MR:end) ; % Melting temperature on the coarse grid
        fv = min(1,max((Hwrv_nv-H_nv)/p.DHv,0)) ;
      
    % Get gas pressure
        pgv = ngv(i,:)*p.Rg.*Temp_nv./(pi*rv(i,:).^2*p.Lv) ; % Vessel gas pressure
    
    % get liquid pressure
        pwv(i,:) = pgv - p.saw./rv(i,:) ; % Vessel water pressure
    
    % sugar concentrations & related quantities
        Cs_p(i,:) = nslc(i,:)./(Vp(i,:)-p.Vp_bound) ;
        Cs_bark(i) = nsb(i)./(Vbark_cell(i)-p.Vbark_bound) ;
        mu(i,:) = 0.09607e-3*exp(2.9*(p.Tc./Temp_nv(:)').^3).*(1+p.e_nu*Cs_p(i,:)/p.rhow.*exp(((Cs_p(i,:)/p.rhow).^p.f_nu)./(p.g_nu*Temp_nv(:)'/p.Tc+p.h_nu))) ;
        mu_bark(i) = 0.09607e-3*exp(2.9*(p.Tc./Text).^3).*(1+p.e_nu*Cs_bark(i)/p.rhow.*exp(((Cs_bark(i)/p.rhow).^p.f_nu)./(p.g_nu*Text/p.Tc+p.h_nu))) ;
        mean_mu = p.mean_nv*mu(i,:)' ; % averaged viscosity
        mean_mu(end) = 0.5*(mu(i,end) + mu_bark(i)) ; % boundary condition
        Tm_bark(i) = p.Tc - Cs_bark(i)*p.Kb/p.rhow ;
        Tm_par(i,:) = p.Tc - Cs_p(i,:)*p.Kb/p.rhow ;
        FH_bark(i) = Tm_bark(i) / Text ; % Frost hardiness based on sugar concentration - BARK
        FH_vac(i,:) = Tm_par(i,:) ./ Temp_nv ; % Frost hardiness based on sugar concentration - VACs
    
    % Compute osmotic and ice pressure
        posm_pv(i,:) = (Cs_p(i,:)'-Cs_v(i,:)')*p.Rg.*Temp_nv' ;
        posm_p(i,:) = Cs_p(i,:)'*p.Rg.*Temp_nv' ;
        pice_v(i,:) = p.rhow*(p.Ew-p.Ei).*log(Temp_nv./p.Tc).*fv ;
        posm_v(i,:) = Cs_v(i,:)'*p.Rg.*Temp_nv' ;
        posm_b(i) = Cs_bark(i)*p.Rg*Text ;
    
    % Compute conductivites
        Kpv(i,:) = p.PC./mu(i,:)*p.Apv*p.muw ;
        Kpp(i,:) = p.RayCond./mean_mu*p.muw ;
end

% Compute global quantities
    Text = p.Tempout(t) ; % Compute external temperature
    dUpp = [zeros(1,p.nv) ; diff(Upp,1,1)/(t(2)-t(1))] ; % P-P flow rate
    Qbx = dUpp(:,end)*p.Nray*p.Nstack ; % bark-xylem flow rate

%% Compute diameter variations

% xylem :
Vtot_par = sum(Vp,2)*p.Nstack*p.Nray; % sum of VAC volume for all vessels and rays
dVpar = Vtot_par-Vtot_par(1); % variations of VACs volume
DDx = sqrt((2*p.Rheart+2*p.L)^2 + 4/pi/p.Lv*dVpar) ; % xylem diameter variations

% bark :
Vbark_total_0=pi*(-p.bt^2 + 2*p.bt*p.Rtree)*p.Lv ; % total intial bark volume
Vfinal = Vbark_total_0*(1-1*p.BFF*(1-(Vbark_cell)/(p.Vbark_cell0))) ; % computed bark volume
D = sqrt(DDx.^2+4/pi/p.Lv*Vfinal) ; % total diameter from volume variations and xylem diameter fluctuations

% expe : 
Delta_D_am = -0.004*(2*p.Rtree*1000)^(1.275) ;

%% load experimental diameter variations

data_expe = importdata('../data_expe.txt') ;

time_expe = data_expe.data(:,3) ; 
d_expe = data_expe.data(:,2) ; 

%% Plots
pvw_ave = sum(pwv,2)/(p.nv)/1000 ;

figure(1)
%plot(t/3600,pvw_ave)
set (gcf,'color','white')
hold on
plot(t/3600,pvw_ave,'k','Linewidth',1.5)
p_save = pvw_ave ; % to save 
th = t/3600 ;
xlabel('time [h]','interpreter','latex','FontSize',15)
ylabel ('Mean vessel pressure [kPa]','interpreter','latex','FontSize',15)
ax = gca;
ax.FontSize = 16; 
ax.TickLabelInterpreter = 'latex' ;
box on
% my_legend = {''} ;
% legend(my_legend,'interpreter','latex','position',[0.4, 0.65, 0.3, 0.25])
% legend boxoff
text(0.025,0.95,'b)','Units','normalized','FontSize',14,'interpreter','latex')
%print -depsc fig3-mean_pv.eps

figure(2)
% Space time plot of vessel water pressure or other quantity
ru = (p.Rheart + p.dru) : p.dru : (p.Rtree-p.dru-p.bt) ;
[X,time]=meshgrid(ru,th) ;
s = surf(X,time,pwv,'Edgecolor','none');
hold on
% map = [0 0 0 
%     1 1 0];
colorbar
%colormap(gray)
view(2)
xlabel('r [m]')
ylabel ('Time [h]')


figure(3)
set (gcf,'color','white')
plot(th,pwv(:,1)/1000,'Color',[0  1 1.0000],'Linewidth',1.5)
hold on
plot(th,pwv(:,end)/1000,'Color',[0.5 0.5 1.0000],'Linewidth',1.5)
xlabel('time [h]','interpreter','latex','FontSize',15)
ylabel ('Vessel water pressure [kPa]','interpreter','latex','FontSize',15)
box on
my_legend = {'Near pith vessel','Near bark vessel'} ;
legend(my_legend,'interpreter','latex','location','best','Position',[0.38 0.3 0.3 0.15])
ax = gca;
ax.FontSize = 16; 
ax.TickLabelInterpreter = 'latex' ;
legend boxoff
text(0.025,0.95,'a)','Units','normalized','FontSize',14,'interpreter','latex')
%print -depsc fig3-pv.eps

figure(4)
set (gcf,'color','white')
plot(th,pp(:,1)/1e6,'Color',[0  1 1.0000],'Linewidth',1.5)
hold on
plot(th,pp(:,end)/1e6,'Color',[0.5 0.5 1.0000],'Linewidth',1.5)
plot(th(1:50:end),pbark(1:50:end)/1e6,'o','Color',[1 0 1.0000],'Linewidth',1)
xlabel('time [h]','interpreter','latex','FontSize',15)
ylabel ('Turgor pressure [MPa]','interpreter','latex','FontSize',15)
box on
my_legend = {'Near pith VACs','Near bark VACs','Bark cells'} ;
legend(my_legend,'interpreter','latex','location','best','Position',[0.38 0.4 0.3 0.2])
ax = gca;
ax.FontSize = 16; 
legend boxoff
ax.TickLabelInterpreter = 'latex' ;
text(0.025,0.95,'c)','Units','normalized','FontSize',14,'interpreter','latex')
%print -depsc fig3-turg.eps

figure(5)
set (gcf,'color','white')
plot(th,posm_p(:,1)/1e6,'Color',[0  1 1.0000],'Linewidth',1.5)
hold on
plot(th,posm_p(:,end)/1e6,'Color',[0.5 0.5 1.0000],'Linewidth',1.5)
plot(th,posm_b/1e6,'Color',[1 0 1.0000],'Linewidth',1.5)
plot(th,-pice_v(:,1)/1e6,'--','Color',[0  1 1.0000],'Linewidth',1.5)
plot(th,-pice_v(:,end)/1e6,'--','Color',[0.5 0.5 1.0000],'Linewidth',1.5)
xlabel('time [h]','interpreter','latex','FontSize',15)
ylabel ('Osmotic/cryostatic pressure [MPa]','interpreter','latex','FontSize',15)
box on
my_legend = {'Near pith $p_{osm}^{vac}$','Near bark $p_{osm}^{vac}$','$p_{osm}^{bark\ cell}$',...
    'Near pith $-p_{ice}^v$','Near bark $-p_{ice}^v$'} ;
%my_legend = {'Near pith VACs','Near bark VACs','Bark cells'} ;
legend(my_legend,'interpreter','latex','location','best','Position',[0.59 0.4 0.3 0.3])
ax = gca;
ax.FontSize = 16; 
legend boxoff
ax.TickLabelInterpreter = 'latex' ;
text(0.025,0.95,'d)','Units','normalized','FontSize',14,'interpreter','latex')
%print -depsc fig3-osm_cryo.eps

figure(6)
set (gcf,'color','white')
plot(th,Cs_p(:,1))
hold on
plot(th,Cs_p(:,end))
plot(th,Cs_bark,'k')
plot ([th(1) th(end)], [5425 5425],'r:')
plot ([th(1) th(end)], [7152 7152],'b:')
xlabel('time [h]','interpreter','latex','FontSize',15)
ylabel ('Sugar concentration [mol/m3]','interpreter','latex','FontSize',15)
box on
my_legend = {'Near pith VACs','Near bark VACs','Bark cells','Saturation limit ($0^\circ C$)','Supersaturation limit ($0^\circ C$)'} ;
legend(my_legend,'interpreter','latex','location','best')
ax = gca;
ax.FontSize = 16; 
ax.TickLabelInterpreter = 'latex' ;
legend boxoff

figure(7)
set (gcf,'color','white')
hold on
plot (th,(D-D(1))*10^+6,'Color',[1.0000  0 1.0000],'Linewidth',1.5)
plot (th,(DDx-DDx(1))*10^+6,'Color',[0 1.0 1.0000],'Linewidth',1.5)
plot ([th(1) th(end)], [Delta_D_am Delta_D_am]*1000,'k--')
xlabel('time after freezing [h]','interpreter','latex','FontSize',15)
ylabel ('Stem diameter changes [$\mu$m]','interpreter','latex','FontSize',15)
box on
% xlim ([0 9])
my_legend = {'$D_{stem}$','$D_{xylem}$','$D_{stem}$ experiment'} ;
legend(my_legend,'interpreter','latex','position',[0.5, 0.5, 0.3, 0.2])
ax = gca;
ax.FontSize = 16; 
ax.TickLabelInterpreter = 'latex' ;
legend boxoff
text(0.085,0.95,'a)','Units','normalized','FontSize',14,'interpreter','latex')
%print -depsc fig4-diam.eps

figure(8)
set (gcf,'color','white')
hold on
plot(th,Text-273.15,'k--')
plot(th,Tm_par(:,1)-273.15,'Color',[0  1 1.0000],'Linewidth',1.5)
plot(th,Tm_par(:,end)-273.15,'Color',[0.5 0.5 1.0000],'Linewidth',1.5)
plot(th,Tm_bark-273.15,'Color',[1 0 1.0000],'Linewidth',1.5)
xlabel('time [h]','interpreter','latex','FontSize',15)
ylabel ('Temperature [$\circ C$]','interpreter','latex','FontSize',15)
box on
my_legend = {'$T_{ext}$','Near pith $T_{m}^{vac}$','Near bark $T_{m}^{vac}$','$T_{m}^{bark\ cell}$'} ;
legend(my_legend,'interpreter','latex','position',[0.4, 0.65, 0.3, 0.25])
ax = gca;
ax.FontSize = 16; 
ax.TickLabelInterpreter = 'latex' ;
legend boxoff
text(0.025,0.95,'e)','Units','normalized','FontSize',14,'interpreter','latex')
%print -depsc fig3-Tm.eps

figure(19)
set (gcf,'color','white')
hold on
plot (th-5.117,(D-D(1))*10^+6,'Color',[1.0000  0 1.0000],'Linewidth',1.5)
plot (th-5.117,(DDx-DDx(1))*10^+6,'Color',[0 0 1],'Linewidth',1.5)
plot (time_expe, d_expe, 'k--')
xlabel('time after freezing [h]','interpreter','latex','FontSize',15)
ylabel ('Stem diameter changes [$\mu$m]','interpreter','latex','FontSize',15)
my_legend = {'$D_{stem}$','$D_{xylem}$','$D_{stem}$ experiment'} ;
legend(my_legend,'interpreter','latex','position',[0.55, 0.4, 0.3, 0.25])
legend boxoff
text(0.065,0.95,'a)','Units','normalized','FontSize',14,'interpreter','latex')
xlim ([0 9])
ax = gca;
ax.FontSize = 16; 
ax.TickLabelInterpreter = 'latex' ;
box on
%print -depsc fig4-diam.eps

%% Gather data to save and print them in a text file
%to_save= [th  p_save D DDx Tm_bark'-273.15 Text-273.15 Qbx] ;

% fid = fopen('data.txt','w')   ;
%     fprintf(fid,'%12s %12s %12s\n','th','pv_mean','D','Dx','Tmb','Text','Qbx');
%     fprintf(fid,'%f %f %f %f %f %f \n', to_save') ;
% fclose(fid) ;

save('pwv.mat','pwv')
save('Cs_v.mat','Cs_v')
save('Cs_p.mat','Cs_p')
save('Cs_bark.mat','Cs_bark')
save('Temp_save.mat','Temp_save')
save('p_mean','p_save')
save('D.mat','D')
save('Dx.mat','DDx')
save('Tm_bark','Tm_bark')
save('Qbx.mat','Qbx')
save('th.mat','th')
save('ru.mat','ru')

cd(olddr)
