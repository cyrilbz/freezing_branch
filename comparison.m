clear all 
close all

%% Program to compare data between cases
%% post_process must have been done first on each case !

root_dr = pwd ; 

 mylist =28:1:36;
% 
% for i=1:length(mylist)
%     dir_names{i} = strcat('\physio_study\PS',num2str(mylist(i))) ;
%     my_legend{i} = strcat('PS',num2str(mylist(i))) ;
% end

dir_names = {'\base_case_2' , ...
    '\base_case_2-Cs_1200', ...
     '\base_case_2-p_col_minus5' };

my_legend = {'base case 2','Cs 1200','$p_{col}=-5MPa$'} ;

%%%% loop over the different directories to load all the data you need
%%%% using 'name.mat' 
%%%% native variables available : t,H,rf,rv,uroot,ufv,upv,vp,pp,nsv,nslc,
%%%% ngv,Vbark_cell,pbark,nsb 
%%%% processed variables :
%%%% th,pwv,Cs_v,Cs_p,Cs_bark,Temp_save,p_mean,D,Dx,Tm_bark,Qbx,ru

for i=1:length(dir_names)
mydr = strcat(root_dr,dir_names{i}) ;
cd(mydr)
th(i) = load('th.mat') ;
ru(i) = load('ru.mat') ;
p_mean(i) = load('p_mean.mat') ;
Cs_bark(i) = load('Cs_bark.mat') ;
pwv(i) = load('pwv.mat') ;
pp(i) = load('pp.mat') ;
pbark(i) = load('pbark.mat') ;
D(i) = load('D.mat') ;
Dx(i) = load('Dx.mat') ;
Qbx(i) = load('Qbx.mat') ;
Tm_bark(i) = load('Tm_bark.mat') ;
Vbark_cell(i) = load('Vbark_cell.mat') ;
end

cd(root_dr)
%%% define some values you want to export
dd=zeros(length(dir_names),1) ;
pmean=zeros(length(dir_names),1) ;
Tmb=zeros(length(dir_names),1) ;

 c=cool(length(dir_names)) ; % color map
 c= [c ; c ] ; % color map
 my_styles= ["-","-","-","--","--","--"] ;
 
 my_styles=[my_styles; my_styles] ;
 
%%% plot for loop and additional processing
for i=1:length(dir_names)
    %%% instant at which you want to export some values
    ind = find(th(i).th==24) ;

    %% diameter comparison
    set (gcf,'color','white')
    buff = D(i).D ;
    buff = (buff - buff(1))*10^6 ; % diameter variations in µm
    
    figure(1) 
    hold on
    f1 = plot(th(i).th,buff,'Color',c(i,:),'Linewidth',1.5,'Linestyle',my_styles(i)) ;
    xlabel('time [h]','interpreter','latex','FontSize',15)
    ylabel('Stem diameter variations [$\mu$ m]','interpreter','latex','FontSize',15)
    %xlim([0 14])
    box on
    hold off
    dd(i) = buff(ind) ; % store diameter shrinkage
    
    %% mean pressure comparison
    buff = p_mean(i).p_save ;
    figure(2)
    set (gcf,'color','white')
    hold on
    f2 = plot(th(i).th,buff,'Color',c(i,:),'Linewidth',1.5,'Linestyle',my_styles(i)) ;
    xlabel('time [h]','interpreter','latex','FontSize',15)
    ylabel('Mean vessel pressure [kPa]','interpreter','latex','FontSize',15)
    box on
    hold off
    pmean(i)=buff(ind);
    
    %% radial pressure gradient comparison
    buff = pwv(i).pwv ;
    set (gcf,'color','white')
    figure(3)
    hold on
    semilogy(ru(i).ru/1000,buff(6470,:)/1000,'Color',c(i,:),'Linewidth',1.5,'Linestyle',my_styles(i)) ;
    set(gca, 'YScale', 'log')
    set (gcf,'color','white')
    xlabel('r [mm]','interpreter','latex','FontSize',15)
    ylabel('Vessel water pressure [kPa]','interpreter','latex','FontSize',15)
    box on
    hold off
    
%     %% bark-xylem flow rate
%     buff = Qbx(i).Qbx ;
%     figure(5)
%     set (gcf,'color','white')
%     hold on
%     f2 = plot(th(i).th,buff,'Color',c(i,:),'Linewidth',1.5,'Linestyle',my_styles(i)) ;
%     %set(gca, 'YScale', 'log')
%     xlabel('$t [h]$','interpreter','latex','FontSize',15)
%     ylabel('$Q_{bx} [m3/s]$','interpreter','latex','FontSize',15)
%    % xlim([0 12])
%     box on
%     hold off
    
    %% bark melting temperature (+ plot Text below this loop)
    buff = Tm_bark(i).Tm_bark-273.15 ;
    figure(4)
    set (gcf,'color','white')
    hold on
    f2 = plot(th(i).th,buff,'Color',c(i,:),'Linewidth',1.5,'Linestyle',my_styles(i)) ;
    xlabel('time [h]','interpreter','latex','FontSize',15)
    ylabel('$T_{m}^{bark\ cell}$ [$^\circ $C]','interpreter','latex','FontSize',15)
    %xlim([5 8])
    box on
    hold off
    set (gcf,'color','white')
    Tmb(i)=buff(ind);
    
%     %% XYLEM diameter comparison
%     buff = Dx(i).DDx ;
%     buff = (buff - buff(1))*10^6 ; % diameter variations in µm
%     figure(6) 
%     hold on
%     f1 = plot(th(i).th,buff,'Color',c(i,:),'Linewidth',1.5,'Linestyle',my_styles(i)) ;
%     xlabel('time [h]','interpreter','latex','FontSize',15)
%     ylabel('Xylem diameter variations [$\mu$ m]','interpreter','latex','FontSize',15)
%     box on
%     hold off
    
%     %% Turgor pressure comparison
%     buff = pp(i).pp ;
%     figure(7)
%     set (gcf,'color','white')
%     hold on
%     f2 = plot(th(i).th,buff(:,end),'Color',c(i,:),'Linewidth',1.5,'Linestyle',my_styles(i)) ;
%     xlabel('time [h]','interpreter','latex','FontSize',15)
%     ylabel('Near bark turgor pressure [kPa]','interpreter','latex','FontSize',15)
%     box on
%     hold off
    
%     %% Turgor pressure comparison
%     buff = pbark(i).pbark ;
%     figure(8)
%     set (gcf,'color','white')
%     hold on
%     f2 = plot(th(i).th,buff/1e6,'Color',c(i,:),'Linewidth',1.5,'Linestyle',my_styles(i)) ;
%     xlabel('time [h]','interpreter','latex','FontSize',15)
%     ylabel('Bark turgor pressure [MPa]','interpreter','latex','FontSize',15)
%     box on
%     hold off
    
%     %% Turgor pressure comparison
%     buff = Vbark_cell(i).Vbark_cell ;
%     figure(9)
%     set (gcf,'color','white')
%     hold on
%     f2 = plot(th(i).th,buff(:)/buff(1),'Color',c(i,:),'Linewidth',1.5,'Linestyle',my_styles(i)) ;
%     xlabel('time [h]','interpreter','latex','FontSize',15)
%     ylabel('Bark cell RELATIVE volume [m3]','interpreter','latex','FontSize',15)
%     box on
%     hold off
end
% save('Tmb_24h.mat','Tmb')
% save('pmean_24h.mat','pmean')
% save('dd_24h.mat','dd')

%% add legends
for i=1:4
    
    figure(i)
    lgnd(i) = legend(my_legend,'interpreter','latex','location','north');
    ax = gca;
    ax.FontSize = 16;       
    ax.TickLabelInterpreter = 'latex' ;
    %lgnd(i).Title.String = {'$k_{ray}=3.63\times 10^{-n}$', '[m$^2$]'} ;
    legend boxoff
    %set(ax,'Color','parula')
end

%% load experimental diameter variations

data_expe = importdata('data_expe.txt') ;

time_expe = data_expe.data(:,3) ; 
d_expe = data_expe.data(:,2) ; 

%% additional plots
Text = max(5 -th(1).th,-10) ; 
Text2 = max(5 -8*th(1).th,-10) ; 
figure(4)
hold on
plot(th(1).th,Text,'k:','DisplayName','$T_{ext} (-1K/h)$')
%plot(th(1).th,Text2,'k-.','DisplayName','$T_{ext} (-8K/h)$')
set(lgnd(4),'color','none')
text(0.075,0.95,'a)','Units','normalized','FontSize',14,'interpreter','latex')
hold off
%print -depsc fig6-Tm_dTCs.eps

figure(2)
text(0.025,0.95,'e)','Units','normalized','FontSize',14,'interpreter','latex')

% h = axes('Parent',gcf,'Position',[0.45 0.25 0.38 0.38] ) ;
% for i=1:length(dir_names)
%     buff = p_mean(i).p_save ;
%     hold on
%     f2 = plot(th(i).th,buff,'Color',c(i,:),'Linewidth',1.5) ;
% end
% box on
% %plot(h,[5 10], [-50 -100]);
% ax_sub = gca;
% ax_sub.TickLabelInterpreter = 'latex' ;
% ax_sub.FontSize = 16; 
% set(h, 'Xlim', [5 5.75], 'Ylim', [15 30]);
%print -depsc fig6-pwv_dTCs.eps

 figure(1)
% %hold on
 text(0.045,0.95,'d)','Units','normalized','FontSize',14,'interpreter','latex')
 ax= gca;
 %print -depsc figS1-diam_kray.eps
 
  figure(3)
% %hold on
 text(0.045,0.95,'f)','Units','normalized','FontSize',14,'interpreter','latex')
 ax= gca;
%print -depsc figS1-pwv_r_kray.eps
%legend boxoff
 
figure(2)
set(lgnd(2),'Position',[0.46, 0.48, 0.3, 0.3])
%print -depsc figS1-mean_pwv_kray.eps
%  %ax.LineStyleCyclingMethod ='withcolor' ;
% h = axes('Parent',gcf,'Position',[0.45 0.45 0.38 0.38] ) ;
% for i=1:length(dir_names)
%     buff = D(i).D ;
%     buff = (buff - buff(1))*10^6 ; % diameter variations in µm
%     hold on
%     plot(h, th(i).th,buff,'Color',c(i,:),'Linewidth',1.5) ;
% end
% box on
% %plot(h,[5 10], [-50 -100]);
% ax_sub = gca;
% ax_sub.TickLabelInterpreter = 'latex' ;
% ax_sub.FontSize = 16; 
% set(h, 'Xlim', [5 5.75], 'Ylim', [-100 0]);

%print -depsc fig6-diam_dTCs.eps

% set(lgnd(1),'color','none')
% hold on
% plot (time_expe, d_expe, 'k--','DisplayName','Experiment')
% hold on
