clear all
close all

% % % This is the main program of a model developped to study winter 
% % % stem pressure and diameter variations during freeze-thaw cycles. 
% % % 
% % % The main.m program computes initial values, launches the temporal 
% % % integration with given ode15s solver options and saves data at 
% % % regular time intervals. This is the program to be launched to run the model.
% % % 
% % % The parameters.m program gather all parameters and useful 
% % % functions used to run the model
% % % 
% % % The dyfun.m program contains all physical equations and 
% % % computes the vector of unknow to be solved.
% % % 
% % % The post_process.m program is used to recompute different 
% % % state variables not saved during the main integration and exploit the results.

%% MAIN CODE
main_start = tic ; % main timer start

p = parameters ; % call the parameters.m code to gather ALL parameters and functions in the structure "p"

tf = 1*24 ; % [h] final time

%% Set initial values
T0 = zeros(p.nc,1) + p.Tinit ; % initial temperature

[Tm, Tirv, Twrv, Hirv, Hwrv, Tirf, Twrf, Hirf, Hwrf] = p.RegParam(p.Cs_v0*ones(p.nc,1)) ; % initial phase change parameters
H0 = p.TtoH(T0, Tirv, Twrv, Hirv, Hwrv, Tirf, Twrf, Hirf, Hwrf) ; % Initial enthalpy
rf0 = zeros(p.nv,1) ; % Initial fiber gas radius [m]
rv0 = p.rgv0*ones(p.nv,1) ; % Initial vessel gas radius [m]
Ufv0 = zeros(p.nv,1) ; % Initial fiber-vessel exchange [m^3]
Uroot0 = zeros(p.nv,1) ; % Initial vessel-root exchange [m^3]
Upv0 = zeros(p.nv,1) ; % Inital parenchyma-vessel exchance [m^3]
Upp0 = zeros(p.nv,1) ; % Initial parenchyma-parenchyma exchange [m3]
Vpa0 = p.Vp0*ones(p.nv,1) ; % Initial parenchyma volume [m^3]
ppa0 = p.pp0*ones(p.nv,1) ; % Initial parenchyma pressure [m^3]
nsv0 = p.Cs_v0*p.Vlv0*ones(p.nv,1)  ; % Initial sugar quantity in vessels [mol]
nslc0 = p.Cs_p0*(p.Vp0-p.Vp_bound)*ones(p.nv,1) ; % Initial sugar quantity in living cells [mol]
ngv0 = p.ngv0*ones(p.nv,1) ; % Initial gas quantity in vessels
%%% p.Vbark_cell0 will be the initial bark cell volume [m^3]
%%% p.pp0 will be the initial bark pressure [Pa]
nsb0 = p.Cs_p0*(p.Vbark_cell0-p.Vbark_bound) ; % Initial bark sugar quantity [mol]

% Initial unknown vector
y0 = [H0 ; rf0 ; rv0 ; Ufv0 ; Uroot0 ; Upv0 ; Upp0 ; Vpa0 ; ppa0 ; nsv0 ; nslc0 ; ngv0 ; p.Vbark_cell0 ; p.pp0 ; nsb0] ; 
t_stops = 0:4:tf ; % vector of instants at which you want to halt the integration
non_negative_vec = [(p.nc+8*p.nv+1):(p.nc+11*p.nv)  (p.nc+11*p.nv+3)] ; % set positions in y vector where you want to force non-negativity
options = odeset('RelTol',1e-13, 'AbsTol',1e-18,'Maxstep',1,'NormControl','on','NonNegative',non_negative_vec) ; %list of solver options

solt = [0] ; % Initial total time vector
solY = y0' ; % Initial total solution matrix

%% Main loop

for ind_loop=1:(length(t_stops)-1) % Main temporal loop 
    % with halts to chek if everything is going well
    loop_start = tic ; % loop timer start
    timespan=[t_stops(ind_loop)*3600:10:t_stops(ind_loop+1)*3600] ; % timespan between halts
    %%%timespan=[t_stops(ind_loop)*3600 t_stops(ind_loop+1)*3600] ; % timespan between halts
    [t, y] = ode15s(@(t,y) dyfun(t,y,p), timespan, y0, options) ; % main call to ode solver (over timespan, starting from y0)
    
    solt = [solt ; t(2:end)] ; % store total time
    solY = [solY ; y(2:end,:)] ; % store complete solution
    y0 = y(end,:) ; % new restart solution
    
    fprintf('Halt at t = %d hours\n', t_stops(ind_loop+1))
    loop_time = toc(loop_start)
    
    %% Save current data
    H = solY(:,1:p.nc) ;
    rf = solY(:,p.nc+1:(p.nc+p.nv)) ;
    rv = solY(:,(p.nc+p.nv+1):(p.nc+2*p.nv)) ;
    Ufv = solY(:,(p.nc+2*p.nv+1):(p.nc+3*p.nv)) ;
    Uroot = solY(:,(p.nc+3*p.nv+1):(p.nc+4*p.nv)) ;
    Upv = solY(:,(p.nc+4*p.nv+1):(p.nc+5*p.nv)) ; % extract parenchyma-vessel exchange
    Upp = solY(:,(p.nc+5*p.nv+1):(p.nc+6*p.nv)) ; % extract parenchyma-parecnhyma exchange
    Vp = solY(:,(p.nc+6*p.nv+1):(p.nc+7*p.nv)) ; % extract parenchyma exchange
    pp = solY(:,(p.nc+7*p.nv+1):(p.nc+8*p.nv)) ; % extract parenchyma pressure
    nsv = solY(:,(p.nc+8*p.nv+1):(p.nc+9*p.nv)) ; % Sugar content in vessels
    nslc = solY(:,(p.nc+9*p.nv+1):(p.nc+10*p.nv)) ; % Sugar content in vessels
    ngv = solY(:,(p.nc+10*p.nv+1):(p.nc+11*p.nv)) ; % Vessel gas content
    Vbark_cell = solY(:,(p.nc+11*p.nv+1)) ; % extract bark volume
    pbark = solY(:,(p.nc+11*p.nv+2)) ; % extract bark pressure
    nsb = solY(:,(p.nc+11*p.nv+3)) ; % extract bark sugar content
    t=solt ;

    save('t.mat','t')
    save('H.mat','H')
    save('rf.mat','rf')
    save('rv.mat','rv')
    save('Ufv.mat','Ufv')
    save('Uroot.mat','Uroot')
    save('Upv.mat','Upv')
    save('Upp.mat','Upp')
    save('Vp.mat','Vp')
    save('pp.mat','pp')
    save('nsv.mat','nsv')
    save('nslc.mat','nslc')
    save('ngv.mat','ngv')
    save('Vbark_cell.mat','Vbark_cell')
    save('pbark.mat','pbark')
    save('nsb.mat','nsb')
end

%% Save final data
    H = solY(:,1:p.nc) ;
    rf = solY(:,p.nc+1:(p.nc+p.nv)) ;
    rv = solY(:,(p.nc+p.nv+1):(p.nc+2*p.nv)) ;
    Ufv = solY(:,(p.nc+2*p.nv+1):(p.nc+3*p.nv)) ;
    Uroot = solY(:,(p.nc+3*p.nv+1):(p.nc+4*p.nv)) ;
    Upv = solY(:,(p.nc+4*p.nv+1):(p.nc+5*p.nv)) ; % extract parenchyma-vessel exchange
    Upp = solY(:,(p.nc+5*p.nv+1):(p.nc+6*p.nv)) ; % extract parenchyma-parecnhyma exchange
    Vp = solY(:,(p.nc+6*p.nv+1):(p.nc+7*p.nv)) ; % extract parenchyma exchange
    pp = solY(:,(p.nc+7*p.nv+1):(p.nc+8*p.nv)) ; % extract parenchyma pressure
    nsv = solY(:,(p.nc+8*p.nv+1):(p.nc+9*p.nv)) ; % Sugar content in vessels
    nslc = solY(:,(p.nc+9*p.nv+1):(p.nc+10*p.nv)) ; % Sugar content in vessels
    ngv = solY(:,(p.nc+10*p.nv+1):(p.nc+11*p.nv)) ; % Vessel gas content
    Vbark_cell = solY(:,(p.nc+11*p.nv+1)) ; % extract bark volume
    pbark = solY(:,(p.nc+11*p.nv+2)) ; % extract bark pressure
    nsb = solY(:,(p.nc+11*p.nv+3)) ; % extract bark sugar content
    t=solt ;

    save('t.mat','t')
    save('H.mat','H')
    save('rf.mat','rf')
    save('rv.mat','rv')
    save('Ufv.mat','Ufv')
    save('Uroot.mat','Uroot')
    save('Upv.mat','Upv')
    save('Upp.mat','Upp')
    save('Vp.mat','Vp')
    save('pp.mat','pp')
    save('nsv.mat','nsv')
    save('nslc.mat','nslc')
    save('ngv.mat','ngv')
    save('Vbark_cell.mat','Vbark_cell')
    save('pbark.mat','pbark')
    save('nsb.mat','nsb')

main_time = toc(main_start)
