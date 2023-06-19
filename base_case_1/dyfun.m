% % This file is part of "Freezing Branch", a code that simulates 
% % freeze induced cell dehydration, pressure and diameter changes 
% % in a tree stem.
% % 
% % Author: Cyril Bozonnet (cyril.bozonnet@inrae.fr) 
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


%% Main function to compute dy/dt
function dydt=dyfun(t,y,p) 

%% Extract all datas from solution vector

    H = y(1:p.nc) ; % Enthalpy
    rf = y(p.nc+1:(p.nc+p.nv)) ; % fiber bubble radius
    rv = y((p.nc+p.nv+1):(p.nc+2*p.nv)) ; % vessel bubble radius
    Ufv = y((p.nc+2*p.nv+1):(p.nc+3*p.nv)) ; % Fiber fessel flow
    Uroot = y((p.nc+3*p.nv+1):(p.nc+4*p.nv)) ; % Root vessel flow
    Upv = y((p.nc+4*p.nv+1):(p.nc+5*p.nv)) ; % Parenchyma-vessel flow
    Upp = y((p.nc+5*p.nv+1):(p.nc+6*p.nv)) ; % Inter-parecnhyma flow 
    Vp = y((p.nc+6*p.nv+1):(p.nc+7*p.nv)) ; % Parenchyma volume
    pp = y((p.nc+7*p.nv+1):(p.nc+8*p.nv)) ; % Parenchyma pressure
    nsv = y((p.nc+8*p.nv+1):(p.nc+9*p.nv)) ; % Vessel sugar content
    nslc = y((p.nc+9*p.nv+1):(p.nc+10*p.nv)) ; % Living cell sugar content
    ngv = y((p.nc+10*p.nv+1):(p.nc+11*p.nv)) ; % Vessel gas content
    Vbark_cell = y((p.nc+11*p.nv+1)) ; % bark cell volume
    pbark = y((p.nc+11*p.nv+2)) ; % bark pressure
    nsb = y((p.nc+11*p.nv+3)) ; % bark sugar content
    
%% Initialize time derivatives
    dHdt = zeros(p.nc,1) ; 
    drfdt = zeros(p.nv,1) ;
    drvdt = zeros(p.nv,1) ;
    dUfv = zeros(p.nv,1) ;
    dUroot = zeros(p.nv,1) ;
    dUpv = zeros(p.nv,1) ;
    dVp = zeros(p.nv,1) ;
    dppdt = zeros(p.nv,1) ;
    dnsvdt = zeros(p.nv,1) ;
    dnslcdt = zeros(p.nv,1) ;
    dngvdt = zeros(p.nv,1) ; 
    dVbark_cell = 0 ;
    dpbark = 0 ;
    dnsbdt = 0 ;
    
%% Macroscale phase change

    % get microscale sugar concentration (C=n/V) ; V being the volume of water in a vessel
        Cs_v = nsv./(pi*(p.Rv^2-rv.^2)*p.Lv) ; % This is on the vessel grid !
    % Upscale it (interpolate on the fine grid + fill end points using nearest value
        Cs_v_nc = interp1(p.ru,Cs_v,p.r,'linear') ;
        Cs_v_nc = fillmissing(Cs_v_nc,'linear','EndValues', 'nearest') ;
    
    % Compute useful quantities ...
    % ...(physical properties, initial and boundary temperature/enthalpy values)
        [Tm, Tirv, Twrv, Hirv, Hwrv, Tirf, Twrf, Hirf, Hwrf] = p.RegParam(Cs_v_nc') ; % get regularized parameters
        Temp = p.HtoT(H, Tirv, Twrv, Hirv, Hwrv, Tirf, Twrf, Hirf, Hwrf) ; % converts enthalpy to temperature
        RHO = p.HtoRHO(H, Hirf, Hwrf, Hirv, Hwrv) ; % converts enthalpy to density
        Text = p.Tempout(t) ; % Compute external temperature
        Hext = p.TtoH(Text, Tirv(end), Twrv(end), Hirv(end), Hwrv(end), Tirf(end), Twrf(end), Hirf(end), Hwrf(end)) ; % Imposed enthalpy on right bc   
        k = p.Difftensor([H ; Hext], [Hirv ; Hirv(end)], [Hwrv ; Hwrv(end)], [Hirf ; Hirf(end)], [Hwrf ; Hwrf(end)]) ; % get diffusivity on cell faces - add boundary point
    
    % Heat equation here below
        flux = p.mean*(k(1:end-1)).*(p.grad_up*Temp) ; % Heat flux leaving each cell from the right (k_face*dT/dr)
        dHdt = p.div*flux +k(1:end-1).*p.grad_c*Temp; % dH/dt = div.(kgradT) + k/r*dT/dr

    % Boundary conditions
        bc_flux = 1/2*(k(end)+k(end-1))*(Text-Temp(p.nc))/p.dr ; % External boundary heat flux
        dHdt(p.nc,1) = k(end-1)/p.r(p.nc)*(Text-Temp(p.nc-1))/2/p.dr ... % Boundary condition : k/r*dT/dr ...
                        + (bc_flux-flux(end-1))/p.dr ; % ... + div(k*gradT)
    
    dHdt = (1./RHO).*dHdt ; % normalize the vector by the density

%% Microscale modelisation
    % Get temperature and enthalpy on a vector of size Nv (Nc = MR * Nv)
        Temp_nv = Temp(1:p.MR:end) ; % Temperature on the coarser grid
        H_nv = H(1:p.MR:end) ; % Enthalpy on the coarser grid
        Hwrf_nv = Hwrf(1:p.MR:end) ; % Fiber water enthalpy on the coarse grid
        Hwrv_nv = Hwrv(1:p.MR:end) ; % Vessel water enthalpy on the coarse grid
        Tm_nv = Tm(1:p.MR:end) ; % Melting temperature on the coarse grid
    
    % Compute vessel volumetric freezing indicators (0: thawed - 1: frozen) 
        fv = min(1,max((Hwrv_nv-H_nv)/p.DHv,0)) ;
    
    % Compute gas pressures (P=nRT/V)
        pgv = ngv*p.Rg.*Temp_nv./(pi*rv.^2*p.Lv) ; % Vessel gas pressure
    
    % Compute water pressure (Laplace's law)
        pwv = pgv - p.saw./rv ; % Vessel water pressure
    
    % Compute sugar related quantities (Cs=n/(V-V_bound))
        Cs_p = nslc./(Vp-p.Vp_bound) ;
        Cs_bark = nsb/(Vbark_cell-p.Vbark_bound) ;
    
    % Compute viscosity (corrective factor that must multiply the original viscosity)
        mu_r=1/p.muw*0.09607e-3*exp(2.9*(p.Tc./Temp_nv).^3).*(1+p.e_nu*Cs_p/p.rhow.*exp(((Cs_p/p.rhow).^p.f_nu)./(p.g_nu*Temp_nv/p.Tc+p.h_nu))) ;
        mu_bark_r = 1/p.muw*0.09607e-3*exp(2.9*(p.Tc./Text).^3).*(1+p.e_nu*Cs_bark/p.rhow.*exp(((Cs_bark/p.rhow).^p.f_nu)./(p.g_nu*Text/p.Tc+p.h_nu))) ;
        mean_mu = [p.mean_nv*mu_r] ;% averaged viscosity
        mean_mu(end) = 0.5*(mu_r(end) + mu_bark_r) ; % boundary condition
    
    % Compute osmotic pressure
        posm_pv = 1*(Cs_p-Cs_v)*p.Rg.*Temp_nv ; % Osmotic pressure difference between parenchyma and vessels
   
    % Compute cryostatic suction
        pice_v = (- 1*p.rhow*(p.Ew-p.Ei).*log(Temp_nv./p.Tc)).*fv ;
            
    % Compute vessel-parenchyma flux (rescaled by Nstack as it is the flux towards one parenchyma)
        dUpv = -p.PC./mu_r*p.Apv/p.Nstack.*(pwv-pice_v-pp+posm_pv) ; % Darcy's law
        dUpv(dUpv<0) = (1-fv(dUpv<0)).*dUpv(dUpv<0) ; % Block ice backflow
    
    % Compute parenchyma ray flux
        dUpp = p.RayCond./mean_mu.*diff([pp ; pbark]) - p.RayCond./mean_mu.*diff([Cs_p*p.Rg.*Temp_nv ;Cs_bark*p.Rg*Text]); % Darcy's law - positive from right to left
    
    % Compute parenchyma & bark volume variation
        dVp = -dUpv + diff([0 ; dUpp]) ; %  Volume variation = Sum of in & out flow rates
        dVbark = -p.Nray*p.Nstack*dUpp(p.nv) ; % Volume variation for the whole bark
        dVbark_cell = dVbark/p.Nbark_cells ; % Volume variation of an individual bark cell

    % Compute vessel radius variations 
        drvdt = -(p.Nstack*dUpv)./rv/(2*pi*p.Lv) ; % vessel gas bubble variations
    
    % Compute elastic modulus
        myB_par = p.Bpar-(p.RWC_plus-(Vp)/(p.Vp0))./(p.RWC_plus-p.RWC_minus)*(p.Bpar - p.Bmpar) ;
        myB_par= max(p.Bmpar,min(p.Bpar,myB_par)) ;
        myB_bark = p.Bbark-(p.RWC_plus_bark-(Vbark_cell)/(p.Vbark_cell0))/(p.RWC_plus_bark-p.RWC_minus_bark)*(p.Bbark - p.Bmbark) ;
        myB_bark= max(p.Bmbark,min(p.Bbark,myB_bark)) ;
    
    % Compute parenchyma turgor pressure variations
        dppdt = myB_par.*dVp./(Vp) ; % Hooke's law for parenchyma cell elasticity
        dpbark = myB_bark/(Vbark_cell)*dVbark_cell ; % Hooke's law for bark elasticity
    
    %% Gather all derivatives before resolution

    dydt = [dHdt ; drfdt ; drvdt ; dUfv ; dUroot ; dUpv ; dUpp ; dVp ; dppdt ; dnsvdt ; dnslcdt ; dngvdt ; dVbark_cell ; dpbark ; dnsbdt] ; 

end
