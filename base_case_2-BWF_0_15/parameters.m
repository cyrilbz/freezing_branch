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


function p = parameters(p) % all parameters are wrapped in a function to pass them to ode15s easily
%% Base case parameters
p.Rtree = 0.0075; % Tree radius [m]
p.reflect = 0.0 ; %root reflection coefficient
p.Rheart= 0.003 ; % 6mm pith radius, constant with diameter - no duramen formation for small enough branches !
p.bt = -2.5513*p.Rtree^2 + 0.2367*p.Rtree - 0.0004 ; % Bark layer thickness [m] - FIT OF EXPERIMENTAL RESULTS
p.L= p.Rtree - p.Rheart - p.bt ; % Sapwood fraction -> constant with diameter

lvd = 2.5e+3 ; % Linear vessel density [nb vessels/m/ray]
p.nv = round(p.L*lvd) ; % Number of vessels/ray
p.MR = 3; % Integer larger than 1 ! Mesh ratio between number of vessels and cells used to discretize the thermal problem
p.nc = p.nv*p.MR ; % number of cells for discretization [-]

% compute here the space vector as well as the discretized radius vector
p.dr = p.L/(p.nc +1) ; % space step
p.dru = p.L/(p.nv +1) ; % space step between functional units
p.r = (p.Rheart + p.dr) : p.dr : (p.Rtree-p.dr-p.bt) ; % space vector (NB: r={Rheart ; Rtree-bt} are boundary "ghost" points)
p.ru = (p.Rheart + p.dru) : p.dru : (p.Rtree-p.dru-p.bt) ;

p.Tc = 273.15 ; % Normal critical temperature for fiber freezing (=0degC)
p.Toutmax = p.Tc + 5 ;   % Maximal outside temperature   
p.Toutmin = p.Tc -10 ;  % Minimal outside temperature
p.Tinit   = p.Tc + 5 ; % Initial temperature 
p.Rc = 1 ; % [K/h] : cooling rate
p.Tmean  = (p.Toutmax + p.Toutmin) / 2 ; % Mean temperature

p.memb_perm_fv= 0.15 ; % Fiber-vessel membrane permeability

p.cinf = 1e+8 ; % [J/K/kg] : regularization parameter in T-H relationship

%% Physical properties (phase change)
p.Ew = 907e+3 ; % [J/kg] : Water enthalpy at Tm
p.Ei = 574e+3 ; % [J/kg] : Ice enthalpy at Tm
p.kw = 10*0.556 ; % [W/m/K] : Water thermal conductivity
p.ki = 10*2.22 ; % [W/m/K] : Ice thermal conductivity
p.rhow = 1000 ; % [kg/m^3] : Water density
p.rhoi = 917 ; % [kg/m^3] : Ice density
p.cw = 4180 ; % [J/K/kg] : water heat capacity
p.ci = 2100 ; % [J/K/kg] : ice heat capacity

%% Physical constants (thermodynamic - fluid mechanics)
p.Ms = 0.3423 ; % Molar mass of sugar [kg/mol]
p.Mg = 0.029 ; % Molar mass of air [kg/mol]
p.Rg = 8.314 ; % Universal gas constant [J/K/mol]
p.Cs_v0 =  50 ; % p.gamma_v/(1-p.gamma_v)*p.rhow/p.Ms ; % Sugar concentration [mol/m^3]
p.Cs_p0 = 200 ; % VAC sugar concentration [mol/m3]
p.siw = 0.033 ; % Ice-water surface tension [N/m]
p.saw = 0.072 ; % Air-water surface tension [N/m]
p.Henry = 0.0 ; %0.0274 ; % Henry's constant for air in water [-]
p.Kb = 1.853 ; % Cryoscopic constant [kg/K/mol]
p.Tm = p.Tc - p.Cs_v0*p.Kb/p.rhow ; % melting temperature with sugar effect [K]

%% Anatomical parameters (microscale modelling)

% Vessel & fibers parameters
p.Rv = 8e-5 ; % Vessel radius [m]
p.Rf = 5e-6 ; % Fiber radius [m]
p.Lv = 1e-3 ; % Vessel length [m]
p.Lf = 1e-3 ; % Fiber length [m]
p.Nf = 16 ; % Number of fibers connected to a vessel
p.rcap = 7.8e-7 ; % radius of pores in fibre-vessel wall [m]
ffv = 0.88 ; % fraction of vessel wall covered by fiber-vessel pits
p.W = 3.64e-6 ; % Fiber-vessel wall thickness [m]
p.muw = 1.6e-3 ; % Water dynamic viscosity at 2°C [Pa.s]
kfv = 3.63e-21 ; % Fiber-vessel wall porosity [m2] (Petty et al 1983, JExpBot)
p.WC = 0 ; %kfv/muw/W ; % Conductivity of fibre-vessel wall [m/s/Pa] 
p.Afv = ffv*2*pi*p.Rv*p.Lv ; % area of fibre vessel-wall [m2]

% Ray parameters
p.Rray = 5e-6 ; % Ray radius [m] => Measured on picture "Noyer-CL-x10-1"
p.l_par = 20e-6 ; % Parenchyma length [m]
p.PC = 1*kfv/p.muw/p.W ; % 1e-14 ; % Pit conductivity [m/s/Pa] 
p.RayCond =10000*kfv/p.muw/p.dru*pi*p.Rray^2 ; %pi*p.Rray^4/(8*muw*p.dru) ; % Ray conductivity using Hagen-Poiseuille [m^3/s/Pa]
p.Apv = (1-ffv)*2*pi*p.Rv*p.Lv ;% npits*pi*p.rpv^2 ; % Area of the vessel-parenchyma pit [m2]
p.trd = 5000 ; % 5000 ; % tangential ray density [1/m]
p.Nray = round(p.trd*2*pi*p.Rtree) ; % Number of rays in a tree section [-]
p.Nstack =p.Apv/(2*p.Rray*p.l_par) ; % Number of rays connected to a vessel [-]
p.BW_ray = 0.1 ; % Percent of bound water
p.Vp0 = (pi*p.Rray^2)*p.l_par ; % Initial water volume
p.Vp_bound = p.BW_ray*p.Vp0 ; % Volume of bound water

% Bark parameters
p.BFF= 0.15 ; %0.75*2*p.Rtree + 0.048 ; % Bark Free water Fraction 
p.Vbark0 = p.BFF*pi*(-p.bt^2 + 2*p.bt*p.Rtree)*p.Lv ; % Initial Bark water volume [m]
p.R_bark_cell = 10e-6 ; % Bark cell radius
p.L_bark_cell = 15e-6 ; % Bark cell length
p.Vbark_cell0 = pi*p.R_bark_cell^2*p.L_bark_cell ; % Initial water volume in bark cells [m3]
p.BW_bark = 0.1 ; % Percent of bound water
p.Vbark_bound = p.BW_bark*p.Vbark_cell0 ; % Volume of bound water

% Sugar effect on viscosity
p.e_nu = 0.730 ;
p.g_nu = 8.345 ;
p.h_nu = -7.042 ;
p.f_nu = 1.10 ;

% Mechanical parameters
p.p_col_par = 0; % Pressure at which VACs collapse
p.p_col_bark = 0 ; % Pressure at which bark cells collapse
p.Bpar =10e6 ; % [Pa] Parenchyma cell wall bulk modulus
p.Bmpar = 1e3 ; % 0.005*p.Bpar ; % [Pa] parenchyma cell membrane buk modulus
p.Bbark = 10e6 ; % [Pa] Bark bulk modulus
p.Bmbark = 1e3 ; % [Pa] bark cell membrane buk modulus

%% Initialization (microscale modelling)
p.Vgv0 = 0.2 ; % Initial vessel gas volume fraction [-]
p.rgv0 = sqrt(p.Vgv0)*p.Rv ; % Initial vessel gas radius [m]
p.Vlv0 = (1-p.Vgv0)*pi*p.Rv^2*p.Lv ; % Initial vessel liquid volume [m]

p.pwv0 = 2e4 ; % p.memb_perm_fv*p.Cs_v0*p.Rg*p.Tinit ; %Initial vessel water pressure 
p.ngv0 = pi*p.rgv0^2*p.Lv/p.Rg/p.Tinit*(p.pwv0+p.saw/p.rgv0) ; % Initial vessel gas density [kg/m^3]

p.pp0 = p.pwv0 +(p.Cs_p0-p.Cs_v0)*p.Rg*p.Tinit; % Initial parenchyma pressure 
p.Vp_TLP = p.Vp0/exp((p.pp0-p.p_col_par)/p.Bpar) ; % Initial Parenchyma volume [m3] 
p.Vb_TLP = p.Vbark_cell0/exp((p.pp0-p.p_col_bark)/p.Bbark) ; % Turgor loss point in bark living cells [m3]
p.Nbark_cells = p.Vbark0/p.Vbark_cell0 ;

p.RWC_tlp = p.Vp_TLP/p.Vp0 ; % RWC at TLP in parenchyma
p.slope_b= 1e10 ; % regularization for B-V relationship
p.RWC_plus = p.RWC_tlp + (p.Bpar-p.Bmpar)/(2*p.slope_b) ;
p.RWC_minus = p.RWC_tlp - (p.Bpar-p.Bmpar)/(2*p.slope_b) ;

p.RWC_tlp_bark = p.Vb_TLP/p.Vbark_cell0 ; % RWC at TLP in bark
p.RWC_plus_bark = p.RWC_tlp_bark + (p.Bbark-p.Bmbark)/(2*p.slope_b);
p.RWC_minus_bark = p.RWC_tlp_bark - (p.Bbark-p.Bmbark)/(2*p.slope_b);

p.Db0 = 2*p.Rtree ; % Initial tree diameter 

%% Regularized parameters (two-steps phase change)
p.Ratio = 1;
p.DHv = (p.Ew-p.Ei)*p.Ratio ; % Latent heat portion attributed to vessels
p.DHf = (p.Ew-p.Ei)*(1-p.Ratio) ; % Latent heat portion attributed to fibers
p.cmixed = p.cw*p.Ratio + p.ci*(1-p.Ratio) ; % Mixed heat capacity
p.kmixed = p.kw*p.Ratio + p.ki*(1-p.Ratio) ; % Mixed conductivity
p.rhomixed = p.rhow*p.Ratio + p.rhoi*(1-p.Ratio) ; % Mixed density

%% Set up anonymous function calls   
% The actual functions are provided below.
% Allows calling functions without passing too many parameters
p.RegParam = @(x) funRegParam(x, p.Tc, p.Kb, p.rhow, p.DHv, p.DHf, p.cinf, p.ci, p.cmixed) ; % x=Csv
p.HtoT = @(x, Tirv, Twrv, Hirv, Hwrv, Tirf, Twrf, Hirf, Hwrf) funHtoT(x, p.ci, p.cw, p.cinf, p.cmixed, Tirv, Twrv, Hirv, Hwrv, Tirf, Twrf, Hirf, Hwrf); % H-->T (x=H)
p.TtoH = @(x, Tirv, Twrv, Hirv, Hwrv, Tirf, Twrf, Hirf, Hwrf) funTtoH(x, p.ci, p.cw, p.cinf, p.cmixed, Tirv, Twrv, Hirv, Hwrv, Tirf, Twrf, Hirf, Hwrf); % T-->H (x=T)
p.Difftensor = @(x, Hirv, Hwrv, Hirf, Hwrf) funDifftensor(x, p.ki, Hirv, p.kw, Hwrv, Hirf, Hwrf, p.kmixed); % H-->k (x=H)
p.HtoRHO = @(x, Hirf, Hwrf, Hirv, Hwrv) funHtoRHO(x, p.rhoi, Hirf, p.rhow, Hwrf, p.rhomixed, Hirv, Hwrv) ; % H --> RHO (x=H)
p.Tempout = @(t) Tempfun1( t, p.Tmean, p.Tinit, p.Toutmax, p.Toutmin, p.Rc); % t-->T (t: time in seconds)

%% Finite difference discretization

% Build d²/dr² operator using second order centered scheme
p.A = - diag(2/(p.dr^2)*ones(p.nc,1)) + diag(1/(p.dr^2)*ones(p.nc-1,1), 1) ...
    + diag(1/(p.dr^2)*ones(p.nc-1,1), -1) ; 
p.A(1,1) = -1/(p.dr^2) ; % Neumann BC on left point (hypothesis : Heartwood is not conductive)

% Build gradient d/dr operator df/dr_(i+1/2) = ((f(i+1)-f(i))/dr
p.grad_up = -diag(ones(p.nc,1)) + diag(ones(p.nc-1,1),1) ; 
p.grad_up = 1/p.dr.*p.grad_up ;

% Build CENTERED 1/r*d/dr operator df/dr_i = (f(i+1)-f(i-1))/(2*dr)
p.grad_c = diag(1 ./ p.r(1:p.nc-1), 1) - diag(1 ./ p.r(2:p.nc), -1);
p.grad_c(1,1) = -1/p.r(1) ; % Neumann BC on left point
p.grad_c = p.grad_c./(2*p.dr) ; 

% Build divergence d/dr operator df/dr_(i) = 1/r*((f(i+1/2)-f(i-1/2))/dr
p.div = diag(1/p.dr*ones(p.nc,1)) -diag(1/p.dr*ones(p.nc-1,1),-1) ; 

% Averaging operator (1/2(f(i)+f(i+1)) for face centered quantities
p.mean = diag(1/2*ones(p.nc,1)) + diag(1/2*ones(p.nc-1,1),1) ;

p.mean_nv = diag(1/2*ones(p.nv,1)) + diag(1/2*ones(p.nv-1,1),1) ; % mean operator on the fluid grid
end

%% Sugar concentration to regularized phase change parameters
function [Tm, Tirv, Twrv, Hirv, Hwrv, Tirf, Twrf, Hirf, Hwrf] = funRegParam(Csv, Tc, Kb, rhow, DHv, DHf, cinf, ci, cmixed)

Tirf = zeros(size(Csv)) ;
Twrf = zeros(size(Csv)) ;

Tm = Tc - Csv*Kb/rhow ; % Melting temeprature (vector)

Tirv = Tm - DHv/(2*cinf) ; % [K] : regularized Vessel ice melting temperature
Twrv = Tm + DHv/(2*cinf) ; % [K] : regularized Vessel water freezing temperature
Hirv = ci*Tirv ; % [J/kg] : regularized Vessel ice enthalpy
Hwrv = Hirv + DHv ; % [J/kg] : regularized Vessel water enthalpy

Tirf(:) = Tc - DHf/(2*cinf) ; % [K] : regularized Fiber ice melting temperature
Twrf(:) = Tc + DHf/(2*cinf) ; % [K] : regularized Fiber water freezing temperature
Hirf = Hwrv + cmixed*(Tirf-Twrv) ; % [J/kg] : regularized Fiber ice enthalpy
Hwrf = Hirf + DHf ; % [J/kg] : regularized Fiber water enthalpy
end

%% Enthalpy - Temperature relationship
% function to compute temperature from enthalpy
function T = funHtoT(H, ci, cw, cinf, cmixed, Tirv, Twrv, Hirv, Hwrv, Tirf, Twrf, Hirf, Hwrf) 
T = zeros(size(H));

if isequal(size(H),size(Tirv)) % check vector sizes compatibility 
    % Ok
else
    error('Reg Parameters must be the same size as H vector')
end

for ind = 1 : length(H)
  if H(ind) <= Hirv(ind)
    T(ind) = H(ind) / ci; % Temperature in full ice state
  elseif (Hirv(ind) < H(ind)) && (H(ind) <= Hwrv(ind))
    T(ind) = Tirv(ind) + (H(ind)-Hirv(ind)) / cinf; % Temperature during vessel freezing
  elseif (Hwrv(ind) < H(ind)) && (H(ind) <= Hirf(ind))
    T(ind) = Twrv(ind) + (H(ind)-Hwrv(ind)) / cmixed ; % Temperature in mixed state
  elseif (Hirf(ind) < H(ind)) && (H(ind) <= Hwrf(ind))
    T(ind) = Tirf(ind) + (H(ind)-Hirf(ind)) / cinf; % Temperature during fiber freezing
  else
    T(ind) = Twrf(ind) + (H(ind)-Hwrf(ind)) / cw; % Temperature in full water state
  end
end
end  % funHtoT.


%% Temperature enthalpy relationship
% function to compute enthalpy from temperature
function [H]= funTtoH(T, ci, cw, cinf, cmixed, Tirv, Twrv, Hirv, Hwrv, Tirf, Twrf, Hirf, Hwrf) 
H = zeros(size(T)) ;

for ind = 1 : length(T)
  if T(ind) <= Tirv(ind)
    H(ind) = T(ind) * ci; % Enthalpy in full ice state
  elseif (Tirv(ind < T(ind)) && (T(ind) <= Twrv(ind)))
    H(ind) = Hirv(ind) + (T(ind) - Tirv(ind)) * cinf;  % Enthalpy during vessel freezing
  elseif (Twrv(ind) < T(ind)) && (T(ind) <= Tirf(ind))
    H(ind) = Hwrv(ind) + (T(ind) - Twrv(ind)) * cmixed; % Enthalpy in mixed state
  elseif (Tirf(ind) < T(ind)) && (T(ind) <= Twrf(ind))
    H(ind) = Hirf(ind) + (T(ind) - Tirf(ind)) * cinf; % Enthalpy during fiber freezing
  else
    H(ind) = Hwrf(ind) + (T(ind)-Twrf(ind))*cw ; % Enthalpy in full water state
  end
end
end  % funTtoH.

%% Enthalpy to Diffusivity relationship
function D = funDifftensor( H, ki, Hirv, kw, Hwrv, Hirf, Hwrf, kmixed)
D = zeros(size(H));

for ind = 1 : length(H)
    if H(ind) <= Hirv(ind)
        D(ind) = ki; % Diffusivity in ful ice state
      elseif (Hirv(ind) < H(ind)) && (H(ind) <= Hwrv(ind))
        D(ind) = ki + (H(ind)-Hirv(ind)) / (Hwrv(ind)-Hirv(ind))*(kmixed-ki); % Diffusivity during vessel freezing
    elseif (Hwrv(ind) < H(ind)) && (H(ind) <= Hirf(ind))
        D(ind) = kmixed ; % Diffusivity in mixed state
    elseif (Hirf(ind) < H(ind)) && (H(ind) <= Hwrf(ind))
        D(ind) = kmixed + (H(ind)-Hirf(ind)) / (Hwrf(ind)-Hirf(ind))*(kw-kmixed); % Diffusivity during fiber freezing
    else
        D(ind) = kw ;
    end
end
end  % Difftensor.

%% Enthalpy to Density relationship
function RHO = funHtoRHO(H, rhoi, Hirf, rhow, Hwrf, rhomixed, Hirv, Hwrv)
RHO = zeros(size(H));

for ind = 1 : length(H)
    if H(ind) <= Hirv(ind)
        RHO(ind) = rhoi; % Density in ful ice state
      elseif (Hirv(ind) < H(ind)) && (H(ind) <= Hwrv(ind))
        RHO(ind) = rhoi + (H(ind)-Hirv(ind)) / (Hwrv(ind)-Hirv(ind))*(rhomixed-rhoi); % Density during vessel freezing
    elseif (Hwrv(ind) < H(ind)) && (H(ind) <= Hirf(ind))
        RHO(ind) = rhomixed ; % Density in mixed state
    elseif (Hirf(ind) < H(ind)) && (H(ind) <= Hwrf(ind))
        RHO(ind) = rhomixed + (H(ind)-Hirf(ind)) / (Hwrf(ind)-Hirf(ind))*(rhow-rhomixed); % Density during
    else
        RHO(ind) = rhow ;
    end
end
end  % HtoRHO

%% Imposed temperature
function T = Tempfun1 (t, Tmean, Tinit, Tmax, Tmin, Rc) 
    % Linear temperature variation starting from t=0s and T= Tinit
    T = max(Tinit -Rc*t/3600,Tmin) ; 
    % Sinusoidal variation
%     A = (Tmax-Tmin)/2 ; % Sinus amplitude
%     phi = asin((Tmean-Tinit)/A) ; % phase lag
%     T = Tmean -  A* sin(2*pi*t/3600/24+phi);
end
