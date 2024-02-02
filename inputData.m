function data = inputData(numberIt)
if nargin == 0
    numberIt = 1;
end

%% canister geometry and physical properties of MH

%canister geometry
base = load('C:\Users\marco\OneDrive\Documenti\Università - Dottorato\Script Matlab\BaseDesign4_2D1Dcomparison.mat');
% data.D      = base.diameter(3);                                            %canister internal diameter
data.D      = 0.0254;
data.MH_nodes    = 50;                                                     %by increasing MH_nodes, also N and r must be increased
data.deltaR = (data.D/2)/(data.MH_nodes-1);
% data.L      = base.length(3);
data.L      = data.D*50;
data.A      = pi*data.D*data.L;
data.V      = ones(data.MH_nodes,1);                                       %canister internal volume
data.V(1)   = pi*data.deltaR^2/4*data.L;
for i = 2:data.MH_nodes-1
    data.V(i) = 2*pi*data.L*data.deltaR^2*(i-1);
end
data.V(data.MH_nodes) = pi*data.L*(data.D^2-4*data.deltaR^2*(.5+(data.MH_nodes-2))^2)/4;
data.P0    = 1e5;                                                          %atmospheric reference pressure
data.fV    = 0.25;                                                         %25% free volume for expansion

%hydride physical properties
% solidMH       = 14*sum(data.V)/(pi*0.11^2/4*0.54);
data.kMH      = 2.0;                                                       %hydride conductivity
data.rhoMH    = 7260;
data.cpMH     = 419;                                                       %hydride specific heat
data.SC       = 3.0;                                                      %stoichiometric coefficient
data.MW_MH    = .432;                                                      %hydride molecular weight
data.MW_H2    = 0.002;                                                     %hydrogen molecular weight
data.wt       = data.SC*data.MW_H2/data.MW_MH;                             %gravimetric density
data.kH       = 0.185;                                                     %hydrogen conductivity
data.cpH      = 14300;                                                     %hydrogen specific heat
data.rhoH     = 0.0899;                                                    %hydrogen density
data.K        = 1e-8;                                                      %LaNi5 permeability

%metal foam properties
data.rho_foam = 2700;
data.cp_foam = 963;
data.k_foam = 10.9;
data.foam_porosity = 0.91;                                                 %metal foam porosity

data.porosity = 0.50;                                                      %hydride porosity

data.rho_foam = (1-data.foam_porosity)*data.rho_foam;
data.rhoM     = data.foam_porosity*(1-data.porosity)*data.rhoMH;

data.kEff     = data.porosity*data.kH + (1-data.porosity)*data.kMH;
% data.kEff     = data.k_foam + data.foam_porosity*data.kMH + data.porosity*data.foam_porosity*data.kH;
data.capEff   = data.porosity*(data.rhoH*data.cpH) + ...
    (1-data.porosity)*(data.rhoMH*data.cpMH);
% data.capEff   = data.rho_foam*data.cp_foam + (data.rhoM + data.rhoH)*data.cpMH + data.rhoH_foam*data.cpH;
data.alfaEff  = data.kEff/data.capEff;                                     %hydride diffusivity

data.m_s      = data.V*data.porosity*data.rhoMH;
data.R        = 8.314;                                                     %ideal gas constant

%% PCM geometry and physical properties
data.N  = 61;       %n° nodi alla PCM (ca. 61)
data.r  = 40;       %n-esimo nodo corrispondente al fronte (ca. 40)
data.E  = 0.11;     %raggio esterno PCM (viene corretto più avanti)
data.r0 = data.D/2; %raggio interno PCM

data.Nreal = data.N+2*data.MH_nodes+1;                                     %teoricamente +1 ma considero già MH
data.rreal = data.r+2*data.MH_nodes+1;

%PCM physical properties
data.rhoPCM    = data.rhoMH*0.3;                                          %PCM density
data.cpPCM     = data.cpMH*4;                                              %PCM specific heat
data.kS        = data.kMH*0.6;                                                     %PCM conductivity (solid)
data.alfaS     = data.kS/(data.rhoPCM*data.cpPCM);
data.kL        = data.kS;                                                  %PCM conductivity (liquid)
data.alfaL     = data.kL/(data.rhoPCM*data.cpPCM);
data.lambdaPCM = 300e3;                                                    %PCM latent heat of fusion
data.beta      = 9.1e-4;                                                  %PCM thermal expansion coefficient
data.Tpcm      = 305;                                                      %PCM solidification temperature
data.mi        = 4e-3 * 1.25;                                                   %PCM dynamic viscosity
data.ni        = data.mi/data.rhoPCM;
data.Pr        = (data.mi*data.cpPCM)/(data.kL);

%% Absorption/desorption coefficients and various
%absorption data
data.deltaH_a = -30478;
data.Ea      = 21170; 
data.Ca       = 59.2;
data.deltaS_a = -108;
data.Tin      = 290;                                                       %Hydrogen initial temperature

%desorption
data.deltaH_d = 30800; 
data.Ed       = 16420; 
data.Cd       = 9.6; 
data.deltaS_d = 108;
% data.E         = sqrt(data.D^2/4 + ...
%     ( ((data.wt*sum(data.m_s)*data.deltaH_d/data.MW_H2)/data.lambdaPCM)/data.rhoPCM )...
%     /(pi*data.L));
data.E         = sqrt(data.D^2/4 + ...
    4*( ((data.wt*sum(data.m_s)*data.deltaH_d/data.MW_H2)/data.lambdaPCM)/data.rhoPCM )...
    /(pi*data.L));
% data.E        = round(data.E,2);
% data.E = 10*(data.D/2);                                                    

%various
data.Dvalve = 5e-3;                                                      %discharge valve surface area
data.Svalve = (pi*data.Dvalve^2/4);
data.sl     = 0.13;                                                        %slope coefficient
data.Pin = 30e5;
data.Pout = 1e5;




% portata limitata (o ciclo WLTP)
data.timeLimit = 100*60*60;                                                %tempo limite di desorbimento in secondi





%%Air
data.air.beta=3.41e-3;
data.air.Tamb=293;
data.air.Pr=0.71;
data.air.ni=1.5e-5;
data.air.rho=1.19;
data.air.k=2.56e-2;
end