%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 2 Directional Coupled Ice Stream Ocean Model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Logan E. Mann // Ice and Climate Group Georgia Tech

clear
clc
close all
yearCon = 3600*24*365;    %Seconds in a year (s)
%% Primary Parameters
p.Ts=273.15-38;    % K
p.G =0.04;         % W m-2

p.xi=10;               % Heinrich Salinity Flux Parameter
p.mf=10/yearCon;       % Submarine Melt Rate (m/s)
p.gamma=260/(3e5*yearCon); % Relaxation Time Parameter

%% Initial Conditions and time step
p.h_init=3000;          % Initial avg ice stream height (m)
p.e_init=0.6;           % Initial Void Ratio
p.Zs_init=1;            % Initial Till Height (m)
p.Tb_init=273.15;         % Initial Basal Temperature (K)
p.L_init=1000e3;         % Initial Grounding Line Position (m)
p.x_init=0.2;             % Initial Ocean Temperature (nondimensional)
p.y_init=0.2;             % Salinity (nondimensional)

p.t_final = 3e5;          %Total time of integration (yr)
p.tspan=[0 yearCon*p.t_final];

%% Ocean Oscillator Parameters
p.T_0=5;                % Temperature of NADW (K)
p.T_A=-15;              % Relaxation Temperature
p.delta=.001;
p.epsilon=-0.01;
p.nu1=0.1;
p.nu2=5;


%% Parameters

% Ice Stream Parameters
p.a = 9.44e8;   % Till empirical coefficient (Pa)
p.ac= 0.1/yearCon;      % Accumulation rate (m sec-1)
p.Ag= 5e-25;    % Glen's law rate factor (Pa-3 s-1)
p.b = 21.7;     % Till empirical exponent (unitless)
p.Ci= 1.94e6;   % Volumetric heat capacity of ice (J K-1 m-3)
p.ec= 0.3;      % Till consolidation threshold (unitless)
p.g = 9.81;     % Acceleration (s-2)
p.hb= 10;       % Thickness of temperate ice layer (m)
p.ki= 2.1;      % Thermal conductivity of ice (J s-1 m-1 K-1)
p.Lf= 3.35e5;   % Specific latent heat of ice (J kg-1)
p.n = 3;        % Glen's law exponent (unitless)
p.W = 40e3;       % Ice stream trunk width (m)
p.ws= 1;        % Till saturation threshold (m)
p.Zs_init= 1;        % Initial effective till layer thickness (m)
p.Zs_min=1e-3;   % minimum till thickness
p.rhoi= 917;      % Ice Density (kg m-3)
p.rhow= 1028;     % Seawater density (kg m-3)
p.Tm= 273.15;    % Melting point of water (K)

% Grounding Line Parameters
p.f=0.4;                % Coulomb Friction coefficient

% Bed Parameters
p.b0=100;               % Ice divide bed height(m)
p.bx=-5e-4;              % Prograde bed slope



%% ODE solve
options=odeset('RelTol',1e-9,'AbsTol',1e-9);
init=[p.h_init p.e_init p.Zs_init p.Tb_init p.L_init p.x_init p.y_init];
%[time,soln]=ode113(@(t,X)IceStreamOcean_Coupled_2Way_RHS(t,X,p),linspace(0,p.tspan(end),10000),init,options);
[time,soln]=ode113(@(t,X)IceStreamOcean_Coupled_2Way_RHS(t,X,p),p.tspan,init,options);

h=soln(:,1);
e=soln(:,2);
Zs=soln(:,3);
Tb=soln(:,4);
L=soln(:,5);
x=soln(:,6);
y=soln(:,7);
%% Diagnostic Equations

% Isostatic Equations
bg=p.b0+p.bx.*L;
hg=-p.rhow/p.rhoi.*bg;

% Basic Model Equations
e=max([e p.ec.*ones(length(time),1)]')';
taub=min([p.a.*exp(-p.b.*e), p.f.*(p.rhoi.*p.g.*h)]')';
taud=p.rhoi.*p.g.*h.^2./L;
ub=(p.Ag./256).*p.W.^(p.n+1)./(4.^p.n.*(p.n+1)+h.^p.n).*max([taud-taub zeros(length(time),1)]')'.^p.n;

% Grounding Line Equations
Qg=0.61.*(8.*p.Ag.*(p.rhoi.*p.g).^p.n)./(4.^p.n.*p.f) .*(1-p.rhoi./p.rhow).^(p.n-1) .* hg.^(p.n+2);   %Grounding Line Flux (Tsai, 2015)
Qd=(2.*p.Ag./(p.n+2) .* h.^2 .* min([taub taud]')'.^p.n);      % Flux from deformation
Qv=ub.*h;                                        % Flux from motion
Q=Qv+Qd;


% Ocean Forcing Pulse
mu1=p.nu1./(1+p.nu1)+p.epsilon*p.nu1;
mu2=p.nu2./(1+p.nu2)+p.epsilon*p.nu2;
mu=mu2-p.delta*p.nu2;

nu=p.nu2*ones(length(time),1);
nu(y-x<=p.epsilon)=p.nu1;

To=max([x.*(p.T_A-p.T_0)+p.T_0, zeros(length(time),1)]')';
%To=x.*(p.T_A-p.T_0)+p.T_0;
Qm=p.mf.*To.*(-bg);
Qg=Qg+Qm;           % Add Melt Flux to GLine flux


%% Plotting
figure(1)
subplot(5,2,1)
hold on
plot(time/yearCon,h/1000)
xlabel('time (yr)')
ylabel('h (km)')

subplot(5,2,2)
hold on
plot(time/yearCon,L/1000)
xlabel('time (yr)')
ylabel('L (km)')

subplot(5,2,3)
hold on
plot(time/yearCon,e)
xlabel('time (yr)')
ylabel('e')

subplot(5,2,4)
hold on
plot(time/yearCon,Q)
xlabel('time (yr)')
ylabel('Q')

subplot(5,2,5)
hold on
plot(time/yearCon,Qg)
xlabel('time (yr)')
ylabel('Qg')

subplot(5,2,6)
hold on
plot(time/yearCon,ub*yearCon/1000)
xlabel('time (yr)')
ylabel('ub (km/yr)')

subplot(5,2,7)
hold on
plot(time/yearCon,hg)
xlabel('time (yr)')
ylabel('h_{g} (m)')

subplot(5,2,8)
hold on
plot(time/yearCon,Tb)
xlabel('time(yr)')
ylabel('Tb')

subplot(5,2,9)
hold on
plot(time/yearCon,bg)
xlabel('time(yr)')
ylabel('b_{g}')



%% Ocean Forcing Figures
figure(2)
subplot(3,1,1)
hold on
plot(time/yearCon, To)
xlabel('time(yr)')
ylabel('T_{o} (^{o}C)')

subplot(3,1,2)
hold on
plot(time/yearCon, L/1000)
xlabel('time (yr)')
ylabel('L (km)')

subplot(3,1,3)
hold on
plot(time/yearCon,Qm*yearCon)
xlabel('time(yr)')
ylabel('Q_{m}(m^{2}/a)')
sgtitle(['G= ' num2str(p.G) ' W/m^{2}, T_{s}= ' num2str(p.Ts-273.15) ' ^{o}C, \xi=' num2str(p.xi) ', m_{f}=' num2str(p.mf*yearCon) ' m/yr^{o}C'])

%% Phase Analysis

dhdt=(h(2:end)-h(1:end-1))./(time(2:end)-time(1:end-1));
dhdt(dhdt>min(dhdt)-min(dhdt)*3/4)=0;

% phi_H=[];
% for i=1:length(dhdt)-1           % Identify the peaks of ub as a value that is not 0 adjacent to a 0
%     if dhdt(i)==0 && dhdt(i+1)~=0
%         phi_H=[phi_H; i+1];
%     end
% end

[~,phi_H]=findpeaks(-dhdt);
phi_H(1)=[];

figure(5)
plot(time(1:end-1)/yearCon,dhdt)
hold on
plot(time(phi_H)/yearCon,dhdt(phi_H),'x')

[~,phi_DO_all]=findpeaks(To);

phi_DO=zeros(length(phi_H),1);
for indx=1:length(phi_H)
    phi_DO_temp=phi_DO_all(phi_DO_all<phi_H(indx));
    if isempty(phi_DO_temp); phi_DO_temp=1; end
    phi_DO(indx)=phi_DO_temp(end);
end
dphi=time(phi_H)-time(phi_DO);   % Phase Difference between Heinrich and DO

T_DO=mean(time(phi_DO_all(2:end))-time(phi_DO_all(1:end-1)));
T_H=mean(time(phi_H(2:end))-time(phi_H(1:end-1)));
%dphi(dphi>T_DO/2)=dphi(dphi>T_DO/2)-T_DO;




%% Plotting
extent=100;
figure(3)
subplot(2,1,1)
yyaxis left
%plot(time/yearCon,ub*yearCon/1000,'b')
plot(time/yearCon/1000,h,'b')
%axis([0 p.t_final 0 0.5])
hold on
%plot(time(phi_H)/yearCon,ub(phi_H),'xk')
plot(time(phi_H)/yearCon/1000,h(phi_H),'xk')

xlabel('time (yr)')
ylabel('h (m)')

yyaxis right
plot(time/yearCon/1000,To,'r','LineWidth',1)
hold on
plot(time(phi_DO)/yearCon/1000,To(phi_DO),'xr')
axis([0 extent 0 2])
ylabel('T_{o}')

subplot(2,1,2)
plot(time(phi_H)/yearCon/1000,dphi/yearCon,'.','MarkerSize', 10)
axis([0 extent min(dphi/yearCon) max(dphi/yearCon)])
sgtitle([num2str((time(phi_H(end))-time(phi_H(end-1)))/(T_DO)) ':' '1' ', \xi=' num2str(p.xi) ', m_{f}=' num2str(p.mf*yearCon) ' m/yr^{o}C'])
xlabel('time (yr)')
ylabel('\Phi_{h} - \Phi_{T_{o}} (Phase Difference in years)')



