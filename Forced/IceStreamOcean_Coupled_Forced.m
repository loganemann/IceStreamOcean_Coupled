%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Ice Stream Model With ELRA Bed Adjustment and Ocean Forcing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Logan E. Mann // Ice and Climate Group Georgia Tech

clear
clc
close all
p.phaShift=0;
yearCon = 3600*24*365;    %Seconds in a year (sec)
%% Primary Parameters
p.Ts=273.15-35;    % K
p.G =0.01;         % W m-2
p.On=true;         % Run with(true) or without(false) ocean forcing

p.xi=0;               % Heinrich Salinity Flux Parameter
p.mf=800/yearCon;       % Submarine Melt Rate (m/s)
%p.mf=200/yearCon;       % Submarine Melt Rate (m/s)


p.gamma=200/(3e5*yearCon); % Relaxation Time Parameter


%% Initial Conditions and time step

p.h_init=1000;          % Initial avg ice stream height (m)
p.e_init=0.6;           % Initial Void Ratio
p.Zs_init=1;            % Initial Till Height (m)
p.Tb_init=273.15;         % Initial Basal Temperature (K)
p.L_init=3000e3;
p.x_o_init=0.2;             % Initial Ocean Temperature (nondimensional)
p.y_o_init=0.2;             % Salinity (nondimensional)

p.t_final = 3e5;          %Total time of integration (yr)
p.tspan=[0 yearCon*p.t_final];


%% Parameters

% Ocean Oscillator Parameters
p.T_0=5;                % Temperature of NADW (K)
p.T_A=-15;              % Relaxation Temperature
%p.T_A=p.Ts-273.15;
p.delta=.001;
p.epsilon=-0.01;
p.nu1=0.1;
p.nu2=5;
p.Td=-200;
%p.Td=-230;

% Ice Stream Parameters
p.a = 9.44e8;   % Till empirical coefficient (Pa)
p.ac= 0.06/yearCon;      % Accumulation rate (m sec-1)
p.Ag= 5e-25;    % Glen's law rate factor (Pa-3 s-1)
p.b = 21.7;     % Till empirical exponent (unitless)
p.Ci= 1.94e6;   % Volumetric heat capacity of ice (J K-1 m-3)
p.ec= 0.3;      % Till consolidation threshold (unitless)
p.g = 9.81;     % Acceleration (s-2)
p.hb= 10;       % Thickness of temperate ice layer (m)
p.ki= 2.1;      % Thermal conductivity of ice (J s-1 m-1 K-1)
p.L = 400e3;      % Ice stream trunk length (m)
p.Lf= 3.35e5;   % Specific latent heat of ice (J kg-1)
p.n = 3;        % Glen's law exponent (unitless)
p.W = 180e3;       % Ice stream trunk width (m)
p.ws= 1;        % Till saturation threshold (m)
p.Zs_init= 1;        % Initial effective till layer thickness (m)
p.Zs_min=1e-3;   % minimum till thickness
p.rhoi= 917;      % Ice Density (kg m-3)
p.rhow= 1000;     % Water density (kg m-3)
p.Tm= 273.15;    % Melting point of water (K)

% Grounding Line Parameters
p.A=1.0e-25;             %Glen's Law Exponent (s-1 Pa-3)
p.C=7.624e6;            % Basal friction coefficient (Pa m^-1/n s^1/n)
p.f=0.4;                % Coulumb Friction coefficient
p.m=1/3;                % Weertman friction law exponent
p.P=0.3/yearCon;        % Time-averaged accumulation rate (m/s)
%p.sigmaP=0.1/yearCon;   % Accumulation rate variance (m/year)
p.alpha=7;              % Interior Ice flux thickness exponent
%p.gamma=3;              % Interior ice flux length exponent
p.theta=0.6;            % Buttressing parameter
p.lambda=p.rhow/p.rhoi;
p.omega=(1/10)*(p.Ag *(p.rhoi*p.g)^p.n+1 * (p.theta*(1-p.lambda^-1))^p.n * (4^p.n*p.C)^-1)^(1/(p.m+1));
p.Lc=500e3;

% Isostatic Parameters
p.b0=500;               % Ice divide bed height(m)
p.bx=-5e-4;              % Prograde bed slope
p.nx=1001;                 % n grid spaces
p.xmax=10000e3;
p.dx=p.xmax/(p.nx-1);             % Grid Spacing
p.x=linspace(-p.xmax,p.xmax,p.nx);

p.sig=100e3;            % Sill width 
p.S_mu=3100e3;              % Sill position (m)
p.B0(p.x<0)=p.b0;       
p.B0(p.x>=0 & p.x<=p.S_mu+200e3)=p.b0+(p.bx/2).*p.x(p.x>=0 & p.x<=p.S_mu+200e3);
p.B0(p.x>p.S_mu+200e3)=p.B0(p.x==p.S_mu+200e3)+p.bx.*(p.x(p.x>p.S_mu+200e3)-p.x(p.x==p.S_mu+200e3));

p.B0(3*floor(length(p.B0)/4):end)=p.B0(3*floor(length(p.B0)/4));
p.B0=p.B0+0.5e8/p.sig/sqrt(2*pi)*exp(-0.5*((p.x-p.S_mu)/p.sig).^2);

p.rhoa=2830;            % Density of aesthenosphere (kg m-3)
p.D=5e24;               % Flexural Rigidity (N m)
p.tau=3000*yearCon;             % Relaxation Time (s)

diag=ones(p.nx,1).*[p.D/p.dx^4, -4*p.D/p.dx^4, p.rhoa*p.g+6*p.D/p.dx^4, -4*p.D/p.dx^4, p.D/p.dx^4];
p.M=spdiags(diag,-2:2,p.nx,p.nx);

% Implement Boundary Conditions
p.M(1,1:3)=[p.rhoa.*p.g 0 0];
p.M(2,1:4)=[0 p.rhoa.*p.g 0 0];
p.M(end-1,(end-3):end)=[0 0 p.rhoa.*p.g 0];
p.M(end,(end-2):end)=[0 0 p.rhoa.*p.g];



%% Calculate Unloaded Bed Topography
bg=interp1(p.x,p.B0,p.L_init);
hg=-p.lambda*bg;
hx=[p.h_init.*ones(length(p.x(p.x<0)),1); interp1([0 p.L_init],[p.h_init hg],p.x(p.x>=0 & p.x<=p.L_init))';  zeros(length(p.x(p.x>p.L_init)),1)];
sig_zz=-p.rhoi.*p.g.*hx;
sig_zz(p.x>p.L_init)=-p.rhow.*p.g.*(-p.B0(p.x>p.L_init));
w=p.M\sig_zz;

p.u0=p.B0-w';
p.u0=p.B0;
%p.u_init=w';
p.u_init=zeros(1,p.nx);

%% ODE solve
options=odeset('RelTol',1e-7,'AbsTol',1e-7);
init=[p.h_init p.e_init p.Zs_init p.Tb_init p.L_init p.x_o_init p.y_o_init p.u_init];
[time,soln]=ode113(@(t,X)IceStreamOcean_Coupled_Forced_RHS(t,X,p),p.tspan,init,options);
%[time,soln]=ode113(@(t,X)IceStreamOcean_Coupled_Forced_RHS(t,X,p),linspace(0,p.tspan(end),5000),init,options);


h=soln(:,1);
e=soln(:,2);
Zs=soln(:,3);
Tb=soln(:,4);
L=soln(:,5);
x_o=soln(:,6);
y_o=soln(:,7);
u=soln(:,8:end);
%% Diagnostic Equations

% Isostatic Equations
B=u+p.u0.*ones(length(time),length(p.u0));
bg=zeros(length(time),1);
for i=1:length(time)
    bg(i,1)=interp1(p.x,B(i,:),L(i));
end
hg=-p.lambda.*bg;

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

bgi=interp1(p.x,p.B0,L')';


% Ocean Forcing Pulse
mu1=p.nu1./(1+p.nu1)+p.epsilon*p.nu1;
mu2=p.nu2./(1+p.nu2)+p.epsilon*p.nu2;
mu=mu2-p.delta*p.nu2;

nu=p.nu2*ones(length(time),1);
nu(y_o-x_o<=p.epsilon)=p.nu1;

To=max([x_o.*(p.T_A-p.T_0)+p.T_0, zeros(length(time),1)]')';
%To= 3.125;

S=B(:,p.x==p.S_mu);

Qm=p.mf.*To.*(-bg);
Qm((S>p.Td & L<p.S_mu) | bg>p.Td)=0;
Qg=Qg+Qm;           % Add Melt Flux to GLine flux



