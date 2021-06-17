function rhs = IceStreamOcean_Stochastic_RHS(t,X,p)
%% Ice Stream RHS

dt=p.dt;
%% Ocean Parameters
T_0=p.T_0;
T_A=p.T_A;
delta=p.delta;
epsilon=p.epsilon;
nu1=p.nu1;
nu2=p.nu2;
gamma=p.gamma;
xi=p.xi;
std=p.std;

%% Read Parameters
yearCon = 3600*24*365;    %Seconds in a year (sec)
% Ice Stream Parameters
a=p.a;
ac=p.ac;
Ag=p.Ag;
b=p.b;
Ci=p.Ci;
ec=p.ec;
g=p.g;
hb=p.hb;
ki=p.ki;
Lf=p.Lf;
n=p.n;
W=p.W;
ws=p.ws;
Zs_init=p.Zs_init;
Zs_min=p.Zs_min;
rhoi=p.rhoi;
rhow=p.rhow;
Tm=p.Tm;

% Grounding Line Parameters
f=p.f;              % Coulumb Friction Coefficient

% Bed Parameters
b0=p.b0;               % Ice divide bed height(m)
bx=p.bx;              % Prograde bed slope


% Ocean Forcing Parameters
mf=p.mf;
%% Read Additional Parameters
Ts=p.Ts;    % K
G = p.G;         % W m-2

%% Set Variables
h=X(1);
e=X(2);
Zs=X(3);
Zs=max([Zs Zs_min]); Zs=min([Zs Zs_init]);
Tb=min([X(4) Tm]);
L=X(5);
x=X(6);
y=X(7);

%% Diagnostic Equations

% Bed Eqs
bg=b0+bx*L;
hg=-rhow/rhoi*bg;

% Basic Model Equations
e=max([e ec]);
taub=a*exp(-b*(e-ec));
taud=rhoi*g*h^2/L;
ub=(Ag/256)*W^(n+1)/(4^n*(n+1)+h^n)*max([taud-taub 0])^n;
%ub=(Ag)*W^(n+1)/(4^n*(n+1)+h^n)*max([taud-taub 0])^n;

% Grounding Line Equations
Qg=0.61*(8*Ag*(rhoi*g)^n)/(4^n*f) *(1-rhoi/rhow)^(n-1) * hg^(n+2);   %Grounding Line Flux (Tsai, 2015)
Qd=(2*Ag/(n+2) * h^2 .* min([taub taud])^n);      % Flux from deformation
Qv=ub*h;                                        % Flux from motion
Q=Qv+Qd;

% Ocean Forcing Pulse
mu1=nu1/(1+nu1)+epsilon*nu1;
mu2=nu2/(1+nu2)+epsilon*nu2;
mu=mu2-delta*nu2+std*randn/sqrt(dt);

if y-x<=epsilon; nu=nu1;
else; nu=nu2; end

To=max(x.*(T_A-T_0)+T_0,0);
%To=x.*(T_A-T_0)+T_0;

Qm=mf*To*(-bg);
Qg=Qg+Qm;           % Add Melt Flux to GLine flux


%% Prognostic Equations
if ((Zs==Zs_min && Tb==273.15 && ((taub*ub) + G + (ki*(Ts-Tb)/h))<0) || (Zs==Zs_min && Tb<273.15))   %till is frozen
    ub=0;
    dedt=0;
    dZsdt=0;
    dTbdt=(1/(hb*Ci))*((taub*ub) + G + (ki*(Ts-Tb)/h));
elseif ((e == ec && Zs==Zs_init && ((taub*ub) + G + (ki*(Ts-Tb)/h))<0) ||  (e == ec && Zs<Zs_init))
    ub=0;
    dedt=0;
    dZsdt=((taub*ub) + G + (ki*(Ts-Tb)/h))/(Lf*rhoi);
    dTbdt=0;
else
    dedt = ((taub*ub) + G + (ki*(Ts-Tb)/h))/(Zs*Lf*rhoi);
    dZsdt = 0;
    dTbdt = 0;
end

dhdt=ac - Qg/L - h*(Q-Qg)/(hg*L);
dLdt=(Q-Qg)/hg;

dxdt=(1-x-nu*x)*gamma;
dydt=(mu*(1-xi*Qg)-nu*y)*gamma;

%% RHS out
rhs=[dhdt; dedt; dZsdt; dTbdt; dLdt; dxdt; dydt];

disp(['Percent Integration Complete:' num2str(100*t./p.tspan(end)) '%']);

end