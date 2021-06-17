function rhs = IceStreamOcean_Coupled_Forced_RHS(t,X,p)
%% Ice Stream RHS
%% Read Parameters
phaShift=p.phaShift;
yearCon = 3600*24*365;    %Seconds in a year (sec)


% Ocean Parameters
T_0=p.T_0;
T_A=p.T_A;
delta=p.delta;
epsilon=p.epsilon;
nu1=p.nu1;
nu2=p.nu2;
gamma=p.gamma;


Td=p.Td;

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
A=p.A;
C=p.C;
f=p.f;
m=p.m;
P=p.P;
%sigmaP=p.sigmaP;
alpha=p.alpha;
%gamma=p.gamma;
theta=p.theta;
omega=p.omega;
lambda=p.lambda;
Lc=p.Lc;

% Isostatic Parameters
b0=p.b0;               % Ice divide bed height(m)
bx=p.bx;              % Prograde bed slope
nx=p.nx;                 % n grid spaces
dx=p.dx;             % Grid Spacing
x=p.x;
u0=p.u0;
rhoa=p.rhoa;            % Density of bedrock (kg m-3)
D=p.D;
tau=p.tau;              %Relaxation time
M=p.M;



sig=p.sig;            % Sill width 
S_mu=p.S_mu;              % Sill position (m)
%% Read Additional Parameters
Ts=p.Ts;    % K
G = p.G;         % W m-2

mf=p.mf;
xi=p.xi;
%% Set Variables
h=X(1);
e=X(2);
Zs=X(3);
Zs=max([Zs Zs_min]); Zs=min([Zs Zs_init]);
Tb=min([X(4) Tm]);
L=X(5);
x_o=X(6);
y_o=X(7);
u=X(8:end);

%% Diagnostic Equations

% Isostatic Rebound
B=u'+u0;
bg=interp1(x,B,L);
hg=-p.lambda*bg;
hx=[0.*ones(length(x(x<0)),1); interp1([0 L],[h hg],x(x>=0 & x<=L))';  zeros(length(x(x>L)),1)];

sig_zz=-rhoi.*g.*hx;
sig_zz(x>L)=-rhow.*g.*(-B(x>L));

% hold off
% plot(p.x,B,'k')
% hold on
% plot([L L],[bg hg],'k')
% plot([0 L],[h hg],'k')
% plot(p.x,zeros(1,length(p.x)),'b')
% axis([0 4000e3 -1000 4500])
% pause(0.01)

if length(sig_zz)~=length(M(:,1))
    disp('help')
    pause(10)
end
    
w=M\sig_zz;


% Basic Model Equations
e=max([e ec]);
taub=min([a*exp(-b*(e)), f*(rhoi*g*h)]);      %Basal Shear Stress (Tsai, 2015) (Power Law)
taud=rhoi*g*h^2/L;
ub=(Ag/256)*W^(n+1)/(4^n*(n+1)+h^n)*max([taud-taub 0])^n;

% Grounding Line Equations
Qg=0.61*(8*Ag*(rhoi*g)^n)/(4^n*f) *(1-rhoi/rhow)^(n-1) * hg^(n+2);   %Grounding Line Flux (Tsai, 2015)
Qd=(2*Ag/(n+2) * h^2 .* min([taub taud])^n);      % Flux from deformation
Qv=ub*h;                                        % Flux from motion
Q=Qv+Qd;


% Ocean Forcing Pulse
mu1=nu1/(1+nu1)+epsilon*nu1;
mu2=nu2/(1+nu2)+epsilon*nu2;
mu=mu2-delta*nu2;

if y_o-x_o<=epsilon; nu=nu1;
else; nu=nu2; end

To=max(x_o.*(T_A-T_0)+T_0,0);

%To=3.125;

S=B(x==S_mu);
if (S<Td && bg<Td) || (L>S_mu && bg<Td); Qm=mf*To*(-bg); else; Qm=0; end
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

dudt=-(u-w)./tau;

dx_o_dt=(1-x_o-nu*x_o)*gamma;
dy_o_dt=(mu*(1-xi*Qg)-nu*y_o)*gamma;



%% RHS out
rhs=[dhdt; dedt; dZsdt; dTbdt; dLdt; dx_o_dt; dy_o_dt; dudt];
disp(['Percent Integration Complete:' num2str(100*t./p.tspan(end)) '%']);
%disp(['Percent Integration Complete:' num2str(100*100*(p.counter/p.iter + t/p.tspan(end)/p.iter^2)), '%'])
%/100+p.i-1
end