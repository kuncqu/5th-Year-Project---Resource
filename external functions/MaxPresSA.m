function [Out] = MaxPresSA(Vinf,Tinf,Tw)

%Function specifically for the Sensitivity Test
%Assumption of air comprising of 80% N2 and 20% O2
%Other Nominal Parameters 
%   Tinf = 300 K
%   Tw = 1000 K
%   Vinf = 7500 m/s
%   Total Number Density = 5e18;

% constants
% dAir = 4.19e-10;
C1 = 1.458e-6;          %Sutherland's Law Coefficient, kg/m-s-K^-0.5
C2 = 110.4;             %Sutherland's Law Coefficient, K
omega = 0.78;           %Temperature Viscosity power law coefficient
k = 1.3806488e-23;      %Boltzmann constant, m2 kg s-2 K-1
Nav = 6.023e23;         %Avogadro Constant 
amu = 1.66053892e-27;
g = 1.4;                %heat capacity ratio
cp = 1004.7;            %specific heat at constant pressure (AIR)
% cp = 1373;            %specific heat at constant pressure
R = cp * (g-1) / g;     %Gas constant, kg m2 s-2 K-1 mol-1
c = sqrt (g*R*Tinf);    %Speed of sound, m s-1
Minf = Vinf / c;        %Mach number

mN2 = 28;               %molar mass of nitrogen molecule, grams/mole
mO2 = 32;               %molar mass of oxigen molecule, grams/mole

%Atmospheric Composition

T_ND = 5.0e18;          %?????total mass of the atmosphere??????

X_N2 = 0.8;             %molar fraction of nitrogen
X_O2 = 0.2;             %molar fraction of oxigen

ND_N2 = T_ND * X_N2;
ND_O2 = T_ND * X_O2;

rho = ND_N2/Nav * mN2/1000 * X_N2 + ND_O2/Nav * mO2/1000 * X_O2;  %density, kg/m^3
P = rho * R * Tinf;
rhow = P/R/Tw;  %Assuming constant pressure in the boundary layer

% Free Stream Stagnation Pressure and Temperature
P01 = P * (1 + 0.5 * (g-1) * norm(Minf)^2) ^ (g/(g-1)); 
T0 = Tinf * (1 + 0.5 * (g-1) * norm(Minf)^2);
mu_Tinf = (C1 * Tinf^(3/2))/(Tinf + C2);
mu_T0 = mu_Tinf * (T0/Tinf)^omega;
rhos = rho*(T0/Tinf)^(1/(g-1));
mu_w = (C1 * Tw^(3/2))/(Tw + C2);

hw = cp * Tw;
h0 = cp * T0;

% Stagnation Pressure after Normal Shock
P02 = P01 * (((g+1)*norm(Minf)^2)/((g-1)*norm(Minf)^2+2))^(g/(g-1)) * ((g+1)/(2*g*norm(Minf)^2-(g-1)))^(1/(g-1));
Cpmax = 2 * (P02 - P) / rho / norm(Vinf) ^ 2;   %pressure coefficient inf-02

%% FMF properties %%

mg = X_N2 * mN2 * amu + X_O2 * mO2 * amu;
vmp = sqrt(2*k*Tinf/mg);
s = norm(Vinf)/vmp;

%Molecular Diameters at 273K from Bird, pg 409

dN2 = 4.17e-10;     %m
dO2 = 4.07e-10;     %m

%Calculating mean Free Path

LN2 = k*Tinf/sqrt(2)/pi()/dN2^2/P;     %mean free path for N2, m
LO2 = k*Tinf/sqrt(2)/pi()/dO2^2/P;     %mean free path for O2, m
L_b = X_N2*LN2 + X_O2*LO2;     %Overall mean free path for 80% N2 and 20% O2

Out = [Cpmax,norm(Minf),s,L_b,rho,g,mu_T0,T0,P,P01,rhos,mu_Tinf,mu_w,rhow,hw,h0];

end

