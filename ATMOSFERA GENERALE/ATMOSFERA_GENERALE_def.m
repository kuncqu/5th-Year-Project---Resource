addpath('C:\Users\GBenedetti\Google Drive\TESI LM\tesi\tps\FOSTRAD\aerothermal-master\aerothermal-master\altmany-export_fig-2bad961');
addpath('C:\Users\GBenedetti\Google Drive\TESI LM\tesi\tps\FOSTRAD\aerothermal-master\aerothermal-master\external functions');
clear all
close all
%%%%%% User Inputs %%%%%%%

SPAZAL=1; %INTERVALLO DI SPAZIATURA PER ALPHA
SPAZH=1; %INTERVALLO DI SPAZIATURA PER H km (lasciare fisso!)
alphamax=20;
alphamin=20;
hmin=0; %km
hmax=300; %km
beta = 0;       % Yaw angle, degrees

limKn_inf=5*10^-4;
limKn_sup=5;
Vinf = 5000;     % Free Stream Velocity, m/s

Lat= 55.8642;   % Latitude
Lon= -4.16;     % Longitude
year= 2010;     % year
yr_day= 1;      % day
UT_sec= 0;      % UT seconds

b=1; %char length

%%% Atmosphere model

[Ti,rhoi] = atmosnrlmsise00([hmin*1e3:SPAZH*1e3:hmax*1e3] , Lat, Lon, year, yr_day, UT_sec);
rhoi(:,[6,9]) = [ ]; %He O N2 O2 Ar H N [1/m^3] gas densities
Ti(:,1)= [];

SN = 1.0;   %Normal Momentum Accommodation coefficient 
ST = 1.0;   %Tangential Momentum Accommodation Coefficient
AT = ST*(2-ST);     %Tangental Energy Accomodation Coefficient
AN = 1 - (1 - SN)^2;    %Normal Energy Accommodation Coefficient
AC = 0.5 * (AN + AT);   %Overall Energy Accommodation Coefficient

%specific heat at constant pressure (AIR)
cp = 1004.7;         

% heat capacity ratio
percentual= rhoi(:,1:7)./  ((rhoi(:,1:7)* ones(7,1)) * ones (1,7)) ; %He O N2 O2 Ar H N

gamma= percentual* [ 5/3 5/3 1.4 1.4 5/3 5/3 5/3]';

%Gas constant, kg m2 s-2 K-1 mol-1
R = cp * (gamma-1) ./ gamma;     

%Temperature Viscosity power law coefficient
omega = 0.78;           

C1 = 1.458e-6;          %Sutherland's Law Coefficient, kg/m-s-K^-0.5
C2 = 110.4;             %Sutherland's Law Coefficient, K

%molecular weight 
mN2 = 28.01340;               %molar mass of nitrogen molecule, grams/mole
mO2 = 31.99880;               %molar mass of oxigen molecule, grams/mole
mO = mO2/2;                   %molar mass of oxigen atom, grams/mole          
mN = mN2/2;                   %molar mass of Nitrogen atom, grams/mole
mAr = 39.9480;                %molar mass of Argon molecule, grams/mole
mHe = 4.0026020;              %molar mass of helium molecule, grams/mole
mH= 1.007940;                 %molar mass of Hydrogen molecule, grams/mole

xN2=percentual(:,3);
xO2=percentual(:,4);
xO= percentual(:,2);
xN= percentual(:,7);
xAr= percentual(:,5);
xHe= percentual(:,1);
xH= percentual(:,6);

%Atmosphere molecular weight
m=mN2.*xN2+mO2.*xO2+mO.*xO+mN.*xN+mAr.*xAr+mHe.*xHe+mH.*xH;

% FMF properties %
amu = 1.66053892e-27;
mg = m.*amu;
Kb = 1.3806488e-23;      %Boltzmann constant, m2 kg s-2 K-1

Nav = 6.023e23;                                     %Avogadro Constant 
f=1.660539040e-27;                                  %conversion factor atomic mass unit - kg


rN = 1;


h=[hmin:1:hmax];%Km

[Z, Z_L, Z_U, T, P, rho, c, g, mu, nu, k, n, n_sum] = atmo(hmax,1,1);

Z(1:hmin,:) = [] ;
T(1:hmin,:) = [] ;
P(1:hmin,:) = [] ;
rho(1:hmin,:) = [] ;
%c(1:hmin,:) = [] ;
g(1:hmin,:) = [] ;


% nu(1:hmin,:) = [] ;
% k(1:hmin,:) = [] ;

% gamma atm comp us76
percentual76= n(:,1:6)./  ((n(:,1:6)* ones(6,1)) * ones (1,6)) ; %N2 O O2 Ar He H
gamma76= percentual76* [ 1.4 5/3 1.4 5/3 5/3 5/3]';

x76N2=percentual76(:,1);
x76O2=percentual76(:,3);
x76O= percentual76(:,2);
x76Ar= percentual76(:,4);
x76He= percentual76(:,5);
x76H= percentual76(:,6);

%Atmosphere molecular weight
j=1;
k=1;
for i=1:hmax+1
    if h(i)<=85;
    m761(k)=m(k);
    k=k+1;
    else     
    m762(j)=mN2.*x76N2(j)+mO2.*x76O2(j)+mO.*x76O(j)+mAr.*x76Ar(j)+mHe.*x76He(j)+mH.*x76H(j);
    j=j+1;
    end
end
m76=[m761 m762]';

mg76 = m76.*amu;



%   Thermal Conductivity Coefficient
k = 2.64638e-3*Ti.^1.5./(Ti+245*10.^(-12./Ti));
k76=2.64638e-3*T.^1.5./(T+245*10.^(-12./T));

%dynamic Viscosity
        if h(end)<=85;
        else
%             mu= [ mu' , (mu(end)*ones ( hmax-86 ,1))' ]';
        mu = 1.458e-6.*T.^1.5./(T+110.4);
            
        end
mu(1:hmin,:) = [] ;


%P00
for i=1:hmax+1
P00(i)=atmo_p(i, Ti(i), sum(rhoi(i,:)));
end
P00=P00';

% rho nmsise00

if h(end)<=85;
rho00 = m.*P00./(R.*Ti);
else
rho00 = rhoi*[mHe mO mN2 mO2 mAr mH mN]'./6.022169e26;
end

%sound speed 
c=sqrt((gamma.*P00)./rho00);

c76=sqrt((gamma76.*P(length(h)-length(n)))./rho(length(h)-length(n)));

%Gas constant, kg m2 s-2 K-1 mol-1
cp = 1004.7;            %specific heat at constant pressure (AIR)
R = cp .* (gamma-1) ./ gamma;   

%mach number
Minf = norm(Vinf) ./ c;  

% Free Stream Stagnation Pressure and Temperature
P01 = P .* (1 + 0.5 .* (gamma-1) .* Minf.^2) .^ (gamma./(gamma-1)); 
T0 = Ti .* (1 + 0.5 .* (gamma-1) .* Minf.^2);

Tw = Ti.*(1+(gamma-1)./(2*(Minf).^2)); %Surface Temperature (original
% formula)


mu_T0 = mu.* (T0./Ti).^omega;
rhos = rho.*(T0./Ti).^(1./(gamma-1));
rhow = P./R./Tw;  %Assuming constant pressure in the boundary layer
mu_w = (C1 .* Tw.^(3/2))./(Tw + C2);
hw = cp .* Tw;
h0 = cp .* T0;



% Stagnation Pressure after Normal Shock
P02 = P01 .* (((gamma+1).*Minf.^2)./((gamma-1).*Minf.^2+2)).^(gamma./(gamma-1)) .* ...
    ((gamma+1)./(2.*gamma.*Minf.^2-(gamma-1))).^(1./(gamma-1));

%Cp max
Cpmax=(2./(gamma.*Minf.^2)).*((P02./P-1)); %Cp maximum

% FMF properties %
mg = m.*amu;
kb = 1.3806488e-23;      %Boltzmann constant, m2 kg s-2 K-1
vmp = sqrt(2.*kb.*Ti./mg);
s = Vinf./vmp;

%mu00
mu00 = 1.458e-6*Ti.^1.5./(Ti+110.4);

%Calculating mean Free Path and Knudsen

l=(mu00./P).*sqrt(pi.*kb.*Ti./(2.*mg));  %mean free path [m]
l76=(mu./P).*sqrt(pi.*kb.*T./(2.*mg76));  %mean free path [m]

Kn = l./b;                             %Knudsen Number
Kn76 = l76./b;

Pr=mu00.*cp./k;                           %Prandtl number
Pr76=mu.*cp./k76; 

Re0 = rho .* Vinf .* rN ./ mu_T0;




    
%P00
for i=1:hmax+1
P00(i)=atmo_p(i, Ti(i), sum(rhoi(i,:)));
end
P00=P00';

%% graphs
figure(1)
subplot(2,2,1)
plot(Ti,Z,'LineWidth',2.5)
hold on
plot(T,Z,'r','LineWidth',2.5)
ylabel('altitude [Km]')
xlabel('Temperature[K]')
legend('NRLMSISE-00','US 76')
hold off

subplot(2,2,2)
semilogx(P00,Z,'b','LineWidth',2.5)
hold on
semilogx(P,Z,'r','LineWidth',2.0)
ylabel('altitude [Km]')
xlabel('Pressure[Pa]')


subplot(2,2,3)
semilogx(rho00,Z,'b','LineWidth',2.5)
hold on
semilogx(rho,Z,'r','LineWidth',2.0) 
ylabel('altitude [Km]')
xlabel('Air density [kg/m^3]')
% legend('NRLMSISE-00','US 76')


subplot(2,2,4)
plot(mu00,Z,'b','LineWidth',2.5)
hold on
plot(mu,Z,'r','LineWidth',2.5)
ylabel('altitude [Km]')
xlabel('Dynamic Viscosity [Pa*s]')
% legend('NRLMSISE-00','US 76')

figure(2)
subplot(2,2,1)
plot(c,Z,'LineWidth',2.5)
hold on
plot([c(1:86); c76],[(1:86)'; (87:301)'] ,'r','LineWidth',2.5)

ylabel('altitude [Km]')
xlabel('c [m/s]')
legend('NRLMSISE-00','US 76')
hold off

subplot(2,2,2)
semilogx(Kn,Z,'LineWidth',2.5)
hold on
semilogx(Kn76,Z,'r','LineWidth',2.0)
ylabel('altitude [Km]')
% legend('NRLMSISE-00','US 76')
xlabel('Kn')


subplot(2,2,3)
plot(Pr,Z,'LineWidth',2.5)
hold on
plot(Pr76,Z,'r','LineWidth',2.5)
ylabel('altitude [Km]')
xlabel('Pr')
% legend('NRLMSISE-00','US 76')

subplot(2,2,4)
plot(m,Z,'LineWidth',2.5)
hold on
plot(m76,Z,'r','LineWidth',2.5)
ylabel('altitude [Km]')
xlabel('gravity [m/s^2]')


print(figure(1),'atmosphere1','-dmeta')
print(figure(2),'atmosphere2','-dmeta')