clc
clear all
%% User Inputs
h=500;          % Final Geometric Altitude[km]
alpha = 10;     % Angle of attack, [deg]
beta = 0;       % Yaw angle, [deg]
Vinf = 900;   % Free Stream Velocity, [m/s]
Lat= 55.8642;   % Latitude
Lon= -4.16;     % Longitude
year= 2010;     % year
yr_day= 1;      % day
UT_sec= 0;      % UT seconds
L= 1300;        % lunghezza caratteristica oggetto
%% Atmosphere parameters

%               1976 Standard Atmosphere Calculator[0-1000 km]
%               Z:          Total Reporting Altitudes[0<=alt<=1000 km][km]{ft}
%               Z_L:        Lower Atmosphere Reporting Altitudes[0<=alt<=86 km][km]{ft}
%               Z_U:        Upper Atmosphere Reporting Altitudes[86<=alt<=1000 km][km]{ft}
%               T:          Temperature array[0<=alt<=1000 km][K]{R}
%               P:          Pressure array[0<=alt<=1000 km][Pa]{in_Hg}
%               rho:        Density array[0<=alt<=1000 km][kg/m^3]{lb/ft^3}
%               c:          Speed of sound array[0<=alt<=86 km][m/s]{ft/s}
%               g:          Gravity array[0<=alt<=1000 km][m/s^2]{ft/s^2}
%               mu:         Dynamic Viscosity array[0<=alt<=86 km][N*s/m^2]{lb/(ft*s)}
%               nu:         Kinematic Viscosity array[0<=alt<=86 km][m^2/s]{ft^2/s}
%               k:          Coefficient of Thermal Conductivity
%                           array[0<=alt<=86 km][W/(m*K)]{BTU/(ft*s*R)}
%               n:          Number Density of individual gases
%                           (N2 O O2 Ar He H)[86km<=alt<=1000km][1/m^3]{1/ft^3}
%               n_sum:      Number Density of total gases
%                           [86km<=alt<=1000km][1/m^3]{1/ft^3}

[Z, Z_L, Z_U, T, P, rho, g, mu, nu, k, n] = atmo(h,1,1);

% atmospheric composition  
%ref. 2001 United States Naval Research Laboratory Mass Spectrometer and Incoherent Scatter Radar Exosphere 
%see http://uk.mathworks.com/help/aerotbx/ug/atmosnrlmsise00.html for description

[Ti,rhoi] = atmosnrlmsise00(0:1000:(h)*10^3, Lat, Lon, year, yr_day, UT_sec);
rhoi(:,[6,9]) = [ ]; %He O N2 O2 Ar H N [1/m^3]
Ti(:,1)= [];

% heat capacity ratio
for i=1:1:length(rhoi-1)

percentual(i,:)= rhoi(i,1:7)./sum(rhoi(i,1:7)); %He O N2 O2 Ar H N

gamma(i)= percentual(i,:)* [ 5/3 5/3 1.4 1.4 5/3 5/3 5/3]';

end

figure
plot(gamma,Z)
title('gamma vs altitude')

%sound speed 
c1=sqrt((gamma'.*P)./rho);

figure
plot(c1,Z)
title('speed sound vs altitude')

%dynamic Viscosity
if h<=86;
    mu=mu;
else 
    for i=87:h+1  %costant values above 86 km   
    mu(i)=mu(87); 
    end
end
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

%% Free Stream values

gamma=gamma(h-1); %heat capacity ratio% 
Vinf=[Vinf * cosd(alpha) * cosd(beta), - Vinf * sind(beta), Vinf * sind(alpha) * cosd(beta)]; % %Free Stream Velocity Vector
normVinf=norm(Vinf); %norm Vinf
Minf=normVinf./c1(h+1); %Mach number
Tinf=Ti(h+1); %free stream temperature
pinf=P(h+1); %free stream pressure
ptot1=pinf*((1+(gamma-1)/2)*Minf^2)^(gamma/(gamma-1)); %free stream stagn. pressure
A=[(1-gamma+2*gamma*Minf^2)/(gamma+1)]*[((gamma+1)^2*Minf^2)/(4*gamma*Minf^2-2*(gamma-1))]^(gamma/(gamma-1));
ptot2=pinf*A;    % stagnation pressure behind a normal shock wave 

% Knudsen Number evaluation

Kb=1.38064852e-23;                                  %Boltzmann constant [J/Kg]
Nav = 6.023e23;                                     %Avogadro Constant 
f=1.660539040e-27;                                  %conversion factor atomic mass unit - kg
l=(mu(h+1)/pinf)*sqrt(pi*Kb*Tinf/(2*(m(h+1)*f)));   %mean free path [m]
Kn = l/L;                                           %Knudsen Number

%% graphs
figure
plot(T,Z)
hold on 
plot(Ti,Z,'r')
legend('T 1976 model','T 2001 model');
hold off

figure
semilogx(P,Z)
title('pressure vs altitude')

figure
semilogx(rho,Z)
title('air density vs altitude')

figure
plot(mu,Z)
title('dynamic viscosity vs altitude')

figure
plot(xN2,Z)
title('N2% vs altitude')

figure
plot(xO2,Z)
title('O2% vs altitude')

figure
semilogx(xO,Z)
title('O% vs altitude')

figure
semilogx(xN,Z)
title('N% vs altitude')

figure
plot(xAr,Z)
title('Ar% vs altitude')

figure
semilogx(xHe,Z)
title('He% vs altitude')

figure
semilogx(xH,Z)
title('H% vs altitude')

