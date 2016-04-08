%JSK

clear all
clc

tic

%%%%% Read in the STL file for the Geometry %%%%%%%%%

[V, F, N] = import_stl_fast('sphere1_ascii.stl',1);
b = 33;            %Characteristic Length of Object

%%%%% Truncating for unwanted rows %%%%%%%
if length(N) > length(F) 
    N = N(1:end-1,:);    
end

%%%%%% Calculating Elemental Areas %%%%%%%%%
A = Area(F,V);
tr = triangulation(F,V);
P = incenter(tr);
CG = COG(P);

%%%%%% Cleaning up duplicates %%%%%%%%
[ N, F, A ] = cleanup( N, F, A );

%%%%%% Modified Newtonian Theory Pressure Calculations %%%%%%%%
%% For continumm regime %%%

%%%%%% User Inputs %%%%%%%

alpha = 10;     % Angle of attack, degrees
beta = 0;       % Yaw angle, degrees
Vinf = 300*2;     % Free Stream Velocity, m/s
Vinf = [Vinf * cosd(alpha) * cosd(beta), - Vinf * sind(beta), Vinf * sind(alpha) * cosd(beta)];
% %Free Stream Velocity Vector
% 
SN = 1.0;   %Normal Momentum Accommodation coefficient
ST = 1.0;   %Tangential Momentum Accommodation Coefficient
AT = ST*(2-ST);     %Tangental Energy Accomodation Coefficient
AN = 1 - (1 - SN)^2;    %Normal Energy Accommodation Coefficient
AC = 0.5 * (AN + AT);   %Overall Energy Accommodation Coefficient
Pr = 0.71;      %Prandtl Number

Tinf = 250;
Tw = Tinf*(1+0.2*(norm(Vinf)/300)^2); %Surface Temperature
% 
% %%%% Modified Newtonian Routine
Out = MaxPresSA(Vinf,Tinf,Tw);
Cpmax = Out(1);
Minf = Out(2);
s = Out(3);
L_b = Out(4);
rho = Out(5);
g = Out(6);
mu_T0 = Out(7);
T0 = Out(8);
p = Out(9);         %Pressure
Ps = Out(10);       %Free Stream Stagnation Pressure
rhos = Out(11);     %Free Stream Stagnation Density
mu_Tinf = Out(12);  %Free Stream Viscosity
mu_w = Out(13);     %Wall Viscosity
rhow = Out(14);     %Density at the Wall
hw = Out(15);       %Wall enthalpy
h0 = Out(16);       %stagnation enthalpy

% 
% %%%%% Calculating Force and Moment Coefficients %%%%%
% 
Kn = L_b/b;      %Knudsen Number
% 
% %%% Parameter Initialization
% 
Cpc = zeros(length(N),1);
Cc = zeros(1,3);
Cfm = zeros(1,3);
Mc = zeros(1,3);
Mfm = zeros(1,3);

dCxc = 0;
dCyc = 0;
dCzc = 0;
dCxfm = 0;
dCyfm = 0;
dCzfm = 0;

Sref = 0;

Q = zeros(length(N),1);
Qfr = zeros(length(N),1);
Stc = zeros(length(N),1);
Stcfr = zeros(length(N),1);

Vinfn = Vinf/norm(Vinf);
rN = 1;
% cp = 1373;   ///////////////////////////////////////////////////////////
cp = 1004.7;

dudx = (1/rN) * sqrt(2*(Ps-p)/rhos);         %Fay and Ridell Parameter

if Kn >= 10
    Qstfm = 0.5 * AC * rho * norm(Vinf)^3;
else
    Re0 = rho * norm(Vinf) * rN / mu_T0;
    Stfm = 1;
    Qsfr = 0.76*(Pr^-0.6)*((rhos*mu_T0)^0.4)*((rhow*mu_w)^0.1)*sqrt(dudx)*(h0-hw);
%     Qsfr = 0.94*((rho*mu_Tinf)^0.4)*((rhow*mu_w)^0.1)*sqrt(dudx)*(h0-hw);

end

% Cpmax = 2.2;////////////////////////////////////////////////////////////
for i = 1:length(N)
       
    theta = acosd(dot(-N(i,:),Vinf)/norm(-N(i,:))/norm(Vinf));
    delta = 90 - theta;     %Local Inclination Angle
    r = P(i,:) - CG;    %Moment Arm
    
    %%%%%%%%%%%%% CONTINUUM AERODYNAMICS %%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if delta > 0
    Cpc(i,1) = Cpmax * (sind(delta))^2;
    Sref = Sref + A(i,1) * sind(delta);
    end
    
    Cpnc = Cpc(i,1) * -N(i,:) * A(i,1);
    Cc = Cc + Cpnc;
    Mc = Mc + cross(r,Cpnc);
    
    %%%%%%%%%%% FREE MOLECULAR AERODYNAMICS %%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    theta = asind (dot(-Vinfn,N(i,:)));
    
    t = [0,0,0];
    Cpfm = 0;
    Ctfm = 0;
    
    if theta == 0
        t = (N(i,:) * dot(Vinfn,N(i,:)) - Vinfn) / sqrt(1 - dot(Vinfn,N(i,:))^2);        
        Ctfm = -(ST*cosd(theta)/s/sqrt(pi)) * (exp(-(s*sind(theta))^2) + sqrt(pi) * s *sind(theta)*(1+erf(s*sind(theta))));
    elseif theta == 90
        Cpfm1 = exp(-(s*sind(theta))^2) * ((2-SN)*s*sind(theta)/sqrt(pi) + SN*sqrt(Tw/Tinf)/2);
        Cpfm2 = (1+erf(s*sind(theta))) * ((2-SN) * ((s*sind(theta))^2 + 0.5) + 0.5*SN*s*sind(theta)*sqrt(pi*Tw/Tinf));
        Cpfm = (Cpfm1+Cpfm2)/s^2;
    elseif theta ~= 0 && theta ~= 90 && theta ~= -90
        t = (N(i,:) * dot(Vinfn,N(i,:)) - Vinfn) / sqrt(1 - dot(Vinfn,N(i,:))^2);
        Cpfm1 = exp(-(s*sind(theta))^2) * ((2-SN)*s*sind(theta)/sqrt(pi) + SN*sqrt(Tw/Tinf)/2);
        Cpfm2 = (1+erf(s*sind(theta))) * ((2-SN) * ((s*sind(theta))^2 + 0.5) + 0.5*SN*s*sind(theta)*sqrt(pi*Tw/Tinf));
        Cpfm = (Cpfm1+Cpfm2)/s^2;
        Ctfm = -(ST*cosd(theta)/s/sqrt(pi)) * (exp(-(s*sind(theta))^2) + sqrt(pi) * s *sind(theta)*(1+erf(s*sind(theta))));
    end
    
    Cpnfm = Cpfm * -N(i,:) * A(i,1);
    Cttfm = Ctfm * t * A(i,1);
    Cfm = Cfm + Cpnfm + Cttfm;
    Mfm = Mfm + cross(r,Cpnfm+Cttfm);
%     
%     
%     %%%%%%%%%%% CONTINUUM AERODYNAMICS HEATING %%%%%%%%%%%
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     
    if Kn < 10
        Stsc = 2.1/Re0;        %Stagnation Point Stanton Number, SCARAB Formulation
        Stsfr = Qsfr/rho/norm(Vinf)/(h0-hw);
    end
%     
%     %%%%%%%%%%% TRANSLATIONAL HEATING %%%%%%%%%%%%%%%%%%%%
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     
    if Kn > 1e-2 && Kn < 10
        Stsc = Stsc/sqrt(1+(Stsc/Stfm)^2);      %Translational Stagnation Point Stanton Number, SCARAB Bridging Function
%         Stsc = Stsc/sqrt(1+Stsc/Stfm);         %Translational Stagnation Point Stanton Number, SESAM Bridging Function
        Stsfr = Stsfr/sqrt(1+(Stsfr/Stfm)^2);      %Translational Stagnation Point Stanton Number, SCARAB Bridging Function
%         Stsfr = Stsfr/sqrt(1+Stsfr/Stfm);         %Translational Stagnation Point Stanton Number, SESAM Bridging Function   
    end
    
    if Kn < 10
    if theta >= 0 && theta <= 90
        Stc(i,1) = Stsc * (0.7*sind(delta)); 
        Stcfr(i,1) = Stsfr * (0.1+0.9*sind(delta)); 
        Q(i,1) = Stc(i,1) * rho * norm(Vinf) * cp * (T0 - Tw);
        Qfr(i,1) = Stcfr(i,1) * rho * norm(Vinf) * cp * (T0 - Tw);
        Total_heat(i,1) = Stc(i,1)*A(i,1);
    end
    end
%     
%     %%%%%%%%%%%%%%%%% FREE MOLECULAR FLOW %%%%%%%%%%%%%%%
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     
    if Kn >= 10
        Q(i,1) = (Qstfm/(s^3)/2/sqrt(pi))*((s^2+g/(g-1)-((g+1)*Tw/2/(g-1)/Tinf))*(exp(-(s*sind(delta))^2)+sqrt(pi)*s*sind(delta)*(1+erf(s*sind(delta))))-0.5*exp(-(s*sind(delta))^2));
    end
%      
end



Cc = Cc/Sref;
Cfm = Cfm/Sref;

Mc = Mc/Sref/b;
Mfm = Mfm/Sref/b;

B2WA = [cosd(alpha)*cosd(beta), -sind(beta), sind(alpha)*cosd(beta); cosd(alpha)*sind(beta), cosd(beta), sind(alpha)*sind(beta);-sind(alpha), 0, cosd(alpha)];

CDc = dot(Cc, B2WA(1,:));
CSc = dot(Cc, B2WA(2,:));
CLc = dot(Cc, B2WA(3,:));

CDfm = dot(Cfm, B2WA(1,:));
CSfm = dot(Cfm, B2WA(2,:));
CLfm = dot(Cfm, B2WA(3,:));

Mlc = Mc(1);
Mmc = Mc(2);
Mnc = Mc(3);

Mlfm = Mfm(1);
Mmfm = Mfm(2);
Mnfm = Mfm(3);

Mlc = dot(Mc, B2WA(1,:));
Mmc = dot(Mc, B2WA(2,:));
Mnc = dot(Mc, B2WA(3,:));

Mlfm = dot(Mfm, B2WA(1,:));
Mmfm = dot(Mfm, B2WA(2,:));
Mnfm = dot(Mfm, B2WA(3,:));

% CD_sp_FM_SC = 0.5*(2-SN+ST)*((4*s^4+4*s^2-1)*erf(s)/2/s^4 + ((2*s^2+1)*exp(-s^2)/sqrt(pi())/s^3)) + 2*SN*sqrt(pi)*sqrt(Tw/Tinf)/3/s;

CDtrans = CDc + (CDfm - CDc) * (sin(pi*(0.5+0.25*log10(Kn))))^3
CStrans = CSc + (CSfm - CSc) * (sin(pi*(0.5+0.25*log10(Kn))))^3
CLtrans = CLc + (CLfm - CLc) * (sin(pi*(0.5+0.25*log10(Kn))))^3

Cltrans = Mlc + (Mlfm - Mlc) * (sin(pi*(0.5+0.25*log10(Kn))))^3
Cmtrans = Mmc + (Mmfm - Mmc) * (sin(pi*(0.5+0.25*log10(Kn))))^3
Cntrans = Mnc + (Mnfm - Mnc) * (sin(pi*(0.5+0.25*log10(Kn))))^3
% 
% 
figure(1)
trisurf(F,V(:,1),V(:,2),V(:,3),Cpc,'EdgeColor','none','LineStyle','none','FaceLighting','phong')
xlabel('X')
ylabel('Y')
zlabel('Z')
axis equal
colorbar

figure(2)
trisurf(F,V(:,1),V(:,2),V(:,3),Q,'EdgeColor','none','LineStyle','none','FaceLighting','phong')
xlabel('X')
ylabel('Y')
zlabel('Z')
axis equal
colorbar

figure(3)
trisurf(F,V(:,1),V(:,2),V(:,3),Qfr,'EdgeColor','none','LineStyle','none','FaceLighting','phong')
xlabel('X')
ylabel('Y')
zlabel('Z')
axis equal
colorbar

toc
