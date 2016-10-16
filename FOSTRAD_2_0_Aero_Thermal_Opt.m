function F = FOSTRAD_2_0_Aero_Thermal_Opt(mainDir, STLname, altitude, Vinf, lref, Sref, alpha, beta, Nsmooth, Bfc, AeroFlag, ThermoFlag)
%  F = FOSTRAD_2_0_Aero_Thermal_Opt(mainDir, STLname, altitude, Vinf, lref, Sref, alpha, beta, Nsmooth, Bfc, AeroFlag, ThermoFlag)
% inputs: description can be found in the FOSTRAD_Controller
% outputs:  F = [h(Hindex(1)) CD_LP5(Hindex(1)) CL_LP5(Hindex(1)) Qav(Hindex(1))];
% where Hindex is the index of the input altitude
% mainDir = '/mnt/lustre/strath/aeroeng/wyb15207/FOSTRAD 2.0/';
% STLname = 'Ellipsoid_coarse_bin.STL';
% altitude = 85; [km]  min~max = 0~999km
% lref = 4;      [m]
% Sref = pi()*1.15^2 [m^2] reference cross section for aerodynamics
% alpha = 20;    [deg]      Pitch Angle
% beta = 5;      [deg]      Yaw Angle
% Nsmooth = 4;   %number of smoothing point pre and post selected altitude%
% % Bfc = BackFaceCulling_opt()  %output from backface culling algorithm
% AeroFlag = 1;   % if == 1 the aerodynamics module is active
% ThermoFlag = 1;   % if == 1 the aero-thermodynamics module is active


warning('off')
tstart = tic;

%Additional required folders
cd(mainDir);
addpath(mainDir);
CadFolder=[mainDir, 'CAD MODELS'];
addpath(CadFolder)
addpath([mainDir,'external functions']);


%% %%%%%%%********************  FOSTRAD  V. 2.0 optimized mod **********************%%%%%%%%%


% ******************* Additional User Inputs ********************
rN =1;     %Nose radius (m)
hmin=altitude; %km              %minimum altitude present in the database: 0km
hmax=altitude; %km              %maximum altitude present in the database: 999km
nHsteps = 1;
h = linspace(hmin,hmax,nHsteps)';  %real vector that will be used for interpolation (Query points)
Hcont = 40;             %Fixed reference continuum altitude
Hfmf = 220;             %Fixed reference Free molecular flow altitude
limKn_inf=1E-4;         %fixed reference Kn continuum limit
limKn_sup=10;           %fixed reference Kn FMF limit
Twi=350;                %fixed wall temperature

% Atmospheric model and GSI constant inputs
SN = 1.0;               %Normal Momentum Accommodation coefficient
ST = 1.0;               %Tangential Momentum Accommodation Coefficient
AT = ST*(2-ST);         %Tangental Energy Accomodation Coefficient
AN = 1 - (1 - SN)^2;    %Normal Energy Accommodation Coefficient
AC = 0.5 * (AN + AT);   %Overall Energy Accommodation Coefficient
cp = 1004.7;            %specific heat at constant pressure (AIR) [J/kgK]

% NRLMSISE-00 DATABASE condition for:
% Year= 2000, Month=  1, Day=  1, Hour= 1.50,
% Time_type = Universal
% Coordinate_type = Geographic
% Latitude=   55.00deg, Longitude=   45.00deg,
% Database format [1000x9]
% Selected output parameters:
% 1  Height, km
% 2  O, cm-3
% 3  N2, cm-3
% 4  O2, cm-3
% 5  Temperature_neutral, K
% 6  He, cm-3 
% 7  Ar, cm-3
% 8  H, cm-3
% 9  N, cm-3
load('NRLMSISE00')
%A database with different condition may be generated and download from:
%     http://ccmc.gsfc.nasa.gov/modelweb/models/nrlmsise00.php


%% Starting Checks

% checks required to guarantee the routine stability
if Nsmooth < 3 && nHsteps < 3
    disp('%%********* WARNING *********%%')
    error('%%********* minimum Nsmooth = 3;  *********%%')
end



if hmin < 0 || hmax > 999
    disp('%%********* WARNING *********%%')
    error('%%********* Altitudes exceeding atmospheric models Database boundaries [0 ~999]km *********%%')
end

if Hcont(1,1) > hmin
    Hcont(1,1) = hmin;
end
if Hfmf(1,1) < hmax
    Hfmf(1,1) = hmax+1;
    if Hfmf(1,1) > 999
        Hfmf(1,1) = 999;
    end
end

% defining the altitude spacing among additional smoothing points
Spacing = (1-(exp(linspace(log(0.001+0.6/Nsmooth),log(1),Nsmooth)))');
Nsmooth = round(Nsmooth);


if nHsteps < 3
    if  Hcont ~= hmin;
        if  hmin >= 140
            Hcont = linspace(Hcont,(Hfmf/2),Nsmooth+1)';
            Hfmf = linspace((Hfmf/2),Hfmf,Nsmooth+1)';
            h = [Hcont; hmin; Hfmf(2:end)];
            h = sort(h);
            Hindex = find(h == hmin);
        else
            Hcont = Hcont + (hmin-Hcont)*flip(Spacing);
            Hfmf = Hfmf - (Hfmf-hmax).*Spacing;
            h = [Hcont; h; Hfmf];
            h = sort(h);
            Hindex = find(h == hmin);
        end
    elseif Hcont == hmin
        Hindex = 1;
        Hcont = linspace(Hcont,(Hfmf/2),Nsmooth+1)';
        Hfmf = linspace((Hfmf/2),Hfmf,Nsmooth+1)';
        h = [Hcont; Hfmf(2:end)];
    end
end

%% %%% Read in the STL file for the Geometry %%%%%%%%%
% [V, F, N] = import_stl_fast(STLname,1);       % ASCII stl
[F,V,N] = stlread(STLname);                     % Binary stl

%%%%% Truncating for unwanted rows %%%%%%%
if length(N) > length(F)
    N = N(1:end-1,:);
end
%%%%%% Calculating Elemental Areas %%%%%%%%%
A= Area(F,V);
Atot = sum(A);
tr = triangulation(F,V);
In = incenter(tr);
CG = COG(In);

%% atmosphere calculation
jj=1;       %initialising indexes
bb=1;



% Rearranging the database inputs to match FOSTRAD inputs
%He O N2 O2 Ar H N [1/m^3] gas densities
rhoi = [NRLMSISE00((Hcont(1,1)+1):(Hfmf(end)+1),6) NRLMSISE00((Hcont(1,1)+1):(Hfmf(end)+1),2:4) NRLMSISE00((Hcont(1,1)+1):(Hfmf(end)+1),7:9)]*1E6;
Ti = NRLMSISE00((Hcont(1,1)+1):(Hfmf(end)+1),5);


%inputs corresponding to the input altitudes interpolated with the
%continuum and FMF reference
hDB = (Hcont(1):1:Hfmf(end))'; %altitudes used for the sample points interpolation on the database
if length(hDB) > 1
    Ti = interp1(hDB,Ti(:,1),h,'spline');
    for i=1:7
        rhotemp(:,i) = interp1(hDB,rhoi(:,i),h,'spline');
    end
    rhoi = rhotemp;
end

percentual= rhoi(:,1:7)./  ((rhoi(:,1:7)* ones(7,1)) * ones (1,7)) ; % percentages He O N2 O2 Ar H N
gamma= percentual* [ 5/3 5/3 1.4 1.4 5/3 5/3 5/3]';

R = cp * (gamma-1) ./ gamma;     %Gas constant, kg m2 s-2 K-1 mol-1
omega = 0.78;           %Temperature Viscosity power law coefficient
C1 = 1.458e-6;          %Sutherland's Law Coefficient, kg/m-s-K^-0.5
C2 = 110.4;             %Sutherland's Law Coefficient, K

%molecular weights
mN2 = 28.01340;               %molar mass of nitrogen molecule, grams/mole
mO2 = 31.99880;               %molar mass of oxigen molecule, grams/mole
mO = mO2/2;                   %molar mass of oxigen atom, grams/mole
mN = mN2/2;                   %molar mass of Nitrogen atom, grams/mole
mAr = 39.9480;                %molar mass of Argon molecule, grams/mole
mHe = 4.0026020;              %molar mass of helium molecule, grams/mole
mH= 1.007940;                 %molar mass of Hydrogen molecule, grams/mole

%fraction disitribution
xN2=percentual(:,3);
xO2=percentual(:,4);
xO= percentual(:,2);
xN= percentual(:,7);
xAr= percentual(:,5);
xHe= percentual(:,1);
xH= percentual(:,6);

m=mN2.*xN2+mO2.*xO2+mO.*xO+mN.*xN+mAr.*xAr+mHe.*xHe+mH.*xH; %Atmosphere molecular weight
rho = rhoi*[mHe mO mN2 mO2 mAr mH mN]'./6.022169e26;
P =  rho .* R .* Ti;

amu = 1.66053892e-27;   % FMF properties %
mg = m.*amu;
Kb = 1.3806488e-23;     %Boltzmann constant, m2 kg s-2 K-1
Nav = 6.023e23;         %Avogadro Constant
f=1.660539040e-27;      %conversion factor atomic mass unit - kg
k = 2.64638e-3*Ti.^1.5./(Ti+245*10.^(-12./Ti));     %   Thermal Conductivity Coefficient
mu = 1.458e-6*Ti.^1.5./(Ti+110.4);                  % dynamic Viscosity
c=sqrt((gamma.*P)./rho);                             %sound speed

Minf = norm(Vinf) ./ c;         %Mach

% Free Stream Stagnation Pressure and Temperature
P01 = P .* (1 + 0.5 .* (gamma-1) .* Minf.^2) .^ (gamma./(gamma-1));
T0 = Ti .* (1 + 0.5 .* (gamma-1) .* Minf.^2);

Tw = Ti.*(1+(gamma-1)./(2*(Minf).^2)); %Surface Temperature (original formula)

mu_T0 = mu.* (T0./Ti).^omega;
rhos = rho.*(T0./Ti).^(1./(gamma-1));
rhow = P./R./Tw;  %Assuming constant pressure in the boundary layer
mu_w = (C1 .* Tw.^(3/2))./(Tw + C2);
hw = cp .* Tw;
h0 = cp .* T0;

% Stagnation Pressure after Normal Shock
P02 = P01 .* (((gamma+1).*Minf.^2)./((gamma-1).*Minf.^2+2)).^(gamma./(gamma-1)) .* ...
    ((gamma+1)./(2.*gamma.*Minf.^2-(gamma-1))).^(1./(gamma-1));

Cpmax=(2./(gamma.*Minf.^2)).*((P02./P-1)); %Cp maximum

% FMF properties %
mg = m.*amu;
kb = 1.3806488e-23;      %Boltzmann constant, m2 kg s-2 K-1
vmp = sqrt(2.*kb.*Ti./mg);
s = Vinf./vmp;

l=(mu./P).*sqrt(pi.*kb.*Ti./(2.*mg));  %mean free path [m]

%Knudsen Number
Kn = l./lref;
Pr=mu.*cp./k; %Prandtl number

%% Heat transfer models

dudx = (1./rN) .* sqrt(2.*(P01-P)./rhos);         %Fay and Ridell Parameter
Qstfm = 0.5 * AC .* rho .* Vinf.^3;
Re0 = rho .* Vinf .* rN ./ mu_T0;
Stfm = 1;
Qsfr = 0.76.*(Pr.^-0.6).*((rhos.*mu_T0).^0.4).*((rhow.*mu_w).^0.1).*sqrt(dudx).*(h0-hw);
Qsfr1 = 0.94.*((rho.*mu).^0.4).*((rhow.*mu_w).^0.1).*sqrt(dudx).*(h0-hw); %fay riddel
Qsfr2=0.94.*((rho.*mu).^0.5).*sqrt(dudx).*(h0-hw); %Van Driest

%%  Parameter Initialization

%Loading the Backface Culling triangles
A = Bfc(:,1);
F = Bfc(:,2:4);
N = Bfc(:,5:7);
tr = triangulation(F,V);            %new triangulation after the BFC has been applied
In = incenter(tr);                  %new triangles centres matrix definition after the BFC has been applied

Cpc = zeros(length(N),1);

if ThermoFlag == 1;
    Stcfr_KDR = zeros(length(N),1);
    Stcfr_FOSTRAD20 =  zeros(length(N),1);
    Q_KDR = zeros(length(N),1);
    Q_FOSTRAD20 = zeros(length(N),1);
    Q_c = zeros(length(N),1);
    Total_heat_c = zeros(length(N),1);
    Q_fm = zeros(length(N),1);
    Stfm = zeros(length(N),1);
    Stfm1 = zeros(length(N),1);
    Total_heat_fm = zeros(length(N),1);
    Q = zeros(length(N),1);
    Qfr = zeros(length(N),1);
    Stc = zeros(length(N),1);
    Stcfr = zeros(length(N),1);
end

%free stream velocity vector
Vinfi = [Vinf * cosd(alpha) * cosd(beta), - Vinf * sind(beta), Vinf * sind(alpha) * cosd(beta)]; % %Free Stream Velocity Vector








%% Computation at different altitudes

%Body to Wind rotation matrix
B2WA = [cosd(alpha)*cosd(beta), -sind(beta), sind(alpha)*cosd(beta); ...
    cosd(alpha)*sind(beta), cosd(beta), sind(alpha)*sind(beta);-sind(alpha), 0, cosd(alpha)];

disp('Simulation started')

for H=1:length(h) %Km
    
    disp(['Altitude: ',num2str(h(H)),'km'])
    
    % Parameter Re-Initialization
    Cc = zeros(1,3);
    Cfm = zeros(1,3);
    Mc = zeros(1,3);
    Mfm = zeros(1,3);
%     Sref = 0;
    
    
    % Recalling variables that have been previously allocated
    Pi=P(H);
    %     Ri=R(H);
    gammai=gamma(H);
    %     ci=c(H);
    Minfi=Minf(H);
    %     P01i=P01(H);
    T0i=T0(H);
    %     mui = mu(H);
    %     mu_T0i = mu_T0(H);
    rhoi=rho(H);
    rhosi=rhos(H);
    %     rhowi=rhow(H);
    %     mu_wi=mu_w(H);
    hwi = hw(H);
    h0i = h0(H);
    %     P02i=P02(H);
    Cpmaxi=Cpmax(H);
    %     mi=m(H);
    %     mgi=mg(H);
    %     vmpi=vmp(H);
    si=s(H);
    %     li=l(H);
    Kni=Kn(H);
    %     Pri=Pr(H);
    Vinfni = Vinfi/Vinf;
    %     dudxi=dudx(H);
    Tii=Ti(H);
    Re0i=Re0(H);
    Qstfmi=Qstfm(H);
    if Kni >= 10
        Qstfmi = Qstfm(H);
    else
        Re0i = Re0(H);
        Qsfri = Qsfr(H);
    end
    
    %% Starting triangles computation cycling through all triangles
    for i = 1:length(N)
        
        theta = acosd(dot(-N(i,:),Vinfi)/norm(-N(i,:))/norm(Vinfi));
%         thetaSref = acosd(dot(-N(i,:),[1 0 0])/norm(-N(i,:)));  %Standard starting reference surface
        delta = 90 - theta;     %Local Inclination Angle
        r = In(i,:) - CG;    %Moment Arm
        
        
        %% ******************* AERODYNAMICS MODULE ********************** %%
        if AeroFlag ==1;
            %%%%%%%%%%%%% CONTINUUM AERODYNAMICS %%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            if delta > 0
                Cpc(i,1) = Cpmaxi * (sind(delta))^2;
%                 Sref = Sref + A(i,1) * sind(delta);    %reference surface computed with respect to the wind direction                      
%                 Sref(i,1) = abs(A(i,1) * sind(90 thetaSref));    % constant reference surface Sref(i,1)              
            end
            
            Cpnc = Cpc(i,1) * -N(i,:) * A(i,1);
            Cc = Cc + Cpnc;
            Mc = Mc + cross(r,Cpnc);
            
            %%%%%%%%%%% FREE MOLECULAR AERODYNAMICS %%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            thetaA = asind (dot(-Vinfni,N(i,:)));
            
            t = [0,0,0];
            Cpfm = 0;
            Ctfm = 0;
            
            if thetaA == 0
                t = (N(i,:) * dot(Vinfi,N(i,:)) - Vinfi) / sqrt(1 - dot(Vinfi,N(i,:))^2);
                Ctfm = -(ST*cosd(thetaA)/si/sqrt(pi)) * (exp(-(si*sind(thetaA))^2) +...
                    sqrt(pi) * si *sind(thetaA)*(1+erf(si*sind(thetaA))));
            elseif thetaA == 90
                Cpfm1 = exp(-(si*sind(thetaA))^2) * ((2-SN)*si*sind(thetaA)/sqrt(pi) +...
                    SN*sqrt(Twi/Tii)/2);
                Cpfm2 = (1+erf(si*sind(thetaA))) * ((2-SN) * ((si*sind(thetaA))^2 + 0.5) +...
                    0.5*SN*si*sind(thetaA)*sqrt(pi*Twi/Tii));
                Cpfm = (Cpfm1+Cpfm2)/si^2;
            elseif thetaA ~= 0 && thetaA ~= 90 && thetaA ~= -90
                t = (N(i,:) * dot(Vinfni,N(i,:)) - Vinfni) / sqrt(1 - dot(Vinfni,N(i,:))^2);
                Cpfm1 = exp(-(si*sind(thetaA))^2) * ((2-SN)*si*sind(thetaA)/sqrt(pi) +...
                    SN*sqrt(Twi/Tii)/2);
                Cpfm2 = (1+erf(si*sind(thetaA))) * ((2-SN) * ((si*sind(thetaA))^2 + 0.5) +...
                    0.5*SN*si*sind(thetaA)*sqrt(pi*Twi/Tii));
                Cpfm = (Cpfm1+Cpfm2)/si^2;
                Ctfm = -(ST*cosd(thetaA)/si/sqrt(pi)) * (exp(-(si*sind(thetaA))^2) + sqrt(pi) *...
                    si *sind(thetaA)*(1+erf(si*sind(thetaA))));
            end
            
            Cpnfm = Cpfm * -N(i,:) * A(i,1);
            Cttfm = Ctfm * t * A(i,1);
            Cfm = Cfm + Cpnfm + Cttfm;
            Mfm = Mfm + cross(r,Cpnfm+Cttfm);
        end %end of Aerodynamics module
        
        %% ********************** Aero-ThermoDynamics Module *********************%%
        if ThermoFlag == 1;
            %     %%%%%%%%%%% CONTINUUM AERODYNAMICS HEATING %%%%%%%%%%%
            %     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %
            % Assigning each triangle a stanton number
            if Kni < limKn_inf
                Stsc(i,1) = 2.1/sqrt(Re0i);             %Stagnation Point Stanton Number, SCARAB Formulation
                Stsfr = Qsfri/rhoi/Vinf/(h0i-hwi); %St Numb. Fay Riddel
                if theta >= 0 && theta <= 90
                    Stc(i,1) = Stsc(i,1) * (0.7*sind(delta)); %da dove spunta sta formula?????
                    Stcfr(i,1) = Stsc(i,1) * (0.1+0.9*cosd(theta)); %modified Lees (SCARAB)
                    Stcfr_KDR(i,1) = Stsc(i,1) * (cosd(theta/2)^5.27); % Kemp Rose Detra
                    Stcfr_FOSTRAD20(i,bb) = Stsc(i,1)*0.74/2.1 * (0.1+0.9*cosd(theta)); % FOSTRAD2.0
                    if Stcfr_FOSTRAD20(i,bb)<0
                        Stcfr_FOSTRAD20(i,bb)=0;
                    end
                else
                    Stc(i,1)=0;
                    Stcfr(i,1)=0;
                    Stcfr_KDR(i,1)=0;
                    Stcfr_FOSTRAD20(i,bb)=0;
                end
                %computing continuum heat transfer (Q [W/m^2]) for different models
                Q(i,1) = Stc(i,1) * rhoi * Vinf* cp * (T0i - Twi);       % Scarab formulation [W/m^2]
                Qfr(i,1) = Stcfr(i,1) * rhoi * Vinf * cp * (T0i - Twi);
                Q_KDR(i,1)=Stcfr_KDR(i,1)* rhoi * Vinf * cp * (T0i - Twi);
                Q_FOSTRAD20(i,1)=Stcfr_FOSTRAD20(i,bb)* rhoi * Vinf * cp * (T0i - Twi);
                Q_c(i,1)=Q_FOSTRAD20(i,1);                               % Heat Transfer [W/m^2] continuum
                Total_heat_c(i,1) = Q_c(i,1)*A(i,1);                     %flusso termico (J/sec)
                Chc=2*Q(i,1)/(rhoi*Vinf^3); %heat transfer coefficient
            end
            % computing the averaged heat flow on the total initial area (pre-backface culling)
            
            if i == length(N) && Kni < limKn_inf
                Qav_c(H) = sum(Total_heat_c(:,1))/Atot; %[W/m^2]
            end
            
            %     %%%%%%%%%%%%%%%%% FMF THERMODYNAMICS HEATING %%%%%%%%%%%%%%%
            %     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %
            if Kni >= limKn_sup
                Stsfm= Qstfmi/rhoi/Vinf/(h0i-hwi);
                %       if theta >= 0 && theta <= 90
                Q_fm(i,1) = (Qstfmi/(si^3)/2/sqrt(pi))*((si^2+gammai/(gammai-1)-((gammai+1)*Twi/2/(gammai-1)/Tii))*...
                    (exp(-(si*sind(delta))^2)+sqrt(pi)*si*sind(delta)*(1+erf(si*sind(delta))))-0.5*exp(-(si*sind(delta))^2));
                %kempra fay riddel,
                Stfm(i,1)=Q_fm(i,1)/(rhoi * Vinf * cp * (T0i - Twi)); %heat transfer coefficient fm
                Stfm1(i,H)=Q_fm(i,1)/((rhoi/2)*Vinf^3);
                Total_heat_fm(i,1) = Q_fm(i,1)*A(i,1);
            end
            if i == length(N)
                Qav_fm(H) = sum(Total_heat_fm(:,1))/Atot; %[W/m^2]
            end
            
        end %End of ThermoDynamics Module
        
        
        
    end %end of triangles analyses
    bb=bb+1; % Altitude indexing for the transitional heating
    
    %%  Computation of Aerodynamics bridging parameters
    if AeroFlag == 1;
        Cc = Cc/Sref;
        Cfm = Cfm/Sref;
        Mc = Mc/Sref/lref;
        Mfm = Mfm/Sref/lref;
        CDc = dot(Cc, B2WA(1,:));
        CSc = dot(Cc, B2WA(2,:));
        CLc = dot(Cc, B2WA(3,:));
        CDfm = dot(Cfm, B2WA(1,:));
        CSfm = dot(Cfm, B2WA(2,:));
        CLfm = dot(Cfm, B2WA(3,:));
        Mlc = dot(Mc, B2WA(1,:));
        Mmc = dot(Mc, B2WA(2,:));
        Mnc = dot(Mc, B2WA(3,:));
        Mlfm = dot(Mfm, B2WA(1,:));
        Mmfm = dot(Mfm, B2WA(2,:));
        Mnfm = dot(Mfm, B2WA(3,:));
        
        %STS Bridging mod
        Cltrans = Mlc + (Mlfm - Mlc) * (sin(pi*(3.5/8+0.125*log10(Kni))))^2;
        Cmtrans = Mmc + (Mmfm - Mmc) * (sin(pi*(3.5/8+0.125*log10(Kni))))^2;
        Cntrans = Mnc + (Mnfm - Mnc) * (sin(pi*(3.5/8+0.125*log10(Kni))))^2;
        CDtrans = CDc + (CDfm - CDc) * (sin(pi*(3.5/8+0.125*log10(Kni))))^2;
        CStrans = CSc + (CSfm - CSc) * (sin(pi*(3.5/8+0.125*log10(Kni))))^2;
        CLtrans = CLc + (CLfm - CLc) * (sin(pi*(3.5/8+0.125*log10(Kni))))^2;
        
        % Matrixes definition [Cx, Kn(H)] for the continuum, FMF, transitional
        if Kn(H)>=limKn_sup
            CDfm_memo(H)=[CDfm(CDfm~=0)'];
            CLfm_memo(H)=[CLfm(CLfm~=0)'];
            CMfm_memo(H)=[Mmfm(Mmfm~=0)'];
            CSfm_memo(H)=[CSfm(CSfm~=0)'];
            if ThermoFlag == 0;
                Knfm(H)=Kn(H);
            end
        elseif Kn(H)<=limKn_inf
            CDc_memo(H)=[CDc(CDc~=0)'];
            CLc_memo(H)=[CLc(CLc~=0)'];
            CMc_memo(H)=[Mmc(Mmc~=0)'];
            CSc_memo(H)=[CSc(CSc~=0)'];
            if ThermoFlag == 0;
                Knc(H)=Kn(H)';
            end
        elseif Kn(H)>=limKn_inf && Kn(H)<=limKn_sup
            CDtrans_memo(H)=[CDtrans(CDtrans~=0)'];
            CLtrans_memo(H)=[CLtrans(CLtrans~=0)'];
            CMtrans_memo(H)=[Cmtrans(Cmtrans~=0)'];
            CStrans_memo(H)=[CStrans(CStrans~=0)'];
            if ThermoFlag == 0;
                Kntrans(H)=Kn(H)';
            end
        end
        
    end
end % altitude analyses ended
close all




%%     %%%%%%%%%%% TRANSLATIONAL HEATING %%%%%%%%%%%%%%%%%%%%
%      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ThermoFlag == 1;
    % indexing for transitional heating
    a=1;
    bb=1;
    c=1;
    for H=1:length(h) %it is necessary to compute the transitional heating after the FM and Continuum have been computed
        if Kn(H) >= limKn_inf && Kn(H) <= limKn_sup
            gr_trans(jj,:)=[h(H) Kn(H) rho(H) T0(H)]; %grandezze fisiche campo transition
            for i=1:length(N)
                theta = acosd(dot(-N(i,:),Vinfi)/norm(-N(i,:))/norm(Vinfi));
                delta = 90 - theta;     %Local Inclination Angle
                r = In(i,:) - CG;       %Moment Arm
                %if theta >= 0 && theta <= 90
                %Stctrans(i,1)= (Stcfr(i,1)+Kn(H)*Stfm1(i,end))/(1+Kn(H));
                Stctrans(i,1)=(Stcfr_FOSTRAD20(i,end)+Kn(H)*max(Stfm1(i,:)))/(1+Kn(H));
                Qtrans(i,1) = Stctrans(i,1) * rho(H) * Vinf * cp * (T0(H) - Twi);
                % Qfr(i,1) = Sttrans(i,1) * rho(H) * Vinf * cp * (T0i - Twi);
                Total_heat(i,1) = Qtrans(i,1)*A(i,1);
                %end
                if i==length(N)
                    Qav_t(H) = sum(Total_heat(:,1))/Atot;
                end
            end
            jj=jj+1;
        end
        
        if Kn(H)>=limKn_sup
            Knfm(a,1)=Kn(H);
            rho_fm(a)=rho(H);
            Stfm1_memo(a)= max(Stfm1(:,H));
            a=a+1;
        elseif Kn(H)<=limKn_inf
            Knc(bb,1)=Kn(H);
            rho_c(bb)=rho(H);
            Stcfr_memo(bb)= max(Stcfr_FOSTRAD20(:,bb));
            bb=bb+1;
        elseif Kn(H)>=limKn_inf && Kn(H)<=limKn_sup
            Kntrans(c,1)=Kn(H);
            Stctranns_memo(c)= max(Stctrans);
            Kn_trans(c,1)=Kn(H);
            c=c+1;
        end
    end
    
    % rearranging the matrixes for the HEATING module


%%    
    Qav_c = Qav_c(Qav_c ~= 0);
    Qav_fm = Qav_fm(Qav_fm ~= 0);
    Qav_t = Qav_t(Qav_t ~= 0);
    Qav = [Qav_c Qav_t Qav_fm]';
    Stcfr_memo=[Stcfr_memo(Stcfr_memo~=0)'];
    Stctranns_memo=[Stctranns_memo(Stctranns_memo~=0)'];
    Stfm1_memo=[Stfm1_memo(Stfm1_memo~=0)'];
    rho_c=rho_c(rho_c~=0)';
    rho_fm=rho_fm(rho_fm~=0)';
    
    %     %% Creating the bridging functions for the transitional heating
    glob_ref=[0.0014 0.0080 0.01 0.0224 0.028 0.0754  0.098 0.1 0.2044 0.227 0.476  0.8826 1 1.219  2.909  3.9908  ; 0.0224 0.0696 0.2615 0.170 0.176 0.4 0.405 0.573 0.665 0.595 0.687 0.802 0.9346 0.838 0.874 0.914 ]; %STS DATA (BELOW) AND SPHERE DATA
    
    [cf1_th]=L5P([Knc; glob_ref(1,:)'; Knfm ],[Stcfr_memo; glob_ref(2,:)'; Stfm1_memo]);
    
    cf_vls_St=coeffvalues(cf1_th); %coefficients p(x) bridging
    
    St_trans_LP5=Stfm1_memo(end) + (Stcfr_memo(1) - Stfm1_memo(end))./((1+(Kn_trans./cf_vls_St(3)).^cf_vls_St(2)).^cf_vls_St(5));
    
    [cf2_th]=L5P([Knc; Kn_trans; Knfm],[Stcfr_memo; St_trans_LP5; Stfm1_memo]);
    cf_vls_St_glob=coeffvalues(cf1_th); %coefficients p(x) global
    
    St_glob=cf_vls_St_glob(4) + (cf_vls_St_glob(1) - cf_vls_St_glob(4))./((1+(Kn./cf_vls_St_glob(3)).^cf_vls_St_glob(2)).^cf_vls_St_glob(5));
    
    Qtrans_LP5 = St_trans_LP5 .* gr_trans(1:end,3) .* Vinf .* cp .* (gr_trans(1:end,4) - Twi); %max heat rate in transitional regime
    Q_c_LP5 =Stcfr_memo.* rho(1:length(Stcfr_memo),1) .* Vinf .* cp .* (T0(1:length(Stcfr_memo),1) - Twi); %max heat rate in continuum regime
    Q_fm_LP5=Stfm1_memo.*(rho_fm./2)*Vinf^3; %max  heat rate in fm regime
    
    % % CHRAD=4.736*10^4.*4.9.^(1.072.*10^6.*(Vinf^(-1.88)).*rho_fm.^(-0.325)).*rho_fm.*1.22.*0.06; % radiative heat transfer coefficient  Tauber Sutton (stagnation point heating relations for earth and mars entries)
    % % Q_fmrad_LP5=CHRAD.*0.5.*rho_fm.*Vinf^3*10e4;%max radiation heat rate in fm regime
    % %
    % % Q_fmglobal=Q_fmrad_LP5+Q_fm_LP5;
    
    Q_global_LP5=[Q_c_LP5; Qtrans_LP5; Q_fm_LP5]; %max heat rate along all regimes
    
    % Plots for Qav
      %Qav vs Kn
%     figure(3)
%     xlabel('Kn')
%     ylabel('Qaveraged')
%     hold on
%     scatter(Knc,Qav_c,'b')
%     hold on
%     scatter(Kntrans,Qav_t ,'g');
%     hold on
%     scatter(Knfm,Qav_fm,'r')
%     hold on
%     plot(Kn,Qav,'k');
%     hold on
%     set(gca,'xscale','log')
%     set(gca,'yscale','log')
%     xlabel('Kn')
%     ylabel('Qav [W/M^2]')
    
    
end %Thermo dynamics module ends


if ThermoFlag == 0;
    Knc=Knc(Knc~=0)';
    Kntrans=Kntrans(Kntrans~=0)';
    Knfm=Knfm(Knfm~=0)';
end


%% ************************ AERODYNAMICS  ****************************%%%%
%      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if AeroFlag == 1;
    %Building the matrixes [St, Kn, Cx, Kn] for  continuum, trans, fm regime
    
    % Knc=Knc(Knc~=0)';
    % Kntrans=Kntrans(Kntrans~=0)';
    % Knfm=Knfm(Knfm~=0)';
    
    CDc_memo=[CDc_memo(CDc_memo~=0)'];
    CDtrans_memo=[CDtrans_memo(CDtrans_memo~=0)'];
    CDfm_memo=[CDfm_memo(CDfm_memo~=0)'];
    
    CLc_memo=[CLc_memo(CLc_memo~=0)'];
    CLtrans_memo=[CLtrans_memo(CLtrans_memo~=0)'];
    CLfm_memo=[CLfm_memo(CLfm_memo~=0)'];
    
    CMc_memo=[CMc_memo(CMc_memo~=0)'];
    CMtrans_memo=[CMtrans_memo(CMtrans_memo~=0)'];
    CMfm_memo=[CMfm_memo(CMfm_memo~=0)'];
    CSfm_memo=[CSfm_memo(CSfm_memo~=0)'];
    
    LDc_memo=CLc_memo./CDc_memo;
    LDtrans_memo=CLtrans_memo./CDtrans_memo;
    LDfm_memo=CLfm_memo./CDfm_memo;
    
    
    
    % creating the aerodynamics Bridging functions and generating the plots
    
    [cf1]=L5P([Knc(end)' (Knc(end)'+Kntrans(1)')/2  Kntrans'  Knfm(1)' ],[CDc_memo(end)'  (CDc_memo(end)'+CDtrans_memo(1)')/2  CDtrans_memo'  CDfm_memo(1)']);
    [cf3]=L5P([Knc' Kntrans' Knfm'],[CDc_memo' CDtrans_memo' CDfm_memo']);
    
    cf_vls=coeffvalues(cf1); %coefficienti p(x) bridging
    cf_vls_global=coeffvalues(cf3); %coefficienti p(x) global
    
    %CD,CL,LD,CM transitonal (LP5 type)
    CDtrans_LP5= CDfm_memo(1) + (CDc_memo(end) - CDfm_memo(1))./((1+(Kntrans./cf_vls(3)).^cf_vls(2)).^cf_vls(5));   
    CLtrans_LP5= CLfm_memo(1) + (CLc_memo(end) - CLfm_memo(1))./((1+(Kntrans./cf_vls(3)).^cf_vls(2)).^cf_vls(5));   
    CMtrans_LP5= CMfm_memo(1) + (CMc_memo(end) - CMfm_memo(1))./((1+(Kntrans./cf_vls(3)).^cf_vls(2)).^cf_vls(5));   
    CStrans_LP5= CSfm_memo(1) + (CSc_memo(end) - CSfm_memo(1))./((1+(Kntrans./cf_vls(3)).^cf_vls(2)).^cf_vls(5));   
    LDtrans_LP5=CLtrans_LP5./CDtrans_LP5;
    
    %CD,CL,LD,CM global (LP5 type)
    CD_LP5= CDfm_memo(1) + (CDc_memo(end) - CDfm_memo(1))./((1+(Kn./cf_vls_global(3)).^cf_vls_global(2)).^cf_vls_global(5));
    CL_LP5= CLfm_memo(1) + (CLc_memo(end) - CLfm_memo(1))./((1+(Kn./cf_vls_global(3)).^cf_vls_global(2)).^cf_vls_global(5));
    CM_LP5= CMfm_memo(1) + (CMc_memo(end) - CMfm_memo(1))./((1+(Kn./cf_vls_global(3)).^cf_vls_global(2)).^cf_vls_global(5));
    CS_LP5= CSfm_memo(1) + (CSc_memo(end) - CSfm_memo(1))./((1+(Kn./cf_vls_global(3)).^cf_vls_global(2)).^cf_vls_global(5));
    LD_LP5=CL_LP5./CD_LP5;
    
    CD_LP5_Kn=[CD_LP5 Kn];
    CL_LP5_Kn=[CL_LP5 Kn];
    CM_LP5_Kn=[CM_LP5 Kn];
    CS_LP5_Kn=[CS_LP5 Kn];
    LD_LP5_Kn=[LD_LP5 Kn];
    
    % Final CD,CL,LD,CM continuum,transitonal,fm
    for i=1:length(CD_LP5_Kn)
        
        if CD_LP5_Kn(i,2)>limKn_inf && CD_LP5_Kn(i,2)<limKn_sup
            CDtrans_DEF(i)=[CD_LP5_Kn(i,1)];
            CLtrans_DEF(i)=[CL_LP5_Kn(i,1)];
            CMtrans_DEF(i)=[CM_LP5_Kn(i,1)];
            CStrans_DEF(i)=[CS_LP5_Kn(i,1)];
            LDtrans_DEF(i)=[LD_LP5_Kn(i,1)];
        end
        
        if CD_LP5_Kn(i,2)<=limKn_inf
            CDcont_DEF(i)=[CD_LP5_Kn(i,1)];
            CLcont_DEF(i)=[CL_LP5_Kn(i,1)];
            CMcont_DEF(i)=[CM_LP5_Kn(i,1)];
            CScont_DEF(i)=[CS_LP5_Kn(i,1)];
            LDcont_DEF(i)=[LD_LP5_Kn(i,1)];
        end
        
        if CD_LP5_Kn(i,2)>=limKn_sup
            CDfm_DEF(i)=[CD_LP5_Kn(i,1)];
            CLfm_DEF(i)=[CL_LP5_Kn(i,1)];
            CMfm_DEF(i)=[CM_LP5_Kn(i,1)];
            CSfm_DEF(i)=[CS_LP5_Kn(i,1)];
            LDfm_DEF(i)=[LD_LP5_Kn(i,1)];
        end
        
    end
    
    
    CDcont_DEF=CDcont_DEF(CDcont_DEF~=0)';
    CLcont_DEF=CLcont_DEF(CLcont_DEF~=0)';
    CMcont_DEF=CMcont_DEF(CMcont_DEF~=0)';
    CScont_DEF=CScont_DEF(CScont_DEF~=0)';
    LDcont_DEF=LDcont_DEF(LDcont_DEF~=0)';
    
    CDtrans_DEF=CDtrans_DEF(CDtrans_DEF~=0)';
    CLtrans_DEF=CLtrans_DEF(CLtrans_DEF~=0)';
    CMtrans_DEF=CMtrans_DEF(CMtrans_DEF~=0)';
    CStrans_DEF=CStrans_DEF(CStrans_DEF~=0)';
    LDtrans_DEF=LDtrans_DEF(LDtrans_DEF~=0)';
    
    
    CDfm_DEF=CDfm_DEF(CDfm_DEF~=0)';
    CLfm_DEF=CLfm_DEF(CLfm_DEF~=0)';
    CMfm_DEF=CMfm_DEF(CMfm_DEF~=0)';
    CSfm_DEF=CSfm_DEF(CSfm_DEF~=0)';
    LDfm_DEF=LDfm_DEF(LDfm_DEF~=0)';
    
%     %  Plots for CD CL
%     %CD vs Kn
%     figure(1)
%     xlabel('Kn')
%     ylabel('CD')
%     hold on
%     scatter(Knc,CDcont_DEF,'b')
%     hold on
%     scatter(Kntrans,CDtrans_DEF,'g');
%     hold on
%     scatter(Knfm,CDfm_DEF,'r')
%     hold on
%     plot(Kn,CD_LP5,'k','LineWidth',0.5);
%     hold on
%     set(gca,'xscale','log')
%     xlabel('Kn')
%     ylabel('C_D')
%     
%     %CL vs Kn
%     figure(2)
%     xlabel('Kn')
%     ylabel('CL')
%     hold on
%     scatter(Knc,CLcont_DEF,'b')
%     hold on
%     scatter(Kntrans,CLtrans_DEF ,'g');
%     hold on
%     scatter(Knfm,CLfm_DEF,'r')
%     hold on
%     plot(Kn,CL_LP5 ,'k');
%     hold on
%     set(gca,'xscale','log')
%     xlabel('Kn')
%     ylabel('C_L')
    
    % FINAL AERODYNAMICS OUTPUT MATRIX CD-CL-L/D-CM-CS
    % OUTPUT_AER=[CD_LP5,CL_LP5,LD_LP5,CM_LP5,CS_LP5];
    
end

%% Simulation finished providing results depending on the used modules


disp(['simulation completed in: ',num2str(toc(tstart)),'s'])

if AeroFlag == 1 && ThermoFlag == 1;
    %     [Altitude CD CL Qav Kn]
%     F = [h(Hindex(1)) CD_LP5(Hindex(1)) CL_LP5(Hindex(1)) Qav(Hindex(1)) Kn(Hindex(1))];  %providing just the single input altitude results

    F = [h CD_LP5 CL_LP5 Qav Kn];     %to have all the altiutudes & Kn    
elseif AeroFlag == 1 && ThermoFlag == 0;
    disp('Qav has not been computed, the module was set to off')
%     F = [h(Hindex(1)) CD_LP5(Hindex(1)) CL_LP5(Hindex(1)) 0 Kn(Hindex(1))];

    F = [h CD_LP5 CL_LP5 zeros(length(h),1) Kn];     %to have all the altiutudes & Kn    
elseif AeroFlag == 0 && ThermoFlag == 1;
    disp('Aerodynamics have not been computed, the module was set to off')
%     F = [h(Hindex(1)) 0 0 Qav(Hindex(1)) Kn(Hindex(1))];

    F = [h zeros(length(h),1) zeros(length(h),1) Qav Kn];     %to have all the altiutudes & Kn
end


end
