clear all
close all
clc
warning('off')
tic
addpath('C:\Users\GBenedetti\Google Drive\TESI LM\tesi\tps\FOSTRAD\FOSTRAD 2.0\altmany-export_fig-2bad961');
addpath('C:\Users\GBenedetti\Google Drive\TESI LM\tesi\tps\FOSTRAD\FOSTRAD 2.0\external functions');
                                                        %% %%%%% FOSTRAD   V. 2.0 %%%%%%%%%
                                         
%%%%% Read in the STL file for the Geometry %%%%%%%%%
% [points,triangles,tri norms]

[V, F, N] = import_stl_fast('C:\Users\GBenedetti\Google Drive\TESI LM\tesi\tps\FOSTRAD\FOSTRAD 2.0\CAD MODELS\shuttle_noRS25_ascii.stl',1);

%[F,V,N] = stlread('C:\Users\GBenedetti\Google Drive\TESI LM\tesi\tps\FOSTRAD\aerothermal-master\aerothermal-master\CAD MODELS\GOCE_simplified_FOSTRAD_coarse.stl');

b =37;          %Characteristic Length of Object

%%%%% Truncating for unwanted rows %%%%%%%
if length(N) > length(F) 
    N = N(1:end-1,:);    
end

%%%%%% Calculating Elemental Areas %%%%%%%%%
A = Area(F,V);


tr = triangulation(F,V);
In = incenter(tr);
CG = COG(In);

%%%%%% Cleaning up duplicates %%%%%%%%
[ N, F, A ] = cleanup( N, F, A );

%%%%%% User Inputs %%%%%%%

SPAZAL=1; %INTERVALLO DI SPAZIATURA PER ALPHA
SPAZH=1; %INTERVALLO DI SPAZIATURA PER H km (lasciare fisso!)
alphamax=1;
alphamin=1;
hmin=60; %km
hmax=300; %km
beta = 20;       % Yaw angle, degrees

limKn_inf=5*10^-5;
limKn_sup=10;
Vinf = 7800;     % Free Stream Velocity, m/s

rN =1.329;     %Nose radius (m)

Lat= 55.8642;   % Latitude
Lon= -4.16;     % Longitude
year= 2010;     % year
yr_day= 1;      % day
UT_sec= 0;      % UT seconds


%%  Atmosphere model

[Ti,rhoi] = atmosnrlmsise00([hmin*1e3:SPAZH*1e3:hmax*1e3] , Lat, Lon, year, yr_day, UT_sec);
rhoi(:,[6,9]) = [ ];    %He O N2 O2 Ar H N [1/m^3] gas densities
Ti(:,1)= [];            %Free stream Temperatures

SN = 1.0;               %Normal Momentum Accommodation coefficient 
ST = 1.0;               %Tangential Momentum Accommodation Coefficient
AT = ST*(2-ST);         %Tangental Energy Accomodation Coefficient
AN = 1 - (1 - SN)^2;    %Normal Energy Accommodation Coefficient
AC = 0.5 * (AN + AT);   %Overall Energy Accommodation Coefficient

cp = 1004.7;            %specific heat at constant pressure (AIR)

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

xN2=percentual(:,3);
xO2=percentual(:,4);
xO= percentual(:,2);
xN= percentual(:,7);
xAr= percentual(:,5);
xHe= percentual(:,1);
xH= percentual(:,6);

m=mN2.*xN2+mO2.*xO2+mO.*xO+mN.*xN+mAr.*xAr+mHe.*xHe+mH.*xH; %Atmosphere molecular weight

amu = 1.66053892e-27;   % FMF properties %
mg = m.*amu;
Kb = 1.3806488e-23;     %Boltzmann constant, m2 kg s-2 K-1

Nav = 6.023e23;         %Avogadro Constant 
f=1.660539040e-27;      %conversion factor atomic mass unit - kg

h=[hmin:1:hmax];%Km
[Z, Z_L, Z_U, T, P, rho, c, g, mu, nu, k, n, n_sum] = atmo(hmax,1,1);

Z(1:hmin,:) = [] ;
T(1:hmin,:) = [] ;
P(1:hmin,:) = [] ;
g(1:hmin,:) = [] ;      %gravity [m/s^2]

% nu(1:hmin,:) = [] ;
% k(1:hmin,:) = [] ;


% Air density
if h(end)<=85;
rho = m.*P./(R.*Ti);
else
rho = rhoi*[mHe mO mN2 mO2 mAr mH mN]'./6.022169e26;
end

k = 2.64638e-3*Ti.^1.5./(Ti+245*10.^(-12./Ti));     %   Thermal Conductivity Coefficient

mu = 1.458e-6*Ti.^1.5./(Ti+110.4);                  % dynamic Viscosity

c=sqrt((gamma.*P)./rho);                             %sound speed 

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

Cpmax=(2./(gamma.*Minf.^2)).*((P02./P-1)); %Cp maximum

% FMF properties %
mg = m.*amu;
kb = 1.3806488e-23;      %Boltzmann constant, m2 kg s-2 K-1
vmp = sqrt(2.*kb.*Ti./mg);
s = Vinf./vmp;

l=(mu./P).*sqrt(pi.*kb.*Ti./(2.*mg));  %mean free path [m]

%Kn=l./rN;                                       %Knudsen Number
Kn = l./b;                               

Pr=mu.*cp./k; %Prandtl number

dudx = (1./rN) .* sqrt(2.*(P01-P)./rhos);         %Fay and Ridell Parameter

Qstfm = 0.5 * AC .* rho .* Vinf.^3;
Re0 = rho .* Vinf .* rN ./ mu_T0;
Stfm = 1;
Qsfr = 0.76.*(Pr.^-0.6).*((rhos.*mu_T0).^0.4).*((rhow.*mu_w).^0.1).*sqrt(dudx).*(h0-hw);
Qsfr1 = 0.94.*((rho.*mu).^0.4).*((rhow.*mu_w).^0.1).*sqrt(dudx).*(h0-hw); %fay riddel
Qsfr2=0.94.*((rho.*mu).^0.5).*sqrt(dudx).*(h0-hw); %Van Driest

%%  Parameter Initialization

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

%% Iteration Alpha

for alpha=alphamin:SPAZAL:alphamax
clear colornumber
clear ID_COL_OR
clear CDATA_OR_NUMBER
clear MATR_COLOR1
clear MATR_COLOR
clear CDATA_FIN
clear p
clear cf1
clear cf3

%free stream velocity vector
Vinfi = [Vinf * cosd(alpha) * cosd(beta), - Vinf * sind(beta), Vinf * sind(alpha) * cosd(beta)]; % %Free Stream Velocity Vector

I=1; 

%% PRE PROCESSING PHASE 1
% pre-selection leeward windward triangles

LEEWARD=0;
WINDWARD=0;


directionX = [1 0 0];
directionY = [0 1 0];
directionZ = [0 0 1];


a=zeros(1,length(N));
W=zeros(length(N),1);

for i=1:length(N)
a(i)=dot((Vinfi/norm(Vinfi)),N(i,:));

if a(i)<=0
    LEEWARD=LEEWARD+1;
    W(i)=1;
else
    WINDWARD=WINDWARD+1;
    W(i)=0;
end

end
W=W';
for j=1:length(N)
if W(j)==0;
N(j,:)= 0;
F(j,:)= 0;
A(j,:)= 0;
end
end

N( ~any(N,2), : ) = [];
F( ~any(F,2), : ) = [];
A( ~any(A,2), : ) = [];
W( ~any(W,2), : ) = [];

%% PRE PROCESSING PHASE 1 - Back-face culling Algorithm

for i=1:1:length(F)
ss = randomRGB();     %generazione random colore in EXA
G(i,:)=hex2rgb(ss);   %generazione matrice per assegnazioni colori random per ogni triangolo
end

myColorMap = G;

myColorMap1=im2double(myColorMap);  %conversione G in double [0 255] ----> [0 1]

colornumber(:,1:3)=[myColorMap(:,1).*256^2 myColorMap(:,2).*256 myColorMap(:,3)];

ID_COL_OR=colornumber(:,1)+colornumber(:,2)+colornumber(:,3); %ID original colors

%generate patch

figure('units','normalized','outerposition',[0 0 1 1])
p=patch('Faces',F,'Vertices',V,'FaceVertexCdata',myColorMap./255,'FaceColor','flat','LineStyle','none');

hold on
rotate(p,directionZ,beta)
rotate(p,directionY,alpha)
axis off
view([-90,0,0])
hold off
print('Front_Shuttle','-dpng','-r0') %Salva immagine %for a lower
% accurate analysis, low computational time
%imageData = export_fig('Front_Shuttle','-png'); %higher accurate analysis, high comp. time

CDATA_OR=get(p,'CData'); %Patch color data
CDATA_OR=im2uint8(CDATA_OR);

%patch color channels
CDATA_OR_R=CDATA_OR(:,:,1);
CDATA_OR_B=CDATA_OR(:,:,2);
CDATA_OR_G=CDATA_OR(:,:,3);

CDATA_OR=[CDATA_OR_R CDATA_OR_B CDATA_OR_G ]; %Patch color data [0 255]

CDATA_OR32=uint32(CDATA_OR);

CDATA_OR_NUMBER(:,1:3)=[CDATA_OR32(:,1).*256^2 CDATA_OR32(:,2).*256 CDATA_OR32(:,3)]; %univoco 

ID_COL_DATA_OR=(CDATA_OR_NUMBER(:,1)+CDATA_OR_NUMBER(:,2)+CDATA_OR_NUMBER(:,3));

%read image from folder

AA= imread('C:\Users\GBenedetti\Google Drive\TESI LM\tesi\tps\FOSTRAD\aerothermal-master\aerothermal-master\Front_Shuttle.png');   

figure('units','normalized','outerposition',[0 0 1 1])

O=image(AA,'CDataMapping','direct'); 
set(gca,'xtick',[],'ytick',[]);
B1=im2double(AA); 

CDATA_FIN=get(O,'CData');       %Image color data
CDATA_FIN32=uint32(CDATA_FIN);  %Image color data [0 255]

clear i
clear k

CDATA_FIN_NUMBER=zeros(size(AA,1),size(AA,2));

            for i=1:size(CDATA_FIN32,1)%rassegna righe px
                for k=1:size(CDATA_FIN32,2) %rassegna colonnne px
                    CDATA_FIN_R=CDATA_FIN32(i,k,1)*256^2;
                    CDATA_FIN_G=CDATA_FIN32(i,k,2)*256;
                    CDATA_FIN_B=CDATA_FIN32(i,k,3);    
                    CDATA_FIN_NUMBER(i,k)=CDATA_FIN_R+CDATA_FIN_B+CDATA_FIN_G;
                end
            end
    
% eliminare elementi uguali CDATA_FIN_NUMBER
clear a
a = unique(CDATA_FIN_NUMBER);
out = [a,histc(CDATA_FIN_NUMBER(:),a)]; % mi dà i numeri asssociati ai colori che compaiono nell'immagine e la loro frequenza (n pixel con quel colore)
CDATA_FIN_NUMBER_UNICI=out(:,1);
CDATA_FIN_NUMBER_UNICI(end, :) = [];    % rimozione colore bianco (sfondo)

clear i
clear y

% creazione matrici per confronto liste
MATR_COLOR=kron(double(ID_COL_DATA_OR),ones(1,length(CDATA_FIN_NUMBER_UNICI)));
CDATA_FIN_NUMBER_UNICI_tr=CDATA_FIN_NUMBER_UNICI';
MATR_COLOR1=zeros(length(F),length(CDATA_FIN_NUMBER_UNICI_tr));

            for i=1:length(F)   
               MATR_COLOR1(i,:)=CDATA_FIN_NUMBER_UNICI_tr;
            end


%comparazione lista colori iniziali con colori matrice

check=MATR_COLOR-double(MATR_COLOR1);
Ch=find(~check);  %cerca elementi di check uguali a zero. Quelli elementi saranno i colori 
                  %utilizzati nei triangolini frontali della mesh di partenza

[row,col]=find(~check);      %row contiene i numeri dei triangolini in vista!!!!       

VIEW_TRIANGLES=row; %triangolini in vista da usare per il calcolo
                    
% Selezione triangolini necessari all'analisi

num_tri=linspace(1,length(F),length(F))';
Fmod=[F num_tri];

INTER_TRIANGLES=intersect(Fmod(:,4),row);

NON_INT_TRIANGLES=setxor(Fmod(:,4),row);

%re-define F,N,A
F_new=F;
N_new=N;
A_new=A;

F_new([NON_INT_TRIANGLES(:)], :) = [];
N_new([NON_INT_TRIANGLES(:)], :) = [];
A_new([NON_INT_TRIANGLES(:)], :) = [];

B2WA = [cosd(alpha)*cosd(beta), -sind(beta), sind(alpha)*cosd(beta); ...
cosd(alpha)*sind(beta), cosd(beta), sind(alpha)*sind(beta);-sind(alpha), 0, cosd(alpha)];

%% iterazione quota 
quota=hmin;
for H=1:SPAZH:hmax-hmin+1 %Km
    
quota=quota+SPAZH;

%% Parameter Re-Initialization


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

% recalling variables

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
    Twi=Tw(H);
    Re0i=Re0(H);
    Qstfmi=Qstfm(H);
    
    if Kni >= 10            
        Qstfmi = Qstfm(H);
    else
        Re0i = Re0(H);
        Qsfri = Qsfr(H);
    end

% iterazione triangles 
 for i = 1:length(N)
 
    theta = acosd(dot(-N(i,:),Vinfi)/norm(-N(i,:))/norm(Vinfi));
    delta = 90 - theta;     %Local Inclination Angle
    r = In(i,:) - CG;    %Moment Arm
    
    %%%%%%%%%%%%% CONTINUUM AERODYNAMICS %%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
        if delta > 0
        Cpc(i,1) = Cpmaxi * (sind(delta))^2;
        Sref = Sref + A(i,1) * sind(delta);

        end
     
        Cpnc = Cpc(i,1) * -N(i,:) * A(i,1);
        Cc = Cc + Cpnc;
        Mc = Mc + cross(r,Cpnc);
    
    %%%%%%%%%%% FREE MOLECULAR AERODYNAMICS %%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
        theta = asind (dot(-Vinfni,N(i,:)));

        t = [0,0,0];
        Cpfm = 0;
        Ctfm = 0;
    
        if theta == 0
            t = (N(i,:) * dot(Vinfi,N(i,:)) - Vinfi) / sqrt(1 - dot(Vinfi,N(i,:))^2);        
            Ctfm = -(ST*cosd(theta)/si/sqrt(pi)) * (exp(-(si*sind(theta))^2) +...
                sqrt(pi) * si *sind(theta)*(1+erf(si*sind(theta))));
        elseif theta == 90
            Cpfm1 = exp(-(si*sind(theta))^2) * ((2-SN)*si*sind(theta)/sqrt(pi) +...
                SN*sqrt(Twi/Tii)/2);
            Cpfm2 = (1+erf(si*sind(theta))) * ((2-SN) * ((si*sind(theta))^2 + 0.5) +...
                0.5*SN*si*sind(theta)*sqrt(pi*Twi/Tii));
            Cpfm = (Cpfm1+Cpfm2)/si^2;
        elseif theta ~= 0 && theta ~= 90 && theta ~= -90
            t = (N(i,:) * dot(Vinfni,N(i,:)) - Vinfni) / sqrt(1 - dot(Vinfni,N(i,:))^2);
            Cpfm1 = exp(-(si*sind(theta))^2) * ((2-SN)*si*sind(theta)/sqrt(pi) +...
                SN*sqrt(Twi/Tii)/2);
            Cpfm2 = (1+erf(si*sind(theta))) * ((2-SN) * ((si*sind(theta))^2 + 0.5) +...
                0.5*SN*si*sind(theta)*sqrt(pi*Twi/Tii));
            Cpfm = (Cpfm1+Cpfm2)/si^2;
            Ctfm = -(ST*cosd(theta)/si/sqrt(pi)) * (exp(-(si*sind(theta))^2) + sqrt(pi) *...
                si *sind(theta)*(1+erf(si*sind(theta))));
        end
    
        Cpnfm = Cpfm * -N(i,:) * A(i,1);
        Cttfm = Ctfm * t * A(i,1);
        Cfm = Cfm + Cpnfm + Cttfm;
        Mfm = Mfm + cross(r,Cpnfm+Cttfm);
 
%     %%%%%%%%%%% CONTINUUM AERODYNAMICS HEATING %%%%%%%%%%%
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     
%         if Kni < limKn_sup
%             Stsc = 2.1/Re0i;        %Stagnation Point Stanton Number, SCARAB Formulation
%             Stsfr = Qsfri/rhoi/Vinf/(h0i-hwi);
%         end
%      
%     %%%%%%%%%%% TRANSLATIONAL HEATING %%%%%%%%%%%%%%%%%%%%
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     
%         if Kni > limKn_inf && Kni < limKn_sup
%             Stsc = Stsc/sqrt(1+(Stsc/Stfm)^2);      %Translational Stagnation Point Stanton Number, SCARAB Bridging Function
%             Stsfr = Stsfr/sqrt(1+(Stsfr/Stfm)^2);      %Translational Stagnation Point Stanton Number, SCARAB Bridging Function
%  
%         end
%     
%         if Kni < limKn_sup
%         if theta >= 0 && theta <= 90
%             Stc(i,1) = Stsc * (0.7*sind(delta)); 
%             Stcfr(i,1) = Stsfr * (0.1+0.9*sind(delta)); 
%             Q(i,1) = Stc(i,1) * rhoi * Vinf * cp * (T0i - Twi);
%             Qfr(i,1) = Stcfr(i,1) * rhoi * Vinf * cp * (T0i - Twi);
%             Total_heat(i,1) = Stc(i,1)*A(i,1);
%             
%             Chc=2*Q(i,1)/(rhoi*Vinf^3); %heat transfer coefficient
%         end
%         end
%     
%     %%%%%%%%%%%%%%%%% FREE MOLECULAR FLOW %%%%%%%%%%%%%%%
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     
        if Kni >= limKn_sup
            Q(i,1) = (Qstfmi/(si^3)/2/sqrt(pi))*((si^2+gammai/(gammai-1)-((gammai+1)*Twi/2/(gammai-1)/Tii))*...
            (exp(-(si*sind(delta))^2)+sqrt(pi)*si*sind(delta)*(1+erf(si*sind(delta))))-...
             0.5*exp(-(si*sind(delta))^2)); %detra,Kemp,riddel formula
        end
%      
end
   
Cc = Cc/Sref;
Cfm = Cfm/Sref;

Mc = Mc/Sref/b;
Mfm = Mfm/Sref/b;


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

% matrici [Cx, Kn] per continuum, trans, fm regime

        if Kni>limKn_inf && Kni<limKn_sup
        CDtrans_memo(H)=[CDtrans(CDtrans~=0)'];
        CLtrans_memo(H)=[CLtrans(CLtrans~=0)'];
        CMtrans_memo(H)=[Cmtrans(Cmtrans~=0)'];
        CStrans_memo(H)=[CStrans(CStrans~=0)'];
        Kntrans(H)=Kni(Kni~=0)';
        end

        if Kni>=limKn_sup
        CDfm_memo(H)=[CDfm(CDfm~=0)'];
        CLfm_memo(H)=[CLfm(CLfm~=0)'];
        CMfm_memo(H)=[Mmfm(Mmfm~=0)'];
        CSfm_memo(H)=[CSfm(CSfm~=0)'];
        Knfm(H)=Kni(Kni~=0)';
        end

        if Kni<=limKn_inf
        CDc_memo(H)=[CDc(CDc~=0)'];
        CLc_memo(H)=[CLc(CLc~=0)'];
        CMc_memo(H)=[Mmc(Mmc~=0)'];
        CSc_memo(H)=[CSc(CSc~=0)'];
        Knc(H)=Kni(Kni~=0)';
        end

I=I+1;
C='ciclo quota km';
disp(C);
disp(I);

end
close all


%% matrici [Cx, Kn] per continuum, trans, fm regime

Knc=Knc(Knc~=0)';
Kntrans=Kntrans(Kntrans~=0)';
Knfm=Knfm(Knfm~=0)';

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


%% generazione bridging LP5 e plot completi 

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

%% Final CD,CL,LD,CM continuum,transitonal,fm
for i=1:length(CD_LP5_Kn)

    if CD_LP5_Kn(i,2)>limKn_inf & CD_LP5_Kn(i,2)<limKn_sup
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

%% grafici finali Cx
%CD vs Kn
figure(1)
xlabel('Kn')
ylabel('CD')
hold on
scatter(Knc,CDcont_DEF,'b')
hold on
scatter(Kntrans,CDtrans_DEF,'g');
hold on
scatter(Knfm,CDfm_DEF,'r')
hold on
plot(Kn,CD_LP5,'k','LineWidth',0.5);
hold on
set(gca,'xscale','log')
xlabel('Kn')
ylabel('C_D')

%CL vs Kn
figure(2)
xlabel('Kn')
ylabel('CL')
hold on
scatter(Knc,CLcont_DEF,'b')
hold on
scatter(Kntrans,CLtrans_DEF ,'g');
hold on
scatter(Knfm,CLfm_DEF,'r')
hold on
plot(Kn,CL_LP5 ,'k');
hold on
set(gca,'xscale','log')
xlabel('Kn')
ylabel('C_L')

%L/D vs Kn
figure(3)
xlabel('Kn')
ylabel('LD')
hold on
scatter(Knc,LDcont_DEF,'b')
hold on
scatter(Kntrans,LDtrans_DEF ,'g');
hold on
scatter(Knfm,LDfm_DEF,'r')
hold on
plot(Kn,LD_LP5,'k');
hold on
set(gca,'xscale','log')
xlabel('Kn')
ylabel('L/D')


%CM vs Kn
figure(4)
xlabel('Kn')
ylabel('CM')
hold on
scatter(Knc,CMcont_DEF,'b')
hold on
scatter(Kntrans,CMtrans_DEF ,'g');
hold on
scatter(Knfm,CMfm_DEF,'r')
hold on
plot(Kn,CM_LP5,'k');
hold on
xlabel('Kn')
ylabel('C_M')
set(gca,'xscale','log')

%CS vs Kn
figure(13)
xlabel('Kn')
ylabel('CS')
hold on
scatter(Knc,CScont_DEF,'b')
hold on
scatter(Kntrans,CStrans_DEF ,'g');
hold on
scatter(Knfm,CSfm_DEF,'r')
hold on
plot(Kn,CS_LP5,'k');
hold on
xlabel('Kn')
ylabel('C_S')
set(gca,'xscale','log')

X='ciclo alpha n.';
disp(X);
disp(alpha);
toc 

%% FINAL AERODYNAMICS OUTPUT MATRIX CD-CL-L/D-CM-CS
OUTPUT_AER=[CD_LP5,CL_LP5,LD_LP5,CM_LP5,CS_LP5];
end
toc