clear all
clc
close all

warning('off')
tic
addpath('C:\Users\GBenedetti\Google Drive\TESI LM\tesi\tps\FOSTRAD\FOSTRAD 2.0\altmany-export_fig-2bad961');
addpath('C:\Users\GBenedetti\Google Drive\TESI LM\tesi\tps\FOSTRAD\FOSTRAD 2.0\external functions');
%%%%% Read in the STL file for the Geometry %%%%%%%%%
% [points,triangles,tri norms]

%%%% OBJECT: GOCE  %%%%
[V, F, N] = import_stl_fast('C:\Users\GBenedetti\Google Drive\TESI LM\tesi\tps\FOSTRAD\FOSTRAD 2.0\CAD MODELS\sphere1_ascii.stl',1);
%[F,V,N] = stlread('C:\Users\GBenedetti\Google Drive\TESI LM\tesi\tps\FOSTRAD\FOSTRAD 2.0\CAD MODELS\GOCE_simplified_FOSTRAD_veryfine.stl');


%%%%% Truncating for unwanted rows %%%%%%%
if length(N) > length(F) 
    N = N(1:end-1,:);    
end

%%%%%% Calculating Elemental Areas %%%%%%%%%
A = Area(F,V);

tr = triangulation(F,V);

In = incenter(tr);
CG = COG(In);

% traslazione sist riferimento default- sistema rif centrato in CG

V(:,1)=V(:,1)-CG(1);
In(:,1)=In(:,1)-CG(1);    
  
V(:,2)=V(:,2)-CG(2);      
In(:,2)=In(:,2)-CG(2);  

V(:,3)=V(:,3)+ abs(CG(3));  
V(:,3)=V(:,3)+abs(CG(3));  

%traslazione sistema rif centrato in CG - sist rif centrato in nose
V(:,1)=V(:,1)+abs(min(In(:,1)))+0.2837;
In(:,1)=In(:,1)+abs(min(In(:,1)))+0.2837;   

%%%%%% Cleaning up duplicates %%%%%%%%
[ N, F, A ] = cleanup( N, F, A );

b=max(In(:,1)); %characteristic length of object


%%% determination x/b, y/b, z/b coordinates

In_xL=[In(:,1)./b (1:1:length(In))']; %coordinate X/L incentro 
In_xL1=[roundx(In_xL(:,1), 1) (1:1:length(In))'];

In_yL=[In(:,2)./b (1:1:length(In))']; %coordinate y/L incentro
In_yL1=[roundx(In_yL(:,1), 1) (1:1:length(In))'];

In_zL=[In(:,3)./b (1:1:length(In))']; %coordinate z/L incentro
In_zL1=[roundx(In_zL(:,1), 1) (1:1:length(In))'];

InLD=[In_xL1(:,1) In_yL1(:,1) In_zL1]; %coordinate /b incentro + numero triangolo corrispondente


%windward centerline coordinates== 0<x/L<1 || y/L=0  || z<0
In_ctrl=InLD(InLD(:,3)<2.336/b,:); %tutti i triangoli con z<0

In_ctrl=In_ctrl((In_ctrl(:,2)==0),:);

triangolixLD=In_xL1(In_xL1(:,1)==0.0 | In_xL1(:,1)==0.1 | In_xL1(:,1)==0.2 | In_xL1(:,1)==0.3 | In_xL1(:,1)==0.4 ...
    | In_xL1(:,1)==0.5 | In_xL1(:,1)==0.6 | In_xL1(:,1)==0.7 | In_xL1(:,1)==0.8 | In_xL1(:,1)==0.9 ...
    | In_xL1(:,1)==1,:);

triangoliyLD=In_yL1(In_yL1(:,1)==0.0 | In_yL1(:,1)==0.05 | In_yL1(:,1)==0.1 | In_yL1(:,1)==0.15 | In_yL1(:,1)==0.2 ...
    | In_yL1(:,1)==0.25 |In_yL1(:,1)==0.3,:);

triangolixLD_ctrl=In_ctrl(In_ctrl(:,1)==0.0 | In_ctrl(:,1)==0.1 | In_ctrl(:,1)==0.2 | In_ctrl(:,1)==0.3 | In_ctrl(:,1)==0.4 ...
    | In_ctrl(:,1)==0.5 | In_ctrl(:,1)==0.6 | In_ctrl(:,1)==0.7 | In_ctrl(:,1)==0.8 | In_ctrl(:,1)==0.9 ...
    | In_ctrl(:,1)==1,:);


%% For continumm regime %%%

%%%%%% User Inputs %%%%%%%


SPAZAL=1; %INTERVALLO DI SPAZIATURA PER ALPHA
SPAZH=1; %INTERVALLO DI SPAZIATURA PER H km (lasciare fisso!)
alphamax=-0.32;
alphamin=-0.32;
beta = 2.01;       % Yaw angle, degrees

hmin=50; %km
hmax=180; %km
Vinf =7810;     % Free Stream Velocity, m/s

Knboundfm=10; %Boundary Kn trans-fm regime
Knboundcont=10^-4;

Lat= 55.8642;   % Latitude
Lon= -4.16;     % Longitude
year= 2013;     % year
yr_day= 294;      % day
UT_sec= 0;      % UT seconds

SN = 1.0;   %Normal Momentum Accommodation coefficient 
ST = 1.0;   %Tangential Momentum Accommodation Coefficient
AT = ST*(2-ST);     %Tangental Energy Accomodation Coefficient
AN = 1 - (1 - SN)^2;    %Normal Energy Accommodation Coefficient
AC = 0.5 * (AN + AT);   %Overall Energy Accommodation Coefficient


%% show object
for i=1:1:length(F)
ss = randomRGB();     %generazione random colore in EXA
G(i,:)=hex2rgb(ss);   %generazione matrice per assegnazioni colori random per ogni triangolo
end
% 
% G=unique(randi([0 255],length(F),3),'rows'); %generazione random colore
% 
% G=uint8(G);         %conversione G in uint8

myColorMap = G;


myColorMap1=im2double(myColorMap);  %conversione G in double [0 255] ----> [0 1]

colornumber(:,1:3)=[myColorMap(:,1).*256^2 myColorMap(:,2).*256 myColorMap(:,3)];

ID_COL_OR=colornumber(:,1)+colornumber(:,2)+colornumber(:,3); %ID original colors

%generate patch

figure('units','normalized','outerposition',[0 0 1 1])
p=patch('Faces',F,'Vertices',V,'FaceVertexCdata',myColorMap./255,'FaceColor','flat','LineStyle','none');
%p=patch('Faces',F,'Vertices',V,'FaceVertexCdata',myColorMap1,'FaceColor','flat','LineStyle','none');

directionX = [1 0 0];
directionY = [0 1 0];
directionZ = [0 0 1];
hold on
rotate(p,directionZ,beta)
rotate(p,directionY,alphamax)
axis off
view([-90,0,0])
hold off
%% atmosphere calculation

[Ti,rhoi] = atmosnrlmsise00([hmin*1e3:SPAZH*1e3:hmax*1e3] , Lat, Lon, year, yr_day, UT_sec);
rhoi(:,[6,9]) = [ ]; %He O N2 O2 Ar H N [1/m^3] gas densities
Ti(:,1)= [];

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

%%% Parameter Initialization

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
Q_fm=zeros(length(N),1);
Stc = zeros(length(N),1);
Stcfr = zeros(length(N),1);
Stfm=zeros(length(N),1);
rN = 0.1; %nose radius


h=[hmin:1:hmax];%Km

[Z, Z_L, Z_U, T, P, rho, c, g, mu, nu, k, n, n_sum] = atmo(hmax,1,1); % other atmosphere parameters

Z(1:hmin,:) = [] ;
T(1:hmin,:) = [] ;
P(1:hmin,:) = [] ;
rho(1:hmin,:) = [] ;
% c(1:hmin,:) = [] ;
g(1:hmin,:) = [] ;

% nu(1:hmin,:) = [] ;
% k(1:hmin,:) = [] ;

%   Thermal Conductivity Coefficient
k = 2.64638e-3*Ti.^1.5./(Ti+245*10.^(-12./Ti));

%dynamic Viscosity
        if h(end)<=85;
        else
            mu= [ mu' , (mu(end)*ones ( hmax-86 ,1))' ]';
        end
mu(1:hmin,:) = [] ;

%sound speed 
c=sqrt((gamma.*P)./rho);

%Gas constant, kg m2 s-2 K-1 mol-1
cp = 1004.7;            %specific heat at constant pressure (AIR)
R = cp .* (gamma-1) ./ gamma;   

%mach number
Minf = norm(Vinf) ./ c;  

% Free Stream Stagnation Pressure and Temperature
P01 = P .* (1 + 0.5 .* (gamma-1) .* Minf.^2) .^ (gamma./(gamma-1)); 
T0 = Ti .* (1 + 0.5 .* (gamma-1) .* Minf.^2);

% Tw = Ti.*(1+((gamma-1)./2).*(Minf.^2)); %Surface Temperature 
% Tw2 = Ti.*(1+0.2.*(norm(Vinf)./c).^2);

%Prandtl number
Pr=mu.*cp./k;

%Reynolds
Re=rho.*Vinf.*b./mu; 

% B=(Pr./2).*((5.*Pr-7)./(5.*Pr+1));
% rec=1-(4.71-4.11.*B-0.601.*Pr).*Re.^(-0.2);     %Recovery factor Seban Formulation
% rec1=1-4.55.*(1-Pr).*Re.^(-0.2);                %Shirokov
% rec3=0.670.*Pr+0.322; %NACA TN 2296 0.65<Pr<1

% Recovery factor 
rect=Pr.^(1/3); %turbulent
recl=sqrt(Pr); %laminar

for y=1:1:(hmax-hmin)+1
    if Re(y)<=5*10^5
        REC(y)=recl(y);         %laminar region
    else
        REC(y)=rect(y);         %turbolent region
    end
end

%Surface Temperature 
Tw=REC'.*(T0-Ti)+Ti;

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

%Calculating mean Free Path and Knudsen

l=(mu./P).*sqrt(pi.*kb.*Ti./(2.*mg));  
% l=(mu(h+1)/pinf)*sqrt(pi*Kb*Tinf/(2*(m(h+1)*f))) 
%mean free path [m]
Kn = l./b;                                        %Knudsen Number
%Kn = l./rN;                               

rN = 1; %nose radius
dudx = (1./rN) .* sqrt(2.*(P01-P)./rhos);         %Fay and Ridell Parameter

Qstfm = 0.5 * AC .* rho .* Vinf.^3;
Re0 = rho .* Vinf .* rN ./ mu_T0;


% Qsfr = 0.76.*(Pr.^-0.6).*((rhos.*mu_T0).^0.4).*((rhow.*mu_w).^0.1).*sqrt(dudx).*(h0-hw);
% Qsfr1 = 0.94.*((rho.*mu).^0.4).*((rhow.*mu_w).^0.1).*sqrt(dudx).*(h0-hw); %fay riddel
% Qsfr2=0.94.*((rho.*mu).^0.5).*sqrt(dudx).*(h0-hw); %Van Driest
Qsfr = 0.76.*(Pr.^-0.6).*((rho.*mu).^0.4).*((rhow.*mu_w).^0.1).*sqrt(dudx).*(h0-hw); %ref  NASA TM 84580 eq6 Fay & Riddel

a=0;
e=0;
v=0;
%%
for alpha=alphamin:SPAZAL:alphamax
 
%free stream velocity vector
Vinfi = [Vinf * cosd(alpha) * cosd(beta), - Vinf * sind(beta), Vinf * sind(alpha) * cosd(beta)]; % %Free Stream Velocity Vector

I=1; 

B2WA = [cosd(alpha)*cosd(beta), -sind(beta), sind(alpha)*cosd(beta); ...
cosd(alpha)*sind(beta), cosd(beta), sind(alpha)*sind(beta);-sind(alpha), 0, cosd(alpha)];


quota=hmin;

    for H=1:SPAZH:(hmax-hmin)+1 %Km

    quota=quota+SPAZH;
    %%% Parameter Initialization

    % Cpc = zeros(length(N),1);
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


        %recalling variables
        Pi=P(H);
    %     Ri=R(H);
        gammai=gamma(H);
    %     ci=c(H);
        Minfi=Minf(H);

    %     P01i=P01(H);
        T0i=T0(H);
        mui = mu(H);
    %     mu_T0i = mu_T0(H);
        rhoi=rho(H);

        rhosi=rhos(H);
        rhowi=rhow(H);
        mu_wi=mu_w(H);
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
        Pri=Pr(H);
        Vinfni = Vinfi/Vinf;
    %     dudxi=dudx(H);
        Tii=Ti(H);
        Twi=(70+273.15); %Tw(H);
        Re0i=Re0(H);
        Qstfmi=Qstfm(H);
        Rei=Re(H);

         if Kni >= Knboundfm%free-m
            Qstfmi = Qstfm(H);
        else
            Re0i = Re0(H);
            Qsfri = Qsfr(H);
        end


     for i = 1:length(N)

        theta = acosd(dot(-N(i,:),Vinfi)/norm(-N(i,:))/norm(Vinfi));
        delta = 90 - theta;     %Local Inclination Angle
        r = In(i,:) - CG;       %Moment Arm


    %     %%%%%%%%%%% CONTINUUM AERODYNAMICS HEATING %%%%%%%%%%%
    %     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %     
             if Kni < Knboundcont
                Stsc = 2.1/sqrt(Re0i);        %Stagnation Point Stanton Number, SCARAB Formulation
                Stsfr = Qsfri/rhoi/Vinf/(h0i-hwi); %St Numb. Fay Riddel
             end

     %legare Stanton number all'i-esimo triangolino....

             if Kni < Knboundfm

                      if theta >= 0 && theta <= 90

                        Stc(i,1) = Stsc * (0.7*sind(delta)); %da dove spunta sta formula?????
                        Stcfr(i,1) = Stsc * (0.1+0.9*cosd(theta)); %modified Lees (SCARAB)
                        Stcfr_KDR(i,1) = Stsc * (cosd(theta/2)^5.27); % Kemp Rose Detra 
                        Stcfr_FOSTRAD20(i,1) = Stsc * (-0.3+0.7*cosd(theta)); % FOSTRAD2.0

                                if Stcfr_FOSTRAD20(i,1)<0
                                   Stcfr_FOSTRAD20(i,1)=0;
                                end
            
                      else
                        Stc(i,1)=0;
                        Stcfr(i,1)=0;
                        Stcfr_KDR(i,1)=0;
                        Stcfr_FOSTRAD20(i,1)=0;
                       end   


            Q(i,1) = Stc(i,1) * rhoi * Vinf* cp * (T0i - Twi);
            Qfr(i,1) = Stcfr(i,1) * rhoi * Vinf * cp * (T0i - Twi); % W/m^2
            Q_KDR(i,1)=Stcfr_KDR(i,1)* rhoi * Vinf * cp * (T0i - Twi);
            Q_FOSTRAD20(i,1)=Stcfr_FOSTRAD20(i,1)* rhoi * Vinf * cp * (T0i - Twi);

            Q_c(i,1)=Qfr(i,1); % Heat Transfer [W/m^2] continuum
            
            
    %         Total_heat(i,1) = Stc(i,1)*A(i,1);
    %         Total_heatfr(i,1) = Stcfr(i,1)*A(i,1); %il flusso termico non può
    %         avere le dimensioni di  un'area, essendo St un coefficiente
    %         adimensionale.....

            Total_heat(i,1) = Q_c(i,1)*A(i,1); %flusso termico (J/sec)
    

       end

    %     %%%%%%%%%%%%%%%%% FREE MOLECULAR FLOW %%%%%%%%%%%%%%%
    %     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %     
        if Kni >= Knboundfm 
        Stsfm= Qstfmi/rhoi/Vinf/(h0i-hwi);

    %       if theta >= 0 && theta <= 90
            Q_fm(i,1) = (Qstfmi/(si^3)/2/sqrt(pi))*((si^2+gammai/(gammai-1)-((gammai+1)*Twi/2/(gammai-1)/Tii))*...
            (exp(-(si*sind(delta))^2)+sqrt(pi)*si*sind(delta)*(1+erf(si*sind(delta))))-0.5*exp(-(si*sind(delta))^2));     
            %kempra fay riddel, rif introvabile 1957

            Stfm(i,1)=Q_fm(i,1)/(rhoi * Vinf * cp * (T0i - Twi)); %heat transfer coefficient fm
            Stfm1(i,1)=Q_fm(i,1)/((rhoi/2)*Vinf^3);

            Total_heat_fm(i,1) = Q_fm(i,1)*A(i,1);   
    %         end

%            figure(7)
%                     graf_3D=trisurf(F,V(:,1),V(:,2),V(:,3), Q_fm,'EdgeColor','none','LineStyle','none','FaceLighting','phong');
%                     colorbar
%                     title('Q_fm ')
%                     hold on
%                     rotate(graf_3D,directionX,beta)
%                     rotate(graf_3D,directionY,alphamax)
%                     xlabel('X')
%                     ylabel('Y')
%                     zlabel('Z')
%                     axis equal
%                     hold off 
        end     


     end

%                 figure(1)
%                 xlabel('kn')
%                 ylabel('St')
%                  if  Kni < Knboundcont
% %                 plot(Kni,max(Stcfr),'*');
%                 plot(Kni,max(Stcfr_SCAR_MOD),'*');
%                 hold on
%                  end
% 
%                 if Kni>Knboundfm
%                    plot(Kni,max(Stfm1),'*');
%                    hold on
%                 end
%                 set(gca,'xscale','log')
% % 
%                     figure(6)
%                     graf_3D=trisurf(F,V(:,1),V(:,2),V(:,3), Q_KDR_MOD,'EdgeColor','none','LineStyle','none','FaceLighting','phong');
%                     colorbar
%                     title('Q_c cont')    
%                     hold on
%                     rotate(graf_3D,directionX,beta)
%                     rotate(graf_3D,directionY,alphamax)
%                     xlabel('X')
%                     ylabel('Y')
%                     zlabel('Z')
%                      view([-90,0,0])
%                     axis equal
%                     hold off

    Q_c1 =[Q_c  (1:length(N))'];
    Q_c_LD1=Q_c1(triangolixLD_ctrl(:,4),:);
    
    Q_c2= [Q_KDR (1:length(N))'];
    Q_c_LD2=Q_c2(triangolixLD_ctrl(:,4),:);
    
    Q_c3= [ Q_FOSTRAD20 (1:length(N))'];
    Q_c_LD3=Q_c3(triangolixLD_ctrl(:,4),:);
   

%        
%                 figure(9)
%                 xlabel('x/L')
%                 ylabel('Heat transfer kW/m^2')
%                 title('Heat transfer along the Windward centerline , M=18.52, Vinf=5.39 km/s, \alpha= 74°, altitude= 74 km')

                jj=[0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1];
                
%                 for j=(1:length(jj))
%                 jjj=triangolixLD_ctrl(triangolixLD_ctrl(:,1)==jj(j),:);
%                 Q_cLD1=Q_c1(jjj(:,4),:);
%                 Q_cLD1=max(Q_cLD1(:,1));
%                     if any(triangolixLD_ctrl(:,1)==jj(j))
%                 hold on
%                 scatter(jj(j),Q_cLD1*10^-6,'b')
%                 hold off
%                     end
%                
%                 end
%                 hold on
%                 
%     
%                 for j=(1:length(jj))
%                 jjj=triangolixLD_ctrl(triangolixLD_ctrl(:,1)==jj(j),:);
%                 Q_cLD2=Q_c2(jjj(:,4),:);
%                 Q_cLD2=max(Q_cLD2(:,1));
%                     if any(triangolixLD_ctrl(:,1)==jj(j))
%                 hold on
%                 scatter(jj(j),Q_cLD2*10^-6,'g')
%                 hold off
%                     end
%                 
%                 end
%                 
%                 jj=[0:0.1:1];
%                 hold on
%                 for j=(1:length(jj))
%                 jjj=triangolixLD_ctrl(triangolixLD_ctrl(:,1)==jj(j),:);
%                 Q_cLD3=Q_c3(jjj(:,4),:);
%                 Q_cLD3=max(Q_cLD3(:,1));
%                     if any(triangolixLD_ctrl(:,1)==jj(j))
%                 hold on
%                 scatter(jj(j),Q_cLD3*10^-3*1.237,'b','LineWidth',2); %Mod SCARAB (4.81 fattore moltiplicativo per usare St glob)
% %                 hold on
% %                 scatter(jj(j),Q_cLD3*10^-3,'g');
%                 hold off
%                     end
%                 
%                 end
                
    
    I=I+1;
    C='ciclo quota n.';
    disp(C);
    disp(I);

         % matrici [St, Kn] per continuum,fm regime
         if Kni<=Knboundcont
%             Stcfr_memo(H)= max(Stcfr);
            Stcfr_memo(H)= max(Stcfr_FOSTRAD20);
            Kn_c(H)=Kni(Kni~=0)';
            rho_c(H)=rhoi(rhoi~=0)';
         end

         if Kni>=Knboundfm
            Stfm1_memo(H)= max(Stfm1);
            Kn_fm(H)=Kni(Kni~=0)';
            rho_fm(H)=rhoi(rhoi~=0)';
         end
     
    end
  
% %     %%%%%%%%%%% TRANSLATIONAL HEATING %%%%%%%%%%%%%%%%%%%%
% %     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if Kni>=Knboundcont 
     
       gr_trans=[h' Kn rho T0]; %grandezze fisiche campo transition
       gr_trans= gr_trans(gr_trans(:, 2) >= Knboundcont & gr_trans(:,2) <= Knboundfm,:);
       H_Knboundcont=gr_trans(1,1);
       H_Knboundfm=gr_trans(end,1);
       
              for H1=1:SPAZH:H_Knboundfm-H_Knboundcont %Km
                 
                 %recalling variables
                 Kni=gr_trans(H1,2);
                 rhoi=gr_trans(H1,3);
                 T0i=gr_trans(H1,4);
                 
                 Sttrans=(max(Stcfr)+ Kni.*max(Stfm1))./(1+ Kni); 
                 rho_tr =rho(H_Knboundcont-hmin:hmax-H_Knboundfm);
                 T0i_tr =T0(H_Knboundcont-hmin:hmax-H_Knboundfm);

                                 for i=1:length(N)
                                        theta = acosd(dot(-N(i,:),Vinfi)/norm(-N(i,:))/norm(Vinfi));
                                        delta = 90 - theta;     %Local Inclination Angle
                                        r = In(i,:) - CG;       %Moment Arm
                                            %if theta >= 0 && theta <= 90
                                                %Stctrans(i,1)= (Stcfr(i,1)+Kni*Stfm1(i,1))/(1+Kni); 
                                                Stctrans(i,1)=(Stcfr_FOSTRAD20(i,1)+Kni*Stfm1(i,1))/(1+Kni);
                                                
                                                Qtrans(i,1) = Stctrans(i,1) * rhoi * Vinf * cp * (T0i - Twi);
                                %               Qfr(i,1) = Sttrans(i,1) * rhoi * Vinf * cp * (T0i - Twi);
                                                Total_heat_trans(i,1) = Qtrans(i,1)*A(i,1);   
                                            %end     

                                 end
                                 
%                          figure(1)
%                          semilogx(Kni,max(Stctrans(:,1)),'*')
%                          hold on
                         

%                          figure(8)
%                             graf_3D=trisurf(F,V(:,1),V(:,2),V(:,3), Qtrans,'EdgeColor','none','LineStyle','none','FaceLighting','phong');
%                             colorbar
%                             title('Qtrans trans')
%                             hold on
%                             rotate(graf_3D,directionX,beta)
%                             rotate(graf_3D,directionY,alphamax)
%                             xlabel('X')
%                             ylabel('Y')
%                             zlabel('Z')
%                             axis equal
%                             hold off

         % matrici [St, Kn] per trans regime             
         if Kni>=Knboundcont && Kni<= Knboundfm
            Stctranns_memo(H1)= max(Stctrans);
            Kn_trans(H1)=Kni(Kni~=0)';
         end   
         
    end
                
    end

            
X='ciclo alpha n.';
disp(X);
disp(alpha);
toc 
end

%% matrici [St, Kn] per continuum, trans, fm regime

Kn_c=Kn_c(Kn_c~=0)';
% Kn_trans=Kn_trans(Kn_trans~=0)';
Kn_trans=gr_trans(:,2);

Kn_fm=Kn_fm(Kn_fm~=0)';
rho_fm=rho_fm(rho_fm~=0)';

Stcfr_memo=[Stcfr_memo(Stcfr_memo~=0)'];
Stctranns_memo=[Stctranns_memo(Stctranns_memo~=0)'];
Stfm1_memo=[Stfm1_memo(Stfm1_memo~=0)'];
rho_c=rho_c(rho_c~=0)';
%% generazione bridging LP5 per St transitional regime
glob_ref=[0.0014 0.0080 0.01 0.0224 0.028 0.0754  0.098 0.1 0.2044 0.227 0.476  0.8826 1 1.219  2.909  3.9908  ; 0.0224 0.0696 0.2615 0.170 0.176 0.4 0.405 0.573 0.665 0.595 0.687 0.802 0.9346 0.838 0.874 0.914 ]; %STS DATA (BELOW) AND SPHERE DATA

[cf1_th]=L5P([Kn_c; glob_ref(1,:)'; Kn_fm ],[Stcfr_memo; glob_ref(2,:)';Stfm1_memo]);

cf_vls_St=coeffvalues(cf1_th); %coefficienti p(x) bridging 

St_trans_LP5=Stfm1_memo(end) + (Stcfr_memo(1) - Stfm1_memo(end))./((1+(Kn_trans./cf_vls_St(3)).^cf_vls_St(2)).^cf_vls_St(5));

[cf2_th]=L5P([Kn_c; Kn_trans; Kn_fm],[Stcfr_memo;St_trans_LP5; Stfm1_memo]);
cf_vls_St_glob=coeffvalues(cf1_th); %coefficienti p(x) global

St_glob=cf_vls_St_glob(4) + (cf_vls_St_glob(1) - cf_vls_St_glob(4))./((1+(Kn./cf_vls_St_glob(3)).^cf_vls_St_glob(2)).^cf_vls_St_glob(5));

Qtrans_LP5 = St_trans_LP5 .* gr_trans(1:end,3) .* Vinf .* cp .* (gr_trans(1:end,4) - Twi); %max heat rate in transitional regime
Q_c_LP5 =Stcfr_memo.* rho(1:length(Stcfr_memo),1) .* Vinf .* cp .* (T0(1:length(Stcfr_memo),1) - Twi); %max heat rate in continuum regime
Q_fm_LP5=Stfm1_memo.*(rho_fm./2)*Vinf^3; %max  heat rate in fm regime

% % CHRAD=4.736*10^4.*4.9.^(1.072.*10^6.*(Vinf^(-1.88)).*rho_fm.^(-0.325)).*rho_fm.*1.22.*0.06; % radiative heat transfer coefficient  Tauber Sutton (stagnation point heating relations for earth and mars entries)
% % Q_fmrad_LP5=CHRAD.*0.5.*rho_fm.*Vinf^3*10e4;%max radiation heat rate in fm regime
% % 
% % Q_fmglobal=Q_fmrad_LP5+Q_fm_LP5;

Q_global_LP5=[Q_c_LP5; Qtrans_LP5; Q_fm_LP5]; %max heat rate along all regimes


%%
figure(15)
xlabel('kn')
ylabel('St')
% semilogx(Kn_c,Stcfr_memo)
hold on
% semilogx(Kn_trans,St_trans_LP5)
% hold on
% semilogx(Kn_trans,St_trans_LP5,'r','LineWidth',1.5)
% hold on
% semilogx(Kn_fm,Stfm1_memo)
hold on
% scatter(glob_ref(1,:),glob_ref(2,:))
hold on
semilogx(Kn,St_glob,'r','LineWidth',1.5)
hold on
plot([Knboundcont Knboundcont],get(gca,'ylim'),'-.b')
hold on
plot([Knboundfm Knboundfm],get(gca,'ylim'),'-.b')
hold on

% figure(16)
% plot(Qtrans_LP5*10^-3,gr_trans(1:end-1,1)) %max heat rate (kW/mq) in trans regime
% hold on
% plot(Q_c_LP5*10^-3,h(1:length(Stcfr_memo)),'r')


h=h';