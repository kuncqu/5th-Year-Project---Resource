%% Control FOSTRAD_2_0 Aero/Thermo
% Contains the inputs necessary to run FOSTRAD_2_0_Aero_Thermal_Opt
 
clear all
close all
clc
 
% FOSTRAD Code main directory, it must contain also the following folders: 
% "CAD MODELS"             -  folder containing the STL that has to be used
% "external functions"     -  Additional functions required by FOSTRAD
mainDir = '**/DirectoryPath**/FOSTRAD 2.0/';
%directory structure:  
%'/**DirectoryPath**/FOSTRAD 2.0/CAD MODELS'   
%'/**DirectoryPath**/FOSTRAD 2.0/external functions'

%Name of thebb=bb+1; stl file, this is a BINARY format that has been created by SolidWorks,
%different files created with other software may NOT work. 
%There is an additional function to read the ASCII version, but it's not 
%compatible with the default ASCII format provided by SolidWorks.
STLname = 'Ellipsoid_coarse_bin.STL';  % it must be Binary



% Altitude that must be simulated. It uses the NRLMSISE-00 at fixed
% latitude, longitude, day and hour of the year. The atmospheric model has
% a range from 0km to 999km, with 1km spacing, an spline interpolation is
% used to compute atmospheric parameters between each integer altitude.
% In this version the altitude must be a SINGLE NUMBER, more altitude can
% be simulated by increasing the Nsmooth (see the descritpion below)
altitude = 80.451; %[km]  min and max: 0km ~ 999km

% Velocity
Vinf = 7810;   %[m/s]

%reference length (it can be used the normalised diameter of the object)
lref = 4;   % [m]

%reference cross section area, used to compute aerodynamic coefficients
Sref = pi()*1^2;  % [m] %ellipsoid cross section

% beta and alpha must be given as single inputs. A nested for loop may be
% used to obtain multi-attitude database as shown in the example below.
alpha = 0;   % [deg]  angle of attack 
beta = 0;    % [deg]  sideslip angle 


% A Nsmooth number of points must be computed in order to provide
% accurate results. the minimum is 3 points, which is going to provide the
% fastest results. Although, it is suggested to use at least 6~7 points to
% simulate a single altitude.
% The use of a high number Nsmooth points (20~200) may be used to obtain altitude range
% analyses. The "altitude" given as an input will generally be the central simulated
% altitude (unless the spacing algorithm generates an altitude equal to the
% one given as input).
Nsmooth = 10;

%Flags used to activate the different modules. 
AeroFlag = 1;       % if Aeroflag = 1 the Aerodynamics module is active
ThermoFlag = 1;     % if Thermoflag = 1 the Aero-Thermodynamics module is active

% Backface culling algorithm - it provides the .stl triangles list updated
% taking into account the wind direction, removing the unnecessary triangles. 
% BackFaceCulling - Preprocessing phase function provides new matrixes for the triangulation
% F = [A_new(:,1) F_new(:,1:3) N_new(:,1:3)];
% A_new - triangle areas
% F_new - definition of the triangles, list of vertices connected 
% N_new - triangle surface normal vectors

Bfc = BackFaceCulling_opt(Vinf, alpha, beta, mainDir, STLname);
 
% FOSTRAD 2.0 Aero-Thermal-optimized code
% The Aerodynamic and Aero-thermodynamic modules have been integrated into
% one script, which has been converted into a function.
% The script is now using a fixed atmospheric database, which description
% can be found within the code.
% The wall temperature is fixed: Twall = 350K
% The nose radious is fixed: rN = 1m
% In the current function the given outputs are the following

F_Aero = FOSTRAD_2_0_Aero_Thermal_Opt(mainDir, STLname, altitude, Vinf, lref, Sref, alpha, beta, Nsmooth, Bfc, AeroFlag, ThermoFlag);
H(:,1) = F_Aero(:,1);         %- simulated altitudes (considering also the Nsmooth points)
CD(:,1) = F_Aero(:,2);        %- Drag coefficients at different altitudes
CL(:,1) = F_Aero(:,3);        %- Lift coefficients at different altitudes
Qav(:,1) = F_Aero(:,4);       %- Heat transfer flux averaged on the entire surface [W/m^2]
Kn(:,1) = F_Aero(:,5);        %- Knudsen number corresponding to each different altitude


%% Example
% 
% 
% %  Mapping the aerodynamics coefficients at different altitudes and
% %  different attitudes
% 
% 
% % alpha = linspace(0,90,19)';    % [deg]  angle of attack 
% beta = linspace(0,90,10)';      % [deg]  side slip angle 
% alpha = 0;
% 
% % initialising indexes
% a = 1;
% b = 1;
% 
% for a = 1:length(alpha)
%     for b = 1:length(beta)
%     
%         % Initialising the triangulation for the current attitude
%         Bfc = BackFaceCulling_opt(Vinf, alpha(a), beta(b), mainDir, STLname);
%         
%         % Starting the simulation with the current inputs
%         F_Aero = FOSTRAD_2_0_Aero_Thermal_Opt(mainDir, STLname, altitude, Vinf, lref, Sref, alpha(a), beta(b), Nsmooth, Bfc, AeroFlag, ThermoFlag);
% 
%         % Storing the results
%         Htested(:,b) = F_Aero(:,1);
%         CD(:,b) = F_Aero(:,2);
%         CL(:,b) = F_Aero(:,3);
%         Qav(:,b) = F_Aero(:,4);
%         Kn(:,b) = F_Aero(:,5);
% 
%     end
% end
% 
% % plotting the results
% 
% figure()
% hold on
% for i = 1:Nsmooth*2+1
%     plot(beta, CD(i,:),'-o')
%     H_kn_cell{i} = [num2str(round(Htested(i),2)),'km - Kn = ',num2str(Kn(i),'%10.2e\n')];
% 
% end
% title('Aerodynamics')
% xlabel('Sideslip [deg]')
% ylabel('Drag Coefficient')
% 
% 
% legend(H_kn_cell, 'Location', 'NorthWest')
% 
% 


