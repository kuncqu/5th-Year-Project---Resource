%% atmospheric database
load('NRLMSISE00')
% NRLMSISE00 condition for:
% Year= 2000, Month=  1, Day=  1, Hour= 1.50,
% Time_type = Universal
% Coordinate_type = Geographic
% Latitude=   55.00, Longitude=   45.00,
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

hmin=50; %km
minHindex = round(hmin-0.5)+1;
hmax=180; %km
maxHindex = round(hmax-0.5)+1;
load('NRLMSISE00')


%He O N2 O2 Ar H N [1/m^3] gas densities
rhoi = [NRLMSISE00(minHindex:maxHindex,6) NRLMSISE00(minHindex:maxHindex,2:4) NRLMSISE00(minHindex:maxHindex,7:9)]*1E6;
Ti = NRLMSISE00(minHindex:maxHindex,5);


