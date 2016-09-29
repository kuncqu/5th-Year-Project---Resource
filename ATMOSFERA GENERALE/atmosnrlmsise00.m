function [Ti,rhoi] = atmosnrlmsise00(Alt , Lat, Lon, year, yr_day, UT_sec)

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

if Lat~=55 || Lon~=45 || year~=2000 || yr_day~=1 || UT_sec~=1.5
    warning('Warning: the conditions are not consistent with atmosnrlmsise00 database');
end

hDiff = diff(Alt);
if max(hDiff)~=min(hDiff)
    error('Error: The Alt vector shall be evenly spaced')
end
if hDiff(1)~=1e3
    error('Error: The interval between contiguous Alt components must be 1km (i.e. 1000)')
end

Alt = Alt/1e3;
hmin=mix(Alt); %km
minHindex = round(hmin-0.5)+1;
hmax=max(Alt); %km
maxHindex = round(hmax-0.5)+1;
load('NRLMSISE00')

% An array of N-by-9 values of densities (kg/m3 or 1/m3) in selected density units.
% The column order is:
% Density of He, in 1/m3
% Density of O, in 1/m3
% Density of N2, in 1/m3
% Density of O2, in 1/m3
% Density of Ar, in 1/m3
% Total mass density, in kg/m3
% Density of H, in 1/m3
% Density of N, in 1/m3
% Anomalous oxygen number density, in 1/m3
rhoi = nan(length(minHindex:maxHindex),9);
rhoi(:,1) = NRLMSISE00(minHindex:maxHindex,6);
rhoi(:,2:4) = NRLMSISE00(minHindex:maxHindex,2:4);
rhoi(:,5) = NRLMSISE00(minHindex:maxHindex,7);
rhoi(:,7:8) = NRLMSISE00(minHindex:maxHindex,8:9);
rhoi = rhoi*1E6;

% Array of N-by-2 values of temperature, in kelvin.
% The first column is exospheric temperature, in kelvin.
% The second column is temperature at altitude, in kelvin.
Ti = [nan(length(minHindex:maxHindex),1) NRLMSISE00(minHindex:maxHindex,5)];


