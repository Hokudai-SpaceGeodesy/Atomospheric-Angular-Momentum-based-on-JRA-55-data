function MakeChiP_IB_06hr(Year, Month, Day, Hour)
% MakeChiP_IB_06hr.m
% R4/07/11: MakeChiW_06hrを改造
% R5/01/11: コメント整理

%%
% Load Data
YMDH = sprintf('%4u%02u%02u%02u', Year, Month, Day, Hour);
InFile = './tmpP/anl_surf125.' + string(YMDH) + '.txt';

LSfile = './data/LandSea.mat';

Fid = fopen(InFile, 'r');
Data = fscanf(Fid, '%f %f %f', [3, inf]);
Data = Data.'; % [Latitude Longtitude P[Pa]]
fclose('all');

%%
% Set Constants from Eubanks et al., 1993.
R = 6.371 * 10 ^ 6; % Earth's mean radius [m]
G =  9.81; % Gravity acceleration [m s^-2]
% Omega = 7.292115 * 10 ^ (-5); % [rad s^-1]
C_A = 2.610 * 10 ^ 35; % of the entire Earth [kg m^2]
Cm = 7.1236 * 10 ^ 37; % of the Mantle [kg m^2]

%Constants for Wind from Eubanks et al., 1993.
Const_XY = -1.098;
Const_Z = 0.753;

%%
% Prepare Calculation
Lats = Data(:,1); % Phi [deg]
Lons = Data(:,2); % Lambda [deg]
P = Data(:,3); % [Pa]

% Load Land-Sea Ratio Data
LSdat = load(LSfile, 'LandSea');
LandSea = LSdat.LandSea; clear LSdat; % Ratio (1:Land 0:Sea)

P = reshape(P, [288. 145]); % Lon * Lat  の2次元配列
Lats = reshape(Lats, [288, 145]);
Lons = reshape(Lons, [288,145]);
LandSea = reshape(LandSea, [288, 145]);

%%
% Inverted Barometer Process
% 緯度+-90degは0になるので、積分は不要

P_Ocean = P .* (1 - LandSea);
P_Ocean_mean = sum(P_Ocean .* cosd(Lats), 'all') ./ sum((1 - LandSea) .* cosd(Lats), 'all');

P = P_Ocean_mean .* (1 - LandSea) + P .* LandSea;

%%
% Calculate Chi_W

F0 = P .* sind(Lats) .* cosd(Lats) .^ 2 .* exp(1i .* deg2rad(Lons));
F0Z = P .* cosd(Lats) .^ 3;

F1 = sum(F0 .* deg2rad(1.25), 1);
F1Z = sum(F0Z .* deg2rad(1.25), 1); %経度について積分

F1(:,1) = F1(:,1) ./ 2;
F1(:,145) = F1(:,145) ./ 2;
F1Z(:,1) = F1Z(:,1) ./ 2;
F1Z(:,145) = F1Z(:,145) ./ 2;

F2 = sum(F1 .* deg2rad(1.25), 2); %緯度について積分
F2Z = sum(F1Z .* deg2rad(1.25), 2);


Chi_rad = Const_XY * R ^ 4 * F2 / (C_A * G); 
Chi = Chi_rad * (180 * 3600) * 1000 / pi; % [rad] -> [mas]

ChiZ_rad = Const_Z * R ^ 4 * F2Z / (Cm * G);
ChiZ = ChiZ_rad * 24 * 3600 * 1000; % [ms] 

%%
%Save Chi

OutFile = './data/ChiP_IB.txt';
MJD = juliandate(datetime(Year, Month, Day), 'modifiedjuliandate');

Fid_out = fopen(OutFile, 'a');
fprintf(Fid_out, '%d %02d %f %f %f\n', MJD, Hour, real(Chi), imag(Chi), ChiZ);
fclose(Fid_out);