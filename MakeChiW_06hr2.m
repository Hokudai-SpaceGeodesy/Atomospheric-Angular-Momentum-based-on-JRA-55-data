function MakeChiW_06hr2(Year, Month, Day, Hour)
% MaleChiW_06hr2.m
% R4/10/31: 地形を考慮する計算の為にMakeChiW_06hr.mを改造
% R4/11/09: 試運転
% R5/01/11: コメントのみ整理


%%
% Load Data
YMDH = sprintf('%4u%02u%02u%02u', Year, Month, Day, Hour);
InFile_u = "./tmpW/anl_p125_ugrd." + string(YMDH) + ".txt";
InFile_v = "./tmpW/anl_p125_vgrd." + string(YMDH) + ".txt";
InFile_surf = "./tmpP/anl_surf125." + string(YMDH) + ".txt";

Fid_u = fopen(InFile_u, 'r');
Fid_v = fopen(InFile_v, 'r');
Fid_surf = fopen(InFile_surf, "r");
Data_u = fscanf(Fid_u, '%f %f %f %f', [4, inf]);
Data_v = fscanf(Fid_v, '%f %f %f %f', [4, inf]);
Data_surf = fscanf(Fid_surf, "%f %f %f %f %f", [5, inf]);
Data_u = Data_u.'; % [Pressure Latitude Longtitude U[m s^-1]]
Data_v = Data_v.'; % [Pressure Latitude Longtitude V[m s^-1]]
Data_surf = Data_surf.';
fclose('all');

%%
% Set Constants from Eubanks et al., 1993.
R = 6.371 * 10 ^ 6; % Earth's mean radius [m]
G =  9.81; % Gravity acceleration [m s^-2]
Omega = 7.292115 * 10 ^ (-5); % [rad s^-1]
C_A = 2.610 * 10 ^ 35; % of the entire Earth [kg m^2]
Cm = 7.1236 * 10 ^ 37; % of the Mantle [kg m^2]

%Constants for Wind from Eubanks et al., 1993.
Const_XY = -1.5913;
Const_Z = 0.998;

%%
% Prepare Calculation
Ps = Data_u(:,1) .* 100; %[hPa] -> [Pa]
Lats = Data_u(:,2); % Phi
Lons = Data_u(:,3); % Lambda
U = Data_u(:,4);
V = Data_v(:,4);
P_surf = Data_surf(:, 3);
U_surf = Data_surf(:, 4);
V_surf = Data_surf(:, 5);

% Reshape Data
Ps = reshape(Ps, [288, 145, 37]); % Lon * Lat * P の3次元配列
Lats = reshape(Lats, [288, 145, 37]);
Lons = reshape(Lons, [288, 145, 37]);
U = reshape(U, [288, 145, 37]);
V = reshape(V, [288, 145, 37]);
P_surf = reshape(P_surf, [288, 145]);
U_surf = reshape(U_surf, [288, 145]);
V_surf = reshape(V_surf, [288, 145]);

clear Data_u Data_v Data_surf
%%
% Calculate Chi_W

F0 = (U .* sind(Lats) + 1i .* V) .* cosd(Lats) .* exp(1i .* deg2rad(Lons)); 
F0Z = (U .* cosd(Lats) .^ 2);
F0_surf = (U_surf .* sind(Lats(:, :, 1)) + 1i .* V_surf) .* cosd(Lats(:, :, 1)) .* exp(1i .* deg2rad(Lons(:, :, 1))); 
F0Z_surf = (U_surf .* cosd(Lats(:, :, 1)) .^ 2);

Pmin = zeros([288, 145]);
F1 = zeros([288, 145]);
F1Z = zeros([288, 145]);

% 以下高度についての積分
for J = 1:288
     for K = 1:145
         for I = 37:-1:1
             if Ps(J, K, I) < P_surf(J, K) && Pmin(J, K) == 0
                Pmin(J, K) = I;
                break
             end
         end

        N = Pmin(J, K) + 1;
        Ps(J, K, N) = P_surf(J, K);
        if N == 1
            F1(J, K) = 0;
            F1Z(J, K) = 0;
        elseif N == 2
            F1(J, K) = (F0(J, K) + F0_surf(J, K)) * (P_surf(J, K) - Ps(J, K, 1)) ./ 2;
            F1Z(J, K) = (F0Z(J, K) + F0Z_surf(J, K)) * (P_surf(J, K) - Ps(J, K, 1)) ./ 2;
        else % N = 3:38
            DeltaP = zeros(N, 1); % Barnes et al., 1983 の足し合わせ
            DeltaP(1) = Ps(J, K, 2) - Ps(J, K, 1);
            DeltaP(N) = Ps(J, K, N) - Ps(J, K, N-1);
            DeltaP(2:N-1) = Ps(J, K, 3:N) - Ps(J, K, 1:N-2);
            F01 = zeros(N, 1); F01Z = zeros(N, 1);
            F01(1:N-1) = F0(J, K, 1:N-1);
            F01Z(1:N-1)  = F0Z(J, K, 1:N-1);
            F01(N) = F0_surf(J, K);
            F01Z(N) = F0_surf(J, K);
            F1(J, K) = sum(F01 .* DeltaP(:, 1) ./ 2, "all"); 
            F1Z(J, K) = sum(F01Z .* DeltaP(:, 1) ./ 2, "all"); 
        end
     end
end
% ここまで

F1(1, :) = F1(1, :) ./ 2;
F1(288, :) = F1(288, :) ./ 2;
F1Z(1, :) = F1Z(1, :) ./ 2;
F1Z(288, :) = F1Z(288, :) ./ 2;

F2 = sum(F1 .* deg2rad(1.25), 1);
F2Z = sum(F1Z .* deg2rad(1.25), 1); % 経度について積分

F2(:, 1) = F2(:, 1) ./ 2;
F2(:, 145) = F2(:, 145) ./ 2;
F2Z(:, 1) = F2Z(:, 1) ./ 2;
F2Z(:, 145) = F2Z(:, 145) ./ 2;

F3 = sum(F2 .* deg2rad(1.25), 2);
F3Z = sum(F2Z .* deg2rad(1.25), 2); % 緯度について積分


%%
% スケール

Chi_rad = Const_XY * R ^ 3 * F3 / (C_A * Omega * G); 
Chi = Chi_rad * (180 * 3600) * 1000 / pi; % [rad] -> [mas]

ChiZ_rad = Const_Z * R ^ 3 * F3Z / (Cm * Omega * G);
ChiZ = ChiZ_rad * 24 * 3600 * 1000; % [ms] 

%%
%Save Chi

OutFile = './data/ChiW2.txt';
% MJD = juliandate(datetime(Year, Month, Day), 'modifiedjuliandate');

Fid_out = fopen(OutFile, 'a');
fprintf(Fid_out, '%d %02d %02d %02d %f %f %f\n', Year, Month, Day, Hour, real(Chi), imag(Chi), ChiZ);
fclose(Fid_out);