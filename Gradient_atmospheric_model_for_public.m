%% Gradient atmospheric model
% 
% Refer to our paper for details on this model:
% Y. Jin and M. Kats, "A gradient model reveals enhanced radiative cooling potential and
%  demonstrates the advantage of broadband emitters," Laser & Photonics Reviews (2025)
%
% 
% Made by Yeonghoon Jin (2024)
% Revised by Yeonghoon Jin (June 08, 2025)
% Please let me know if you find any problems (jin269@wisc.edu, mkats@wisc.edu)


clc;
clear;

% Define parameters
h = 6.62607004e-34;     % m^2 kg/s
k_B = 1.38064852e-23;   % m^2 kg/s^2 K
c = 299792458;          % speed of light: 299 792 459 (m/s)


%% [1] Set a mid-infrared wavelength range (Default = 3 to 25 um)
% Please do not change this range, unless you are familiar with this code.
% The atmospheric transmittance data we have has the same wavelength
% range (3-25 um), so unless you have your own atmosphric transmittance,
% there is no other choice.

lambda_st = 3e-6; % unit: "m". Not um.
lambda_end = 25e-6; 
interval = 10e-9;
lambda = [lambda_st:interval:lambda_end];


%% [2] Choose the location where you want to use the atmospheric transmittance
% So far, we have prepared altitude-dependent atmospheric
% transmittance spectra at 7 locations (and times) on Earth.
% Choose any locations among the bottom location list.
% Copy and paste one, for example, 'Singapore, May 1st 2023';

Location = 'Singapore, May 1st 2023';

% [Location List]
%   Singapore, May 1st 2023
%   LA, Aug 1st 2023
%   New York, Aug 1st 2023
%   Atacama desert, Dec 1st 2023
%   Phoenix, Aug 1st 2023
%   Cairo, Aug 1st 2023
%   Houston, Aug 1st 2023


% We divided the atmosphere into 14 layers.
% If you want to use different number of atmosphric layers,
% you have to edit "atmospheric_Transmittance_Library.mat" first.
altitude_km = [0.5, 1, 1.5, 2, 2.5, 3, 4, 5, 10, 15, 20, 25, 30, 50]; % in km

for t = 1:length(altitude_km)
    altitude = num2str(altitude_km(t));
    T_atm_temp{t} = getAtmTransmittance(Location,altitude,lambda); 
end

% Altitude-dependent atmospheric temperature (in K) at each atmospheric layer.
% These were extracted from NASA PSG, and refer to
% [Y. Jin and M. Kats, A gradient model reveals enhanced radiative cooling potential
% and demonstrates the advantage of broadband emitters, Laser & Photonics Review, 2025]
% for detailed explanation on how these values were extracted.


if Location == "Singapore, May 1st 2023"
    Temp_atm = [301.5, 296, 294, 290, 287, 282, 279, 273,...
        247, 205, 200, 215, 225, 260];    
elseif Location == "LA, Aug 1st 2023"
    Temp_atm = [300.5, 297, 294.5, 292.2, 287.6, 282.5, 277.5, 269.5,...
        242, 211.5, 206.5, 218, 226, 250];
elseif Location == "New York, Aug 1st 2023"
    Temp_atm = [295.3, 289, 284, 280.5, 276, 272.3, 270, 263.5,...
        235, 221, 216.4, 221.5, 226, 255.7];
elseif Location == "Atacama desert, Dec 1st 2023"
    Temp_atm = [288.4, 282.5, 280.8, 276.5, 272.2, 268.1, 264.2, 257,...
        227, 199.2, 208.4, 220.4, 228, 257];
elseif Location == "Phoenix, Aug 1st 2023"
    Temp_atm = [308.5, 303, 296, 290, 285.6, 281, 276.7, 267,...
        242.5, 208.2, 205.7, 218.8, 225.5, 252];
elseif Location == "Cairo, Aug 1st 2023"
    Temp_atm = [306.5, 300, 293, 287.3, 287, 285.7, 282, 273.7,...
        250, 210.8, 202.6, 218.8, 223.4, 250.7];
elseif Location == "Houston, Aug 1st 2023"
    Temp_atm = [304.5, 299, 295.8, 290.8, 286, 281.2, 277.3, 272.6,...
        246.5, 207.4, 207.5, 218.3, 224.6, 248.5];
end

%% [3] Solar absorption
% Import solar AM1.5 Spectrum
fid = fopen('AM_1.5_Global horizontal_based on G197-14_300 nm_to_2500 nm_2201 steps.txt');
dataimport = fscanf(fid, '%g %g', [2 inf]); 
fclose(fid);

% import data
Solar_AM15 = dataimport(2,:);

% Enter the solar irradiance.
% If you want to assume no solar absorption, enter 0 in "Solar_Power".
Solar_Power = 1000; % Unit: W/m2

% Power compensation
Solar_AM15_compensated = Solar_AM15/sum(Solar_AM15) * Solar_Power;

% Import the absorptance spectrum of the device you want to use
% The txt file should be 2 columns [wavelengths, solar absorptance].
% The wavelengths should be from 0.3 to 2.5 um, and the wavelength
% interval must be 0.01 um (10 nm), which is 2201 steps.

fid = fopen('Example_Solar_absorptance_spectrum_from_0.3_to_2.5 um.txt');
dataimport = fscanf(fid, '%g %g', [2 inf]); 
fclose(fid);

Solar_abs_spectrum = dataimport(2,:);

Solar_abs = Solar_abs_spectrum .* Solar_AM15_compensated;
Solar_abs_power = sum(Solar_abs);


%% [4] Set the emitter temperature range you want to calculate
Temp_emitter = 273 + [10:1:35]; % Type a temperature range in degC
aoi_interval = 5;             % Interval of angle of incidence 
aoi = [0:aoi_interval:85];    % Angle of incidence from 0 to 85 deg


%% [5-1] If you want to choose ideal thermal emitters, which are defined right below, 
%  you can choose the emissivity here, and you can ignore [5-2] below.
%  But if you want to import your own emissvity, go to [5-2] directly.

% 1: Broadband emitter (BE) - emissivity of 1 in the mid-infrared (e.g.,3-25 um)
% 2: Selective emitter (SE) - emissivity of 1 at 8-13 um, and 0 at elsewhere

% Selective emitter wavelength range (8-13 um)
lam_st = int16((8e-6 - lambda_st)/interval + 1);
lam_end = int16((13e-6 - lambda_st)/interval + 1);


% Emissivity of the ideal BE
% The emissivity is 1 for all incident/emission angles
e_BE = 1 * ones(length(aoi),length(lambda)); 

% Emissivity of the ideal SE
e_SE = 0 * ones(length(aoi),length(lambda));

for j = 1:length(aoi)
    for i = lam_st:lam_end
        e_SE(j,i) = 1;
    end
end

% Choose the emissivity ("e_BE" or "e_SE")
e_sample = e_BE;


%% [5-2] Import your own emissivity

% [Case 1] Assume the emissivity is not angle-dependent.
% This assumption allows you to roughly estimate the cooling performance of
% your sample, if you have an emissivity at normal incident angle only.
% 
% (1) import your emissivity txt file, and (2) enter 0 in the angle_dependency.

% The txt file should be 2 columns [wavelengths, emissivity].
% The wavelength range should be 3-25 um in 0.01 um intervals (2201 steps).

angle_dependency = 0; % 0 for Case 1, 1 for Case 2.

% Assuming that the imported emissivity applies to all incident/emission angles.
% This means that your imported emissivity is maintained at all angles.
for j = 1:length(aoi)
    if angle_dependency == 0
        fid = fopen('Example_Emissivity_spectrum_from_3_to_25 um_10 nm step.txt');
        dataimport = fscanf(fid, '%g', [2 inf]);
        fclose(fid);
        
        e_sample(j,:) = dataimport(2,:);
    else
        
% [Case 2] Emissivity is angle-dependent.
% If you want to consider the angle-dependency, (1) enter 1 in the angle_dependency.
% (2) Prepare the angle-dependent emissivity in a single txt file.
% Please refer to an example txt file:
% I set the angle interval of 5 deg, and consider the angle
% from 0 deg (normal incident) to 85 deg, for example,
% [wavelengths, 0 deg, 5 deg, ..., 85 deg] -> 19 columns.

        fid = fopen('Example_angle_dependent_Emissivity_0_to_85_5_deg step.txt');
        dataimport = fscanf(fid, '%g', [19 inf]);
        fclose(fid);

        % The 1st column is wavelength so we don't need to import this.    
        for i = 2:19
            e_sample(j,:) =  dataimport(i,:);
        end
    end
end

%% [6] Set the non-radiative heat transfer coefficient (hc)
h_c = 6; % unit: W/m2/K

for nn = 1:length(Temp_emitter)
        
    %% Calculation of Patm: The power radiated from the atmosphere and absorbed by the emitter    
    for pp = 1:length(Temp_atm)
        I_BB_atm = 2*h*c^2./lambda.^5 * 1./(exp((h*c)./(lambda * k_B * Temp_atm(pp))) - 1); % unit: W/m^2/m/sr
        
        if pp == 1 % The first atmospheric layer from the ground             
            T_atm = T_atm_temp{pp}; % Atmospheric transmittance
            
            for i = 1:length(aoi)
                e_atm = 1 - (T_atm).^(1/cosd(aoi(i)));
                P_atm_temp = 2 * pi * sind(aoi(i)) * cosd(aoi(i)) * I_BB_atm .* e_sample(i,:) .* e_atm;
                P_atm_tempC(i) = sum(P_atm_temp) * interval; % unit: W/m^2/sr
            end

            PP_atm(pp) = sum(P_atm_tempC) * aoi_interval *pi/2/90; % unit: W/m2                

        else
            % (pp)th atmospheric layer
            T_atm = T_atm_temp{pp};

            for i = 1:length(aoi) 
                e_atm = 1 - (T_atm).^(1/cosd(aoi(i)));
                P_atm_temp = 2 * pi * sind(aoi(i)) * cosd(aoi(i)) * I_BB_atm .* e_sample(i,:) .* e_atm;                
                P_atm_tempC(i) = sum(P_atm_temp) * interval; % unit: W/m^2/sr
            end

            % (pp-1)th atmospheric layer          
            T_atm = T_atm_temp{pp-1};
            for i = 1:length(aoi) 
                e_atm = 1 - (T_atm).^(1/cosd(aoi(i)));
                P_atm_temp = 2 * pi * sind(aoi(i)) * cosd(aoi(i)) * I_BB_atm .* e_sample(i,:) .* e_atm;                
                P_atm_tempP(i) = sum(P_atm_temp) * interval; % unit: W/m^2/sr
            end
            PP_atm(pp) = (sum(P_atm_tempC) - sum(P_atm_tempP)) * aoi_interval *pi/2/90; % unit: W/m2
        end
    end

    P_atm(nn) = sum(PP_atm);

    %% Calculation of P_rad: Power radiated from the emitter
    I_BB = 2*h*c^2./lambda.^5 * 1./(exp((h*c)./(lambda * k_B * Temp_emitter(nn))) - 1); % unit: W/m^2/sr/m

    for i = 1:length(aoi)
        P_rad_temp = 2 * pi * sind(aoi(i)) * cosd(aoi(i)) * I_BB .* e_sample(i,:); 
        P_rad_temp2(i) = sum(P_rad_temp) * interval; % Intergral as a function of the wavelength (W/m^2/sr)
    end

    P_rad(nn) = sum(P_rad_temp2) * aoi_interval *pi/2/90; % Intergral as a function of the aoi (W/m^2)

    %% Calculation of P_sun (already calculated at the top)
    P_sun(nn) = Solar_abs_power;   % Unit: W/m^2

    %% Calculation of non-radiative heat exchange (e.g. convection)    
    P_con(nn) = h_c * (Temp_atm(1) - Temp_emitter(nn)); % "Temp_atm(1)" is the ambient temperature.

    %% Net cooling power
    P_cool(nn) = P_rad(nn) - P_sun(nn) - P_atm(nn) - P_con(nn); % unit: W/m2
end

figure('color','w');
plot(Temp_emitter-273,P_cool);
xlabel('Emitter temperature in C');
ylabel('Net cooling power (W/m2)');
