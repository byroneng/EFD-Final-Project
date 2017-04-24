%% EFD - Final Project
% Goal~ Analyze data collected at the Sage Brush and Playa sites during the 
% Fall 2012 MATERHORN field campaign using fundamental turbulence methods. 
%
% Hypothesis ~ The impact of the surface roughness increase at the Sage
% brush site compared to the smooth Playa site will impact the boundary
% layer turbulence characteristics in the following ways. 1) Increased
% dissipation 2) increased mixing 3) larger latent heat flux 4)
%
% 1st Mat 2017
% Written by Byron Eng, Matt Moody, & Travis Morrison

%clear vars/close windows & command line
clear all; close all; clc;
%% Load/Preparing Data 

%Tower heights
playa_z = [25.5 19.4 10.4 5.3 2.02 0.61]; %Playa tower sensor heights [m]
sagebrush_z = [18.6 10.15 5.87 2.04 0.55]; %Sagebrush tower sensor heights [m]

%%%%%%%% ONLY LOOKING AT OCT. 14 %%%%%%%%%%%%%%%%%%%%%%%%%
% Playa T, RH, and U for Playa on October 14th in descending height
Data1 = load('./SonicData/Playa Data1415/CSV_1341.Playa_1HZ_2012_10_14_0000.dat');
Data2 = load('./SonicData/Playa Data1415/CSV_1341.Playa_20HZ_2012_10_14_0000.dat');
playa.oct14.Hz1.T = Data1(:,5:2:end-1);
playa.oct14.Hz1.RH = Data1(:,6:2:end);
playa.oct14.Hz20.Ux = rearrangeHeights_playa(Data2(:,5:5:30));
playa.oct14.Hz20.Uy = rearrangeHeights_playa(Data2(:,6:5:31));
playa.oct14.Hz20.Uz = rearrangeHeights_playa(Data2(:,7:5:32));
clear Data1 Data2;

% Sagebrush T, RH, and U for Playa on October 15th in descending height
Data1 = load('./SonicData/Sagebrush Data1415/CSV_6989.Sage_1HZ_2012_10_14_0000.dat');
Data2 = load('./SonicData/Sagebrush Data1415/CSV_6989.Sage_20HZ_2012_10_14_0000.dat');
sagebrush.oct14.Hz1.T = Data1(:,5:2:end-1);
sagebrush.oct14.Hz1.RH = Data1(:,6:2:end);
sagebrush.oct14.Hz20.Ux = rearrangeHeights_sagebrush(Data2(:,5:5:25));
sagebrush.oct14.Hz20.Uy = rearrangeHeights_sagebrush(Data2(:,6:5:26));
sagebrush.oct14.Hz20.Uz = rearrangeHeights_sagebrush(Data2(:,7:5:27));
clear Data1 Data2;
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%rotate the velocities, calculate fluctuations from the mean, wind speed
%and wind direction
N = length(sagebrush.oct14.Hz20.Ux(:,1)); %number of data points
freq = 20; % [Hz]
period = N/freq; %period of interest [seconds]
detrend_chunk = 60*30*20; %number of points we detrend from (30 min)

[ sagebrush.oct14.Hz20.Ux, sagebrush.oct14.Hz20.Uy, sagebrush.oct14.Hz20.Uz,...
    sagebrush.oct14.Hz20.Ux_fluct, sagebrush.oct14.Hz20.Uy_fluct,...
    sagebrush.oct14.Hz20.Uz_fluct, sagebrush.oct14.Hz20.ws_wd ] = ...
    rotate_data( freq,period,N, sagebrush.oct14.Hz20.Ux, sagebrush.oct14.Hz20.Uy, ...
    sagebrush.oct14.Hz20.Uz, length(sagebrush_z), detrend_chunk);
%%
N = length(playa.oct14.Hz20.Ux(:,1));
[ playa.oct14.Hz20.Ux, playa.oct14.Hz20.Uy, playa.oct14.Hz20.Uz,...
    playa.oct14.Hz20.Ux_fluct, playa.oct14.Hz20.Uy_fluct,...
    playa.oct14.Hz20.Uz_fluct, playa.oct14.Hz20.ws_wd ] = ...
    rotate_data( freq,period,N, playa.oct14.Hz20.Ux, playa.oct14.Hz20.Uy,...
    playa.oct14.Hz20.Uz,length(playa_z),detrend_chunk );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Load data on the 18th

% % Playa T, RH, and U for Playa on October 18th in descending height
% Data1 = load('./PlayaDataOct18_2012/CSV_1341.Playa_1HZ_2012_10_18_0000.dat');
% Data2 = load('./PlayaDataOct18_2012/CSV_1341.Playa_20HZ_2012_10_18_0000.dat');
% playa.oct18.Hz1.T = Data1(:,5:2:end-1);
% playa.oct18.Hz1.RH = Data1(:,6:2:end);
% playa.oct18.Hz20.Ux = rearrangeHeights_playa(Data2(:,5:5:30));
% playa.oct18.Hz20.Uy = rearrangeHeights_playa(Data2(:,6:5:31));
% playa.oct18.Hz20.Uz = rearrangeHeights_playa(Data2(:,7:5:32));
% clear Data1 Data2;
% 
% % Sagebrush T, RH, and U for Playa on October 18th in descending height
% Data1 = load('./SagebrushDataOct18_2012/CSV_6989.Sage_1HZ_2012_10_18_0000.dat');
% Data2 = load('./SagebrushDataOct18_2012/CSV_6989.Sage_20HZ_2012_10_18_0000.dat');
% sagebrush.oct18.Hz1.T = Data1(:,5:2:end-1);
% sagebrush.oct18.Hz1.RH = Data1(:,6:2:end);
% sagebrush.oct18.Hz20.Ux = rearrangeHeights_sagebrush(Data2(:,5:5:25));
% sagebrush.oct18.Hz20.Uy = rearrangeHeights_sagebrush(Data2(:,6:5:26));
% sagebrush.oct18.Hz20.Uz = rearrangeHeights_sagebrush(Data2(:,7:5:27));
% clear Data1 Data2;
%%  1.
%Simple Time Averaging: calculate (a) 30-minute averages of u, v, w, and T, as well as (b) ws, wd,
% ?u, ?v, ?w, w0T0,s, u?, Hs, tke, L, and w? (if appropriate), where ws and wd are average wind speed
% and direction, respectively. Describe the stability during the analyzed period.



%% 2.
% Probability Distributions: For a representative 30-minute averaging period, generate a CDF and PDF
% for u, v, w, and T and report the skewness and kurtosis.


%%
% 3. Autocorrelation: Calculate the autocorrelation of at least one 30-minute period. What does this
% indicate?


%% 4. 
% Dissipation: Using Taylor?s frozen turbulence hypothesis, calculate the dissipation rate of turbulent
% kinetic energy for several 30-minute periods. Calculate the Kolmogorov length scale.

%% 5.
%  Turbulence Spectra: For at least one 30-minute averaging period, calculate the following turbulent
% energy spectra: Suu, Svv, Sww, ST T , as well as the following cospectra: Suw, Svw, Swt. What do
% these spectra indicate about the analyzed boundary layer? Is there an inertial subrange?


