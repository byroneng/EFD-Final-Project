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

%Options:
reloaddata = true; %Reload all the data. Set this to true for first run
recalc = true; %Recalculate things like PDFs, CDFs, etc.

if reloaddata
%clear vars/close windows & command line
clear all; close all; clc;
%% Load/Preparing Data 

%Tower heights
% playa_z = [25.5 19.4 10.4 5.3 2.02 0.61]; %Playa tower sensor heights [m]
% sagebrush_z = [18.6 10.15 5.87 2.04 0.55]; %Sagebrush tower sensor heights [m]

%%%%%%%% ONLY LOOKING AT OCT. 14 %%%%%%%%%%%%%%%%%%%%%%%%%
% Playa T, RH, and U for Playa on October 14th in descending height
% Data1 = load('./SonicData/Playa Data1415/CSV_1341.Playa_1HZ_2012_10_14_0000.dat');
% Data2 = load('./SonicData/Playa Data1415/CSV_1341.Playa_20HZ_2012_10_14_0000.dat');
% t =reshape(Data1(:,1:4),[],4);
% 
% playa.oct14.Hz1.T = Data1(1:end,5:2:end-1);
% playa.oct14.Hz1.RH = Data1(:,6:2:end);
% playa.oct14.Hz20.Ux = rearrangeHeights_playa(Data2(1:end,5:5:30));
% playa.oct14.Hz20.Uy = rearrangeHeights_playa(Data2(:,6:5:31));
% playa.oct14.Hz20.Uz = rearrangeHeights_playa(Data2(:,7:5:32));
% clear Data1 Data2;
% 
% % Sagebrush T, RH, and U for Playa on October 15th in descending height
% Data1 = load('./SonicData/Sagebrush Data1415/CSV_6989.Sage_1HZ_2012_10_14_0000.dat');
% Data2 = load('./SonicData/Sagebrush Data1415/CSV_6989.Sage_20HZ_2012_10_14_0000.dat');
% sagebrush.oct14.Hz1.T = Data1(:,5:2:end-1);
% sagebrush.oct14.Hz1.RH = Data1(:,6:2:end);
% sagebrush.oct14.Hz20.Ux = rearrangeHeights_sagebrush(Data2(:,5:5:25));
% sagebrush.oct14.Hz20.Uy = rearrangeHeights_sagebrush(Data2(:,6:5:26));
% sagebrush.oct14.Hz20.Uz = rearrangeHeights_sagebrush(Data2(:,7:5:27));
% clear Data1 Data2;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%rotate the velocities, calculate fluctuations from the mean, wind speed
%and wind direction
% N = length(sagebrush.oct14.Hz20.Ux(:,1)); %number of data points
% freq = 20; % [Hz]
% period = N/freq; %period of interest [seconds]
% detrend_chunk = 60*30*20; %number of points we detrend from (30 min)
% 
% [ sagebrush.oct14.Hz20.Ux, sagebrush.oct14.Hz20.Uy, sagebrush.oct14.Hz20.Uz,...
%     sagebrush.oct14.Hz20.Ux_fluct, sagebrush.oct14.Hz20.Uy_fluct,...
%     sagebrush.oct14.Hz20.Uz_fluct, sagebrush.oct14.Hz20.ws_wd ] = ...
%     rotate_data( freq,period,N, sagebrush.oct14.Hz20.Ux, sagebrush.oct14.Hz20.Uy, ...
%     sagebrush.oct14.Hz20.Uz, length(sagebrush_z), detrend_chunk);
% 
% N = length(playa.oct14.Hz20.Ux(:,1));
% [ playa.oct14.Hz20.Ux, playa.oct14.Hz20.Uy, playa.oct14.Hz20.Uz,...
%     playa.oct14.Hz20.Ux_fluct, playa.oct14.Hz20.Uy_fluct,...
%     playa.oct14.Hz20.Uz_fluct, playa.oct14.Hz20.ws_wd ] = ...
%     rotate_data( freq,period,N, playa.oct14.Hz20.Ux, playa.oct14.Hz20.Uy,...
%     playa.oct14.Hz20.Uz,length(playa_z),detrend_chunk );
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Tower heights
playa_z = [25.5 19.4 10.4 5.3 2.02 0.61]; %Playa tower sensor heights [m]
sagebrush_z = [18.6 10.15 5.87 2.04 0.55]; %Sagebrush tower sensor heights [m]

%Load data on the 18th
% Playa T, RH, and U for Playa on October 18th in descending height
Data1 = load('./PlayaDataOct18_2012/CSV_1341.Playa_1HZ_2012_10_18_0000.dat');
Data2 = load('./PlayaDataOct18_2012/CSV_1341.Playa_20HZ_2012_10_18_0000.dat');
playa.oct18.Hz1.T = Data1(:,5:2:end-1);
playa.oct18.Hz1.RH = Data1(:,6:2:end);
playa.oct18.Hz20.Ux = rearrangeHeights_playa(Data2(:,5:5:30));
playa.oct18.Hz20.Uy = rearrangeHeights_playa(Data2(:,6:5:31));
playa.oct18.Hz20.Uz = rearrangeHeights_playa(Data2(:,7:5:32));
clear Data1 Data2;

% Sagebrush T, RH, and U for Playa on October 18th in descending height
Data1 = load('./SagebrushDataOct18_2012/CSV_6989.Sage_1HZ_2012_10_18_0000.dat');
Data2 = load('./SagebrushDataOct18_2012/CSV_6989.Sage_20HZ_2012_10_18_0000.dat');
sagebrush.oct18.Hz1.T = Data1(:,5:2:end-1);
sagebrush.oct18.Hz1.RH = Data1(:,6:2:end);
sagebrush.oct18.Hz20.Ux = rearrangeHeights_sagebrush(Data2(:,5:5:25));
sagebrush.oct18.Hz20.Uy = rearrangeHeights_sagebrush(Data2(:,6:5:26));
sagebrush.oct18.Hz20.Uz = rearrangeHeights_sagebrush(Data2(:,7:5:27));
clear Data1 Data2;
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%rotate the velocities, calculate fluctuations from the mean, wind speed
%and wind direction
num_chuck = 42; 

N = length(sagebrush.oct18.Hz20.Ux(:,1)); %number of data points
freq = 20; % [Hz]
period = N/freq; %period of interest [seconds]
detrend_chunk = 60*30*20; %number of points we detrend from (30 min)

[ sagebrush.oct18.Hz20.Ux, sagebrush.oct18.Hz20.Uy, sagebrush.oct18.Hz20.Uz,...
    sagebrush.oct18.Hz20.Ux_fluct, sagebrush.oct18.Hz20.Uy_fluct,...
    sagebrush.oct18.Hz20.Uz_fluct, sagebrush.oct18.Hz20.ws_wd ] = ...
    rotate_data( freq,period,N, sagebrush.oct18.Hz20.Ux, sagebrush.oct18.Hz20.Uy, ...
    sagebrush.oct18.Hz20.Uz, length(sagebrush_z), detrend_chunk);

N = length(playa.oct18.Hz20.Ux(:,1));

[ playa.oct18.Hz20.Ux, playa.oct18.Hz20.Uy, playa.oct18.Hz20.Uz,...
    playa.oct18.Hz20.Ux_fluct, playa.oct18.Hz20.Uy_fluct,...
    playa.oct18.Hz20.Uz_fluct, playa.oct18.Hz20.ws_wd ] = ...
    rotate_data( freq,period,N, playa.oct18.Hz20.Ux, playa.oct18.Hz20.Uy,...
    playa.oct18.Hz20.Uz,length(playa_z),detrend_chunk );

<<<<<<< Updated upstream
end %if reloaddata
    
=======
>>>>>>> Stashed changes
%%  1.
%Simple Time Averaging: calculate (a) 30-minute averages of u, v, w, and T, as well as (b) ws, wd,
% ?u, ?v, ?w, w0T0,s, u?, Hs, tke, L, and w? (if appropriate), where ws and wd are average wind speed
% and direction, respectively. Describe the stability during the analyzed period.



%% 2.
% Probability Distributions: For a representative 30-minute averaging period, generate a CDF and PDF
% for u, v, w, and T and report the skewness and kurtosis.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           Calculate PDF and CDF
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if recalc

chunk = 30*60*1;
t_start = chunk*num_chuck;
t_end = t_start+chunk;
N_1hz = chunk;
%std
sagebrush.oct18.Hz1.T_std = std(sagebrush.oct18.Hz1.T(t_start:t_end,:),0,1,'omitnan');
playa.oct18.Hz1.T_std = std(playa.oct18.Hz1.T(t_start:t_end,:),0,1,'omitnan');
%mean
sagebrush.oct18.Hz1.T_mean = mean(sagebrush.oct18.Hz1.T(t_start:t_end,:),1,'omitnan');
playa.oct18.Hz1.T_mean= mean(playa.oct18.Hz1.T(t_start:t_end,:),1,'omitnan');
%compute pdf and cdf
for z = 1:length(sagebrush_z)
    [Nx,sagebrush.oct18.Hz1.T_edges(:,z)]=histcounts(sagebrush.oct18.Hz1.T(:,z),N_1hz,'Normalization','pdf');
    [sagebrush.oct18.Hz1.T_pdf(:,z)] = PDF(sagebrush.oct18.Hz1.T_std(z),...
        sagebrush.oct18.Hz1.T_mean(z),sagebrush.oct18.Hz1.T_edges(:,z));
    sagebrush.oct18.Hz1.T_cdf(:,z) = .5.*(1+erf((sagebrush.oct18.Hz1.T_edges(:,z)- ...
        sagebrush.oct18.Hz1.T_mean(z))./(sagebrush.oct18.Hz1.T_std(z)*2^.5)));
end
for z = 1:length(playa_z)
    [Nx,playa.oct18.Hz1.T_edges(:,z)]=histcounts(playa.oct18.Hz1.T(:,z),N_1hz,'Normalization','pdf');
    [playa.oct18.Hz1.T_pdf(:,z)] = PDF(playa.oct18.Hz1.T_std(z),...
        playa.oct18.Hz1.T_mean(z),playa.oct18.Hz1.T_edges(:,z));
    playa.oct18.Hz1.T_cdf(:,z) = .5.*(1+erf((playa.oct18.Hz1.T_edges(:,z)- ...
        playa.oct18.Hz1.T_mean(z))./(playa.oct18.Hz1.T_std(z)*2^.5)));
end

chunk = 30*60*20;
t_start = chunk*num_chuck;
t_end = t_start+chunk;
N_20hz = chunk;
%std
sagebrush.oct18.Hz20.Ux_std(:) = std(sagebrush.oct18.Hz20.Ux(t_start:t_start+chunk,:),0,1,'omitnan');
sagebrush.oct18.Hz20.Uy_std(:) = std(sagebrush.oct18.Hz20.Uy(t_start:t_start+chunk,:),0,1,'omitnan');
sagebrush.oct18.Hz20.Uz_std(:) = std(sagebrush.oct18.Hz20.Uz(t_start:t_start+chunk,:),0,1,'omitnan');
playa.oct18.Hz20.Ux_std(:) = std(playa.oct18.Hz20.Ux(t_start:t_start+chunk,:),0,1,'omitnan');
playa.oct18.Hz20.Uy_std(:) = std(playa.oct18.Hz20.Uy(t_start:t_start+chunk,:),0,1,'omitnan');
playa.oct18.Hz20.Uz_std(:) = std(playa.oct18.Hz20.Uz(t_start:t_start+chunk,:),0,1,'omitnan');
%mean
sagebrush.oct18.Hz20.Ux_mean(:) = mean(sagebrush.oct18.Hz20.Ux(t_start:t_start+chunk,:),1,'omitnan');
sagebrush.oct18.Hz20.Uy_mean(:) = mean(sagebrush.oct18.Hz20.Uy(t_start:t_start+chunk,:),1,'omitnan');
sagebrush.oct18.Hz20.Uz_mean(:) = mean(sagebrush.oct18.Hz20.Uz(t_start:t_start+chunk,:),1,'omitnan');
playa.oct18.Hz20.Ux_mean(:) = mean(playa.oct18.Hz20.Ux(t_start:t_start+chunk,:),1,'omitnan');
playa.oct18.Hz20.Uy_mean(:) = mean(playa.oct18.Hz20.Uy(t_start:t_start+chunk,:),1,'omitnan');
playa.oct18.Hz20.Uz_mean(:) = mean(playa.oct18.Hz20.Uz(t_start:t_start+chunk,:),1,'omitnan');

%compute pdf and cdf
for z = 1:length(sagebrush_z)
    %U
    [Nx,sagebrush.oct18.Hz20.Ux_edges(:,z)]=histcounts(sagebrush.oct18.Hz20.Ux(:,z),N_20hz,'Normalization','pdf');
    [sagebrush.oct18.Hz20.Ux_pdf(:,z)] = PDF(sagebrush.oct18.Hz20.Ux_std(z),...
        sagebrush.oct18.Hz20.Ux_mean(z),sagebrush.oct18.Hz20.Ux_edges(:,z));
    sagebrush.oct18.Hz20.Ux_cdf(:,z) = .5.*(1+erf((sagebrush.oct18.Hz20.Ux_edges(:,z)- ...
        sagebrush.oct18.Hz20.Ux_mean(z))./(sagebrush.oct18.Hz20.Ux_std(z)*2^.5)));
    %V
    [Nx,sagebrush.oct18.Hz20.Uy_edges(:,z)]=histcounts(sagebrush.oct18.Hz20.Uy(:,z),N_20hz,'Normalization','pdf');
    [sagebrush.oct18.Hz20.Uy_pdf(:,z)] = PDF(sagebrush.oct18.Hz20.Uy_std(z),...
        sagebrush.oct18.Hz20.Uy_mean(z),sagebrush.oct18.Hz20.Uy_edges(:,z));
    sagebrush.oct18.Hz20.Uy_cdf(:,z) = .5.*(1+erf((sagebrush.oct18.Hz20.Uy_edges(:,z)- ...
        sagebrush.oct18.Hz20.Uy_mean(z))./(sagebrush.oct18.Hz20.Uy_std(z)*2^.5)));
    %W
    [Nx,sagebrush.oct18.Hz20.Uz_edges(:,z)]=histcounts(sagebrush.oct18.Hz20.Uz(:,z),N_20hz,'Normalization','pdf');
    [sagebrush.oct18.Hz20.Uz_pdf(:,z)] = PDF(sagebrush.oct18.Hz20.Uz_std(z),...
        sagebrush.oct18.Hz20.Uz_mean(z),sagebrush.oct18.Hz20.Uz_edges(:,z));
    sagebrush.oct18.Hz20.Uz_cdf(:,z) = .5.*(1+erf((sagebrush.oct18.Hz20.Uz_edges(:,z)- ...
        sagebrush.oct18.Hz20.Uz_mean(z))./(sagebrush.oct18.Hz20.Uz_std(z)*2^.5)));
    
end

for z = 1:length(playa_z) %U
    [Nx,playa.oct18.Hz20.Ux_edges(:,z)]=histcounts(playa.oct18.Hz20.Ux(:,z),N_20hz,'Normalization','pdf');
    [playa.oct18.Hz20.Ux_pdf(:,z)] = PDF(playa.oct18.Hz20.Ux_std(z),...
        playa.oct18.Hz20.Ux_mean(z),playa.oct18.Hz20.Ux_edges(:,z));
    playa.oct18.Hz20.Ux_cdf(:,z) = .5.*(1+erf((playa.oct18.Hz20.Ux_edges(:,z)- ...
        playa.oct18.Hz20.Ux_mean(z))./(playa.oct18.Hz20.Ux_std(z)*2^.5)));
    %V
    [Nx,playa.oct18.Hz20.Uy_edges(:,z)]=histcounts(playa.oct18.Hz20.Uy(:,z),N_20hz,'Normalization','pdf');
    [playa.oct18.Hz20.Uy_pdf(:,z)] = PDF(playa.oct18.Hz20.Uy_std(z),...
        playa.oct18.Hz20.Uy_mean(z),playa.oct18.Hz20.Uy_edges(:,z));
    playa.oct18.Hz20.Uy_cdf(:,z) = .5.*(1+erf((playa.oct18.Hz20.Uy_edges(:,z)- ...
        playa.oct18.Hz20.Uy_mean(z))./(playa.oct18.Hz20.Uy_std(z)*2^.5)));
    %W
    [Nx,playa.oct18.Hz20.Uz_edges(:,z)]=histcounts(playa.oct18.Hz20.Uz(:,z),N_20hz,'Normalization','pdf');
    [playa.oct18.Hz20.Uz_pdf(:,z)] = PDF(playa.oct18.Hz20.Uz_std(z),...
        playa.oct18.Hz20.Uz_mean(z),playa.oct18.Hz20.Uz_edges(:,z));
    playa.oct18.Hz20.Uz_cdf(:,z) = .5.*(1+erf((playa.oct18.Hz20.Uz_edges(:,z)- ...
        playa.oct18.Hz20.Uz_mean(z))./(playa.oct18.Hz20.Uz_std(z)*2^.5)));
end

end %if recalc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                         Plot CDF and PDF
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure()
subplot(4,2,1)
for z = 1:length(sagebrush_z)
    plot(sagebrush.oct18.Hz20.Ux_edges(:,z),sagebrush.oct18.Hz20.Ux_pdf(:,z));
    hold on
end
set(gca,'fontsize', 15)
ylabel('PDF $u$', 'interpreter','latex','fontsize',20);
legend('18.6 m','10.15','5.87','2.04','0.55')
title('Sagebrush', 'interpreter','latex','fontsize',25)
xlabel('u (ms$^{-1}$)', 'interpreter','latex','fontsize',20);


subplot(4,2,2)
for z = 1:length(playa_z)
    plot(playa.oct18.Hz20.Ux_edges(:,z),playa.oct18.Hz20.Ux_pdf(:,z));
    hold on
end
set(gca,'fontsize', 15)
title('Playa', 'interpreter','latex','fontsize',25)
legend('25.5 m','19.4','10.4','5.3','2.02','0.61')
xlabel('u (ms$^{-1}$)', 'interpreter','latex','fontsize',20);

subplot(4,2,3)
for z = 1:length(sagebrush_z)
    plot(sagebrush.oct18.Hz20.Uy_edges(:,z),sagebrush.oct18.Hz20.Uy_pdf(:,z));
    hold on
end
set(gca,'fontsize', 15)
ylabel('PDF $v$', 'interpreter','latex','fontsize',20);
legend('18.6 m','10.15','5.87','2.04','0.55')
xlabel('v (ms$^{-1}$)', 'interpreter','latex','fontsize',20);

subplot(4,2,4)
for z = 1:length(playa_z)
    plot(playa.oct18.Hz20.Uy_edges(:,z),playa.oct18.Hz20.Uy_pdf(:,z));
    hold on
end
set(gca,'fontsize', 15)
legend('25.5 m','19.4','10.4','5.3','2.02','0.61')
xlabel('v (ms$^{-1}$)', 'interpreter','latex','fontsize',20);

subplot(4,2,5)
for z = 1:length(sagebrush_z)
    plot(sagebrush.oct18.Hz20.Uz_edges(:,z),sagebrush.oct18.Hz20.Uz_pdf(:,z));
    hold on
end
set(gca,'fontsize', 15)
ylabel('PDF $w$', 'interpreter','latex','fontsize',20);
legend('18.6 m','10.15','5.87','2.04','0.55')
xlabel('w (ms$^{-1}$)', 'interpreter','latex','fontsize',20); 

subplot(4,2,6)
for z = 1:length(playa_z)
    plot(playa.oct18.Hz20.Uz_edges(:,z),playa.oct18.Hz20.Uz_pdf(:,z));
    hold on
end
set(gca,'fontsize', 15)
legend('25.5 m','19.4','10.4','5.3','2.02','0.61')
xlabel('w (ms$^{-1}$)', 'interpreter','latex','fontsize',20);

subplot(4,2,7)
for z = 1:length(sagebrush_z)
    plot(sagebrush.oct18.Hz1.T_edges(:,z),sagebrush.oct18.Hz1.T_pdf(:,z));
    hold on
end
set(gca,'fontsize', 15)
ylabel('PDF $T$', 'interpreter','latex','fontsize',20);
legend('18.6 m','10.15','5.87','2.04','0.55')
xlabel('T ($^\circ C$)', 'interpreter','latex','fontsize',20);

subplot(4,2,8)
for z = 1:length(playa_z)
    plot(playa.oct18.Hz1.T_edges(:,z),playa.oct18.Hz1.T_pdf(:,z));
    hold on
end
set(gca,'fontsize', 15)
legend('25.5 m','19.4','10.4','5.3','2.02','0.61')
xlabel('T ($^\circ C$)', 'interpreter','latex','fontsize',20);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                         Plot CDF 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure()
subplot(4,2,1)
for z = 1:length(sagebrush_z)
    plot(sagebrush.oct18.Hz20.Ux_edges(:,z),sagebrush.oct18.Hz20.Ux_cdf(:,z));
    hold on
end
set(gca,'fontsize', 15)
ylabel('PDF $u$', 'interpreter','latex','fontsize',20);
legend('18.6 m','10.15','5.87','2.04','0.55')
title('Sagebrush', 'interpreter','latex','fontsize',25)
xlabel('u (ms$^{-1}$)', 'interpreter','latex','fontsize',20);


subplot(4,2,2)
for z = 1:length(playa_z)
    plot(playa.oct18.Hz20.Ux_edges(:,z),playa.oct18.Hz20.Ux_cdf(:,z));
    hold on
end
set(gca,'fontsize', 15)
title('Playa', 'interpreter','latex','fontsize',25)
legend('25.5 m','19.4','10.4','5.3','2.02','0.61')
xlabel('u (ms$^{-1}$)', 'interpreter','latex','fontsize',20);

subplot(4,2,3)
for z = 1:length(sagebrush_z)
    plot(sagebrush.oct18.Hz20.Uy_edges(:,z),sagebrush.oct18.Hz20.Uy_cdf(:,z));
    hold on
end
set(gca,'fontsize', 15)
ylabel('PDF $v$', 'interpreter','latex','fontsize',20);
legend('18.6 m','10.15','5.87','2.04','0.55')
xlabel('v (ms$^{-1}$)', 'interpreter','latex','fontsize',20);

subplot(4,2,4)
for z = 1:length(playa_z)
    plot(playa.oct18.Hz20.Uy_edges(:,z),playa.oct18.Hz20.Uy_cdf(:,z));
    hold on
end
set(gca,'fontsize', 15)
legend('25.5 m','19.4','10.4','5.3','2.02','0.61')
xlabel('v (ms$^{-1}$)', 'interpreter','latex','fontsize',20);

subplot(4,2,5)
for z = 1:length(sagebrush_z)
    plot(sagebrush.oct18.Hz20.Uz_edges(:,z),sagebrush.oct18.Hz20.Uz_cdf(:,z));
    hold on
end
set(gca,'fontsize', 15)
ylabel('PDF $w$', 'interpreter','latex','fontsize',20);
legend('18.6 m','10.15','5.87','2.04','0.55')
xlabel('w (ms$^{-1}$)', 'interpreter','latex','fontsize',20); 

subplot(4,2,6)
for z = 1:length(playa_z)
    plot(playa.oct18.Hz20.Uz_edges(:,z),playa.oct18.Hz20.Uz_cdf(:,z));
    hold on
end
set(gca,'fontsize', 15)
legend('25.5 m','19.4','10.4','5.3','2.02','0.61')
xlabel('w (ms$^{-1}$)', 'interpreter','latex','fontsize',20);

subplot(4,2,7)
for z = 1:length(sagebrush_z)
    plot(sagebrush.oct18.Hz1.T_edges(:,z),sagebrush.oct18.Hz1.T_cdf(:,z));
    hold on
end
set(gca,'fontsize', 15)
ylabel('PDF $T$', 'interpreter','latex','fontsize',20);
legend('18.6 m','10.15','5.87','2.04','0.55')
xlabel('T ($^\circ C$)', 'interpreter','latex','fontsize',20);

subplot(4,2,8)
for z = 1:length(playa_z)
    plot(playa.oct18.Hz1.T_edges(:,z),playa.oct18.Hz1.T_cdf(:,z));
    hold on
end
set(gca,'fontsize', 15)
legend('25.5 m','19.4','10.4','5.3','2.02','0.61')
xlabel('T ($^\circ C$)', 'interpreter','latex','fontsize',20);
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                         Plot u and T vertical profiles
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
num_chuck = 42; 
chunk = 30*60*20;
t_start = chunk*num_chuck;
t_end = t_start+chunk;
figure()
subplot(4,2,1)
plot(mean(sagebrush.oct18.Hz20.Ux(t_start:t_end,:),1,'omitnan'),sagebrush_z,'r--*');
hold on
plot(mean(playa.oct18.Hz20.Ux(t_start:t_end,:),1),playa_z,'b--o');
set(gca,'fontsize', 15)
xlabel('$\overline{u}$ (ms$^{-1}$)', 'interpreter','latex','fontsize',20);
%legend('Sagebrush','Playa')
ylabel('z (m)', 'interpreter','latex','fontsize',20);
axis([0 2 0.55 25.5])
grid on

num_chuck = 42; 
chunk = 30*60*20;
t_start = chunk*num_chuck;
t_end = t_start+(20*60);
subplot(4,2,3)
for i = t_start:t_end
    patchline(sagebrush.oct18.Hz20.Ux_fluct(i,:),sagebrush_z,'edgecolor','r','edgealpha',0.025);
    hold on
end
hold on
for i = t_start:t_end
    patchline(playa.oct18.Hz20.Ux_fluct(i,:),playa_z,'edgecolor','b','edgealpha',0.025);
    hold on
end
set(gca,'fontsize', 15)
xlabel('u$^\prime$ (ms$^{-1}$)', 'interpreter','latex','fontsize',20);
%legend('Sagebrush','Playa')
ylabel('z (m)', 'interpreter','latex','fontsize',20);
axis([-2 2 0.55 25.5])
grid on

subplot(4,2,5)
for i = t_start:t_end
    patchline(sagebrush.oct18.Hz20.Uy_fluct(i,:),sagebrush_z,'edgecolor','r','edgealpha',0.025);
    hold on
end
hold on
for i = t_start:t_end
    patchline(playa.oct18.Hz20.Uy_fluct(i,:),playa_z,'edgecolor','b','edgealpha',0.025);
    hold on
end
set(gca,'fontsize', 15)
ylabel('z (m)', 'interpreter','latex','fontsize',20);
%legend('Sagebrush','Playa')
xlabel('v$^\prime$ (ms$^{-1}$)', 'interpreter','latex','fontsize',20);
axis([-2 2 0.55 25.5])
grid on

subplot(4,2,7)
for i = t_start:t_end
    patchline(sagebrush.oct18.Hz20.Uz_fluct(i,:),sagebrush_z,'edgecolor','r','edgealpha',0.025);
    hold on
end
hold on
for i = t_start:t_end
    patchline(playa.oct18.Hz20.Uz_fluct(i,:),playa_z,'edgecolor','b','edgealpha',0.025);
    hold on
end
set(gca,'fontsize', 15)
xlabel('w$^\prime$ (ms$^{-1}$)', 'interpreter','latex','fontsize',20);
%legend('Sagebrush','Playa')
ylabel('z (m)', 'interpreter','latex','fontsize',20);
axis([-2 2 0.55 25.5])
grid on

chunk = 30*60*1;
t_start = chunk*num_chuck;
t_end = t_start+chunk;
subplot(1,2,2)
plot(mean(sagebrush.oct18.Hz1.T(t_start:t_end,:),1,'omitnan'),sagebrush_z,'r--*');
hold on
plot(mean(playa.oct18.Hz1.T(t_start:t_end,:),1),playa_z,'b--o');
set(gca,'fontsize', 15)
xlabel('$\overline{T}$ ($^\circ$C)', 'interpreter','latex','fontsize',20);
legend('Sagebrush','Playa')
ylabel('z (m)', 'interpreter','latex','fontsize',20);
axis([15 20 0.55 25.5])
grid on

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                         Plot momentum correlations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

chunk1 = 30*60*1;
t_start1 = chunk1*num_chuck;
t_end1 = t_start1+chunk1;

sagebrush.oct18.Hz1.T_fluct = sagebrush.oct18.Hz1.T(t_start1:t_end1,:) -...
    mean(sagebrush.oct18.Hz1.T(t_start1:t_end1,:),1);

playa.oct18.Hz1.T_fluct = playa.oct18.Hz1.T(t_start1:t_end1,:) -...
    mean(playa.oct18.Hz1.T(t_start1:t_end1,:),1);
h = 3; 
figure()
subplot(3,2,1)
plot(sagebrush.oct18.Hz20.Uy_fluct(t_start:t_end,h),sagebrush.oct18.Hz20.Ux_fluct(t_start:t_end,h),'r*');
hold on
plot(playa.oct18.Hz20.Uy_fluct(t_start:t_end,h),playa.oct18.Hz20.Ux_fluct(t_start:t_end,h),'bo');
set(gca,'fontsize', 15)
xlabel('$v^\prime$ (ms$^{-1}$)', 'interpreter','latex','fontsize',20);
legend('Sagebrush','Playa')
ylabel('$u^\prime$ (ms$^{-1}$)', 'interpreter','latex','fontsize',20);
axis([-2 2 -3 3])
plot([-5 5],[0 0],'k','LineWidth', 2)
plot([0 0],[-5 5],'k','LineWidth', 2)
grid on

subplot(3,2,2)
plot(sagebrush.oct18.Hz20.Uz_fluct(t_start:t_end,h),sagebrush.oct18.Hz20.Ux_fluct(t_start:t_end,h),'r*');
hold on
plot(playa.oct18.Hz20.Uz_fluct(t_start:t_end,h),playa.oct18.Hz20.Ux_fluct(t_start:t_end,h),'bo');
set(gca,'fontsize', 15)
xlabel('$w^\prime$ (ms$^{-1}$)', 'interpreter','latex','fontsize',20);
legend('Sagebrush','Playa')
ylabel('$u^\prime$ (ms$^{-1}$)', 'interpreter','latex','fontsize',20);
axis([-2 2 -3 3])
plot([-5 5],[0 0],'k','LineWidth', 2)
plot([0 0],[-5 5],'k','LineWidth', 2)
grid on

subplot(3,2,3)
plot(sagebrush.oct18.Hz20.Uz_fluct(t_start:t_end,h),sagebrush.oct18.Hz20.Uy_fluct(t_start:t_end,h),'r*');
hold on
plot(playa.oct18.Hz20.Uz_fluct(t_start:t_end,h),playa.oct18.Hz20.Uy_fluct(t_start:t_end,h),'bo');
set(gca,'fontsize', 15)
xlabel('$w^\prime$ (ms$^{-1}$)', 'interpreter','latex','fontsize',20);
legend('Sagebrush','Playa')
ylabel('$v^\prime$ (ms$^{-1}$)', 'interpreter','latex','fontsize',20);
axis([-2 2 -3 3])
plot([-5 5],[0 0],'k','LineWidth', 2)
plot([0 0],[-5 5],'k','LineWidth', 2)
grid on

subplot(3,2,4)
plot(sagebrush.oct18.Hz1.T_fluct(:,h),sagebrush.oct18.Hz20.Uz_fluct(t_start:20:t_end,h),'r*');
hold on
plot(playa.oct18.Hz1.T_fluct(:,h),playa.oct18.Hz20.Uz_fluct(t_start:20:t_end,h),'bo');
set(gca,'fontsize', 15)
xlabel('$T^\prime$ ($^\circ$C)', 'interpreter','latex','fontsize',20);
legend('Sagebrush','Playa')
ylabel('$w^\prime$ (ms$^{-1}$)', 'interpreter','latex','fontsize',20);
axis([-2 2 -3 3])
plot([-5 5],[0 0],'k','LineWidth', 2)
plot([0 0],[-5 5],'k','LineWidth', 2)
grid on

subplot(3,2,5)
plot(sagebrush.oct18.Hz1.T_fluct(:,h),sagebrush.oct18.Hz20.Ux_fluct(t_start:20:t_end,h),'r*');
hold on
plot(playa.oct18.Hz1.T_fluct(:,h),playa.oct18.Hz20.Ux_fluct(t_start:20:t_end,h),'bo');
set(gca,'fontsize', 15)
xlabel('$T^\prime$ ($^\circ$C)', 'interpreter','latex','fontsize',20);
legend('Sagebrush','Playa')
ylabel('$u^\prime$ (ms$^{-1}$)', 'interpreter','latex','fontsize',20);
axis([-2 2 -3 3])
plot([-5 5],[0 0],'k','LineWidth', 2)
plot([0 0],[-5 5],'k','LineWidth', 2)
grid on

subplot(3,2,6)
plot(sagebrush.oct18.Hz1.T_fluct(:,h),sagebrush.oct18.Hz20.Uy_fluct(t_start:20:t_end,h),'r*');
hold on
plot(playa.oct18.Hz1.T_fluct(:,h),playa.oct18.Hz20.Uy_fluct(t_start:20:t_end,h),'bo');
set(gca,'fontsize', 15)
xlabel('$T^\prime$ ($^\circ$C)', 'interpreter','latex','fontsize',20);
legend('Sagebrush','Playa')
ylabel('$v^\prime$ (ms$^{-1}$)', 'interpreter','latex','fontsize',20);
axis([-2 2 -3 3])
plot([-5 5],[0 0],'k','LineWidth', 2)
plot([0 0],[-5 5],'k','LineWidth', 2)
grid on


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           Calculate 3rd and 4th moments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if recalc

%Skewness and Kurtosis corr on u,v,and w
chunk = 30*60*20;
t_start = chunk*num_chuck;
t_end = t_start+chunk;

for z = 1:length(sagebrush_z)
    sagebrush.oct18.Hz20.Ux_skew(:,z) = skewness(sagebrush.oct18.Hz20.Ux(t_start:t_end,z));
    sagebrush.oct18.Hz20.Uy_skew(:,z) = skewness(sagebrush.oct18.Hz20.Uy(t_start:t_end,z));
    sagebrush.oct18.Hz20.Uz_skew(:,z) = skewness(sagebrush.oct18.Hz20.Uz(t_start:t_end,z));
    sagebrush.oct18.Hz20.Ux_kurt(:,z) = kurtosis(sagebrush.oct18.Hz20.Ux(t_start:t_end,z));
    sagebrush.oct18.Hz20.Uy_kurt(:,z) = kurtosis(sagebrush.oct18.Hz20.Uy(t_start:t_end,z));
    sagebrush.oct18.Hz20.Uz_kurt(:,z) = kurtosis(sagebrush.oct18.Hz20.Uz(t_start:t_end,z));
end

for z = 1:length(playa_z)
    playa.oct18.Hz20.Ux_skew(:,z) = skewness(playa.oct18.Hz20.Ux(t_start:t_end,z));
    playa.oct18.Hz20.Uy_skew(:,z) = skewness(playa.oct18.Hz20.Uy(t_start:t_end,z));
    playa.oct18.Hz20.Uz_skew(:,z) = skewness(playa.oct18.Hz20.Uz(t_start:t_end,z));
    playa.oct18.Hz20.Ux_kurt(:,z) = kurtosis(playa.oct18.Hz20.Ux(t_start:t_end,z));
    playa.oct18.Hz20.Uy_kurt(:,z) = kurtosis(playa.oct18.Hz20.Uy(t_start:t_end,z));
    playa.oct18.Hz20.Uz_kurt(:,z) = kurtosis(playa.oct18.Hz20.Uz(t_start:t_end,z));
end

%Skewness and Kurtosis corr on T
chunk = 30*60*1;
t_start = chunk*num_chuck;
t_end = t_start+chunk;
for z = 1:length(sagebrush_z)
    sagebrush.oct18.Hz1.T_skew(:,z) = skewness(sagebrush.oct18.Hz1.T(t_start:t_end,z));
    sagebrush.oct18.Hz1.T_kurt(:,z) = kurtosis(sagebrush.oct18.Hz1.T(t_start:t_end,z));
end
for z = 1:length(playa_z)
    playa.oct18.Hz1.T_skew(:,z) = skewness(playa.oct18.Hz1.T(t_start:t_end,z));
    playa.oct18.Hz1.T_kurt(:,z) = kurtosis(playa.oct18.Hz1.T(t_start:t_end,z));
end
%%
% 3. Autocorrelation: Calculate the autocorrelation of at least one 30-minute period. What does this
% indicate?
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           Calculate Autocorrelations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Auto corr on u v and w
chunk = 30*60*20;
t_start = chunk*num_chuck;
t_end = t_start+chunk;
%tic
for z = 1:length(sagebrush_z)
    sagebrush.oct18.Hz20.Ux_autocorr(:,z) = AutoCorrelation(sagebrush.oct18.Hz20.Ux(t_start:t_end,z));
    sagebrush.oct18.Hz20.Uy_autocorr(:,z) = AutoCorrelation(sagebrush.oct18.Hz20.Uy(t_start:t_end,z));
    sagebrush.oct18.Hz20.Uz_autocorr(:,z) = AutoCorrelation(sagebrush.oct18.Hz20.Uz(t_start:t_end,z));
end


for z = 1:length(playa_z)
    playa.oct18.Hz20.Ux_autocorr(:,z) = AutoCorrelation(playa.oct18.Hz20.Ux(t_start:t_end,z));
    playa.oct18.Hz20.Uy_autocorr(:,z) = AutoCorrelation(playa.oct18.Hz20.Uy(t_start:t_end,z));
    playa.oct18.Hz20.Uz_autocorr(:,z) = AutoCorrelation(playa.oct18.Hz20.Uz(t_start:t_end,z));
end

%Auto corr on T
chunk = 30*60*1;
t_start = chunk*num_chuck;
t_end = t_start+chunk;
for z = 1:length(sagebrush_z)
    sagebrush.oct18.Hz1.T_autocorr(:,z) = AutoCorrelation(sagebrush.oct18.Hz1.T(t_start:t_end,z));
end
for z = 1:length(playa_z)
    playa.oct18.Hz1.T_autocorr(:,z) = AutoCorrelation(playa.oct18.Hz1.T(t_start:t_end,z));
end
%toc

end %if recalc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           Plot Autocorrelations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x_axis1hz = linspace(0,15,length(sagebrush.oct18.Hz1.T_autocorr(:,1)));
x_axis20hz = linspace(0,15,length(sagebrush.oct18.Hz20.Ux_autocorr(:,1)));

figure()
subplot(4,2,1)
plot(x_axis20hz,sagebrush.oct18.Hz20.Ux_autocorr);
set(gca,'fontsize', 15)
ylabel('$\rho_u$', 'interpreter','latex','fontsize',20);
legend('18.6 m','10.15','5.87','2.04','0.55')
title('Sagebrush', 'interpreter','latex','fontsize',25)
axis([0 15 -1 1]); 

subplot(4,2,2)
plot(x_axis20hz,playa.oct18.Hz20.Ux_autocorr);
set(gca,'fontsize', 15)
title('Playa', 'interpreter','latex','fontsize',25)
legend('25.5 m','19.4','10.4','5.3','2.02','0.61')
axis([0 15 -1 1]); 

subplot(4,2,3)
plot(x_axis20hz,sagebrush.oct18.Hz20.Uy_autocorr);
set(gca,'fontsize', 15)
ylabel('$\rho_v$', 'interpreter','latex','fontsize',20);
legend('18.6 m','10.15','5.87','2.04','0.55')
axis([0 15 -1 1]); 

subplot(4,2,4)
plot(x_axis20hz,playa.oct18.Hz20.Uy_autocorr);
set(gca,'fontsize', 15)
legend('25.5 m','19.4','10.4','5.3','2.02','0.61')
axis([0 15 -1 1]); 

subplot(4,2,5)
plot(x_axis20hz,sagebrush.oct18.Hz20.Uz_autocorr);
set(gca,'fontsize', 15)
ylabel('$\rho_w$', 'interpreter','latex','fontsize',20);
legend('18.6 m','10.15','5.87','2.04','0.55')
axis([0 15 -1 1]); 

subplot(4,2,6)
plot(x_axis20hz,playa.oct18.Hz20.Uz_autocorr);
set(gca,'fontsize', 15)
legend('25.5 m','19.4','10.4','5.3','2.02','0.61')
axis([0 15 -1 1]); 

subplot(4,2,7)
plot(x_axis1hz,sagebrush.oct18.Hz1.T_autocorr);
set(gca,'fontsize', 15)
ylabel('$\rho_T$', 'interpreter','latex','fontsize',20);
legend('18.6 T','10.15','5.87','2.04','0.55')
xlabel('time (min)', 'interpreter','latex','fontsize',20);
axis([0 15 -1 1]); 

subplot(4,2,8)
plot(x_axis1hz,playa.oct18.Hz1.T_autocorr);
set(gca,'fontsize', 15)
legend('25.5 m','19.4','10.4','5.3','2.02','0.61')
xlabel('time (min)', 'interpreter','latex','fontsize',20);
axis([0 15 -1 1]); 



%% 4. 
% Dissipation: Using Taylor's frozen turbulence hypothesis, calculate the dissipation rate of turbulent
% kinetic energy for several 30-minute periods. Calculate the Kolmogorov length scale.
chunk = 30*60*20;
t_start = chunk*num_chuck;
t_end = t_start+chunk;

IDX = 3;
diss_rate_playa     = zeros(IDX);
kol_length_playa    = zeros(IDX);
diss_rate_sagebrush = zeros(IDX);
kol_length_sagebrush = zeros(IDX);
for i=1:IDX
    Vel_mag_playa = sqrt(playa.oct18.Hz20.Ux(t_start:t_end,5).^2 ...
        + playa.oct18.Hz20.Uy(t_start:t_end,5).^2);
    [diss_rate_playa(i) kol_length_playa(i)] = taylors(Vel_mag_playa);
    
    Vel_mag_sagebrush = sqrt(sagebrush.oct18.Hz20.Ux(t_start:t_end,5).^2 ...
        + sagebrush.oct18.Hz20.Uy(t_start:t_end,5).^2);
    [diss_rate_sagebrush(i) kol_length_sagebrush(i)] = taylors(Vel_mag_sagebrush);
    t_start = t_end;
    t_end =t_start+chunk;
end

%% 5.
%  Turbulence Spectra: For at least one 30-minute averaging period, calculate the following turbulent
% energy spectra: Suu, Svv, Sww, ST T , as well as the following cospectra: Suw, Svw, Swt. What do
% these spectra indicate about the analyzed boundary layer? Is there an inertial subrange?

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           Calculate Turbulence Spectra
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Turbulence Spectra Suu, Svv, Sww, Stt
% Cospectra Suw, Svw, Swt
% Calculated using 2 meter height sensor
chunk = 30*60*20;
t_start = chunk*num_chuck;
t_end = t_start+chunk;

% Spectra u, v, w
[sagebrush_Ux_spectra_ww, sagebrush_Ux_spectra_ff]...
    = spectra(sagebrush.oct18.Hz20.Ux(t_start:t_end,5));
[sagebrush_Uy_spectra_ww, sagebrush_Uy_spectra_ff]...
    = spectra(sagebrush.oct18.Hz20.Uy(t_start:t_end,5));
[sagebrush_Uz_spectra_ww, sagebrush_Uz_spectra_ff]...
    = spectra(sagebrush.oct18.Hz20.Uz(t_start:t_end,5));

[playa_Ux_spectra_ww, playa_Ux_spectra_ff]...
    = spectra(playa.oct18.Hz20.Ux(t_start:t_end,5));
[playa_Uy_spectra_ww, playa_Uy_spectra_ff]...
    = spectra(playa.oct18.Hz20.Uy(t_start:t_end,5));
[playa_Uz_spectra_ww, playa_Uz_spectra_ff]...
    = spectra(playa.oct18.Hz20.Uz(t_start:t_end,5));

% Spectra on T
chunk = 30*60*1;
t_start = chunk*num_chuck;
t_end = t_start+chunk;

[playa_T_spectra_ww, playa_T_spectra_ff]...
    = spectra(playa.oct18.Hz1.T(t_start:t_end,5));
[sagebrush_T_spectra_ww, sagebrush_T_spectra_ff]...
    = spectra(sagebrush.oct18.Hz1.T(t_start:t_end,5));

% Cospectra uw, vw

chunk = 30*60*20;
t_start = chunk*num_chuck;
t_end = t_start+chunk;

[sagebrush_UxUz_spectra_ww, sagebrush_UxUz_spectra_ff]...
    = cospectra(sagebrush.oct18.Hz20.Ux(t_start:t_end,5),...
    sagebrush.oct18.Hz20.Uz(t_start:T_end,5));
[sagebrush_UyUz_spectra_ww, sagebrush_Uy_spectra_ff]...
    = spectra(sagebrush.oct18.Hz20.Uy(t_start:t_end,5),...
    sagebrush.oct18.Hz20.Uz(t_start:T_end,5));

[playa_UxUz_spectra_ww, playa_UxUz_spectra_ff]...
    = spectra(playa.oct18.Hz20.Ux(t_start:t_end,5),...
    playa.oct18.Hz20.Uz(t_start:T_end,5));
[playa_UyUz_spectra_ww, playa_UyUz_spectra_ff]...
    = spectra(playa.oct18.Hz20.Uy(t_start:t_end,5),...
    playa.oct18.Hz20.Uz(t_start:T_end,5));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           Plot Spectra
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure();
loglog(playa_Ux_spectra_ww,playa_Ux_spectra_ff,'r--')
hold on
loglog(sagebrush_Ux_spectra_ww,sagebrush_Ux_spectra_ff,'b-')
ww_slope = 0.05:0.1:3;
loglog(ww_slope, 10^-6.5 * ww_slope.^(-5/3),'k-')
set(gca,'fontsize', 15)
ylabel('$S_{uu}$', 'interpreter','latex','fontsize',20);
title('$S_{uu}$', 'interpreter','latex','fontsize',25)
xlabel('f (Hz)', 'interpreter','latex','fontsize',25)
legend('Playa','Sagebrush','(-5/3) slope')

figure();
loglog(playa_Uy_spectra_ww,playa_Uy_spectra_ff,'r--')
hold on
loglog(sagebrush_Uy_spectra_ww,sagebrush_Uy_spectra_ff,'b-')
ww_slope = 0.05:0.1:3;
loglog(ww_slope, 10^-6.5 * ww_slope.^(-5/3),'k-')
set(gca,'fontsize', 15)
ylabel('$S_{vv}$', 'interpreter','latex','fontsize',20);
title('$S_{vv}$', 'interpreter','latex','fontsize',25)
xlabel('f (Hz)', 'interpreter','latex','fontsize',25)
legend('Playa','Sagebrush','(-5/3) slope')

figure();
loglog(playa_Uz_spectra_ww,playa_Uz_spectra_ff,'r--')
hold on
loglog(sagebrush_Uz_spectra_ww,sagebrush_Uz_spectra_ff,'b-')
ww_slope = 0.05:0.1:3;
loglog(ww_slope, 10^-6.5 * ww_slope.^(-5/3),'k-')
set(gca,'fontsize', 15)
ylabel('$S_{ww}$', 'interpreter','latex','fontsize',20);
title('$S_{ww}$', 'interpreter','latex','fontsize',25)
xlabel('f (Hz)', 'interpreter','latex','fontsize',25)
legend('Playa','Sagebrush','(-5/3) slope')

figure();
loglog(playa_T_spectra_ww,playa_T_spectra_ff,'r--')
hold on
loglog(sagebrush_T_spectra_ww,sagebrush_T_spectra_ff,'b-')
ww_slope = 0.05:0.1:3;
loglog(ww_slope, 10^-6.5 * ww_slope.^(-5/3),'k-')
set(gca,'fontsize', 15)
ylabel('$S_{TT}$', 'interpreter','latex','fontsize',20);
title('$S_{TT}$', 'interpreter','latex','fontsize',25)
xlabel('f (Hz)', 'interpreter','latex','fontsize',25)
legend('Playa','Sagebrush','(-5/3) slope')


figure();
loglog(playa_UxUz_spectra_ww,playa_UxUz_spectra_ff,'r--')
hold on
loglog(sagebrush_UxUz_spectra_ww,sagebrush_UxUz_spectra_ff,'b-')
set(gca,'fontsize', 15)
ylabel('$S_{uw}$', 'interpreter','latex','fontsize',20);
title('$S_{uw}$', 'interpreter','latex','fontsize',25)
xlabel('f (Hz)', 'interpreter','latex','fontsize',25)
legend('Playa','Sagebrush')

figure();
loglog(playa_UyUz_spectra_ww,playa_UyUz_spectra_ff,'r--')
hold on
loglog(sagebrush_UyUz_spectra_ww,sagebrush_UyUz_spectra_ff,'b-')
set(gca,'fontsize', 15)
ylabel('$S_{uw}$', 'interpreter','latex','fontsize',20);
title('$S_{uw}$', 'interpreter','latex','fontsize',25)
xlabel('f (Hz)', 'interpreter','latex','fontsize',25)
legend('Playa','Sagebrush')
