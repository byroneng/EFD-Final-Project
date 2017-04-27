function [diss_rate kol_length] = taylors(uu2)
F = 20;     % Sampling Frequency
Fn = F/2;   % Nyquist Frequency

tt2=1:length(uu2);
tt2=tt2/F;
upert=uu2-mean(uu2);

% Dissipation rate and Kolmogorov Scale
% x=aveu1*t; -> dx  = aveu1 * dt % Taylor Frozen Hypothesis

vis = 1.628*10^-5/0.7364;

dt = tt2(2) - tt2(1);
dupert_dt = (upert(2:length(upert)) - upert(1:length(upert)-1))/(mean(uu2)*dt);
diss_rate = 15*vis * mean(dupert_dt.^2); % Dissipation Rate
kol_length = (vis^3 / diss_rate)^(1/4) * 10^3; % Kolmogorov Length Scale (mm) * 10^3 % Kolmogorov Length Scale (mm) 
end


















