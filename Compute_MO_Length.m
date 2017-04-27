function [ L ] = Compute_MO_Length( rho, tau, Q, Theta, u_star)
%Compute the MO length
%If you do not have friction velocity, set u_star = 0

k = 0.4; %Von Karman
g = 9.8; %gravitational constant on earth

if u_star == 0
    u_star = sqrt(tau./rho); %shear velocity
end

L = - (u_star.^3.*Theta)./(k.*g.*Q); %MO length 
end

