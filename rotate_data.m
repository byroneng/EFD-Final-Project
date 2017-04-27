function [ u_rot, v_rot, w_rot, u_flut, v_flut, w_flut, out ] = rotate_data(f,period,Npoints, u, v, w, num_heights, detrend_chunk )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here





%Initializing the rotated variables for speed-up reasons.
mu = zeros(Npoints,num_heights);
mv = zeros(Npoints,num_heights);
mw = zeros(Npoints,num_heights);
u_rot = zeros(Npoints,num_heights);
v_rot = zeros(Npoints,num_heights);
w_rot = zeros(Npoints,num_heights);
u_flut = zeros(Npoints,num_heights);
v_flut = zeros(Npoints,num_heights);
w_flut = zeros(Npoints,num_heights);

% Beginning of the rotation process:

for z = 1:num_heights
    %validate data
    for i = 1:Npoints
        if u(i,z) < -100 || u(i,z) > 100
            u(i,z) = NaN;
        end
         if v(i,z) < -100 || v(i,z) > 100
            v(i,z) = NaN;
         end
         if w(i,z) < -100|| w(i,z) > 100
            w(i,z) = NaN;
        end
        
    end
    cnt = 1;
    for c = 1:floor((Npoints/detrend_chunk))
        % ... calculate the averages
        mu(cnt:cnt+detrend_chunk,z) = mean(u(cnt:cnt+detrend_chunk,z),'omitnan');
        mv(cnt:cnt+detrend_chunk,z) = mean(v(cnt:cnt+detrend_chunk,z),'omitnan');
        mw(cnt:cnt+detrend_chunk,z) = mean(w(cnt:cnt+detrend_chunk,z),'omitnan');
        
        % ... calculate the wind speed within the window
        out.meanU(cnt:cnt+detrend_chunk,z) = sqrt(mu(cnt:cnt+detrend_chunk,z).^2 + mv(cnt:cnt+detrend_chunk,z).^2 + mw(cnt:cnt+detrend_chunk,z).^2);
   
        % ... calculate wind direction by assuming that the original
       
        out.dir(cnt:cnt+detrend_chunk,z) = mod(270 - (atan2(mv(cnt:cnt+detrend_chunk,z),mu(cnt:cnt+detrend_chunk,z)) * 180/pi),360);
        
        % ... calculate the first rotation angle for DR - rotate into mean wind
        alpha = atan2(mv(cnt:cnt+detrend_chunk,z),mu(cnt:cnt+detrend_chunk,z)).*180/pi;
        
        % ... rotate the data
        u1(cnt:cnt+detrend_chunk,z) = u(cnt:cnt+detrend_chunk,z).*cosd(alpha)+v(cnt:cnt+detrend_chunk,z).*sind(alpha);
        v1(cnt:cnt+detrend_chunk,z) = -u(cnt:cnt+detrend_chunk,z).*sind(alpha)+v(cnt:cnt+detrend_chunk,z).*cosd(alpha);
        w1(cnt:cnt+detrend_chunk,z) = w(cnt:cnt+detrend_chunk,z);
        
        % ... caluculate the second rotation angle - rotate to make meanW = 0
        beta(cnt:cnt+detrend_chunk,z) = atan2(nanmean(w1(cnt:cnt+detrend_chunk,z)),nanmean(u1(cnt:cnt+detrend_chunk,z))).*180/pi;
        
        % ... components in the rotated coordinate system
        u_rot(cnt:cnt+detrend_chunk,z) = u1(cnt:cnt+detrend_chunk,z).*cosd(beta(cnt:cnt+detrend_chunk,z))+w1(cnt:cnt+detrend_chunk,z).*sind(beta(cnt:cnt+detrend_chunk,z));
        v_rot(cnt:cnt+detrend_chunk,z) = v1(cnt:cnt+detrend_chunk,z);
        w_rot(cnt:cnt+detrend_chunk,z) = -u1(cnt:cnt+detrend_chunk,z).*sind(beta(cnt:cnt+detrend_chunk,z))+w1(cnt:cnt+detrend_chunk,z).*cosd(beta(cnt:cnt+detrend_chunk,z));
        % End of rotation
        %------------------------------------------------------
        % ... calculate the deviations from the mean by detrending the signal 
        u_flut(cnt:cnt+detrend_chunk,z)= detrend_nan(u_rot(cnt:cnt+detrend_chunk,z));
        v_flut(cnt:cnt+detrend_chunk,z)= detrend_nan(v_rot(cnt:cnt+detrend_chunk,z));
        w_flut(cnt:cnt+detrend_chunk,z) = detrend_nan(w_rot(cnt:cnt+detrend_chunk,z));
        cnt = cnt+detrend_chunk;
    end 
    %------------------------------------------------------

end

end



