function [ww_new ff_new] = cospectra(uu2,vv2)

tbin = 0.1;         % Time bins (minutes)
tbin = tbin * 60;   % seconds

uu2 = uu2 - mean(uu2); % perturbation u'
vv2 = vv2 - mean(vv2); % perturbation v'

u2=uu2;
v2=vv2;

F=20;       % Sampling Frequency
Fn = F/2;   % Nyquist Frequency
tt2=1:length(uu2);
tt2=tt2/F;
T=max(tt2);
nbins=T/tbin;

ti=T/nbins;

k=1;
while ti<=max(tt2)-T/nbins
      jj=1;
for i=1:length(tt2)
    if tt2(i)>=ti-T/nbins && tt2(i)<=ti+T/nbins
        u2(jj)=uu2(i);
        v2(jj)=vv2(i);
        jj=jj+1;
    end
end
    k=k+1;
    ti=ti+T/nbins;
    lll=length(u2);
    N=length(u2);
    IDX = floor(N/2);
    p = (real(fft(u2)) .* real(fft(v2))  + imag(fft(v2)) .* imag(fft(v2)))/((N/2)^2);
    p = p(1:IDX);
    L=length(p);
    ww(k,1:L) = Fn * (0:IDX-1)/ (IDX-1);
    ffiuu2(k,1:L)=p;
end
s=size(ffiuu2);
for m=1:s(2)
    ff_average(m)=mean(ffiuu2(:,m));
    ww_average(m)=mean(ww(:,m));
end
%Post process averaging to smooth the data
%smooth out ff_average
binL = 1; %#pts
count = 1;
binC = 0;
newIndex = 1;
while count < length(ff_average)
    if ww_average(count) <= 3*10^-3
        ff_sum = ff_average(count);
        ww_sum = ww_average(count);
        binL = 1;
        binC = binC+1;
    elseif ww_average(count)>3*10^-3 && ww_average(count) <= 2*10^-2
        ff_sum = ff_sum + ff_average(count);
        ww_sum = ww_sum + ww_average(count);
        binL = 10;
        binC = binC+1;          
    elseif ww_average(count)>2*10^-2 && ww_average(count) <= 2*10^-1
        ff_sum = ff_sum + ff_average(count);
        ww_sum = ww_sum + ww_average(count);
        binL = 20;
        binC = binC+1;
    elseif ww_average(count)>2*10^-1 && ww_average(count) <= 2*10^0
        ff_sum = ff_sum + ff_average(count);
        ww_sum = ww_sum + ww_average(count);
        binL = 50;
        binC = binC+1;
    elseif ww_average(count)>2*10^0 && ww_average(count) <= 10^2
        ff_sum = ff_sum + ff_average(count);
        ww_sum = ww_sum + ww_average(count);
        binL = 100;
        binC = binC+1;
    end
    
    if binC == binL
        ff_new(newIndex) = ff_sum/binL;
        ww_new(newIndex) = ww_sum/binL;
        binC = 0;
        newIndex = newIndex + 1;
        ff_sum = 0;
        ww_sum = 0;
    end
    count = count + 1;
end


















