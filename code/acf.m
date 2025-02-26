function [r,s]=acf(x1,dt)

avgv = (x1(2:end)-x1(1:end-1))/dt;
r = xcorr(avgv, ceil(sqrt(length(x1))),'unbiased');
s = dt*[0:1:length(r)-1];

figure
plot(s,r)
xlabel('s[s]')
ylabel('R [m^2]')

end
