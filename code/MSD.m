function[msd,s1] = MSD(position,dt,nPart)

for part=1:nPart
for n = 0:1:sqrt(size(position,1))
msd(n+1,part) = mean((position(n+1:end,part)-position(1:end-n,part)).^2);
end
s1(part) = dt*[0:1:(size(msd,1)-1)];
end


figure
hold on;
for part=1:nPart
plot(s1(part),msd(:,part))
end

xlabel('s[s]')
ylabel('MSD [m^2]')

end
