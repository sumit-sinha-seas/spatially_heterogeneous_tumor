function vec = distPBC3D(vec,L)

% Calculate the half box size in each direction
hL = L/2.0;

% Distance vector should be in the range -hLx -> hLx and -hLy -> hLy
% Therefore, we need to apply the following changes if it's not in this range:

for dim=1:3
if (vec(dim) > hL)
vec(dim) = vec(dim)-L;
elseif (vec(dim) < -hL)
vec(dim) = vec(dim)+L;
end
end


end
