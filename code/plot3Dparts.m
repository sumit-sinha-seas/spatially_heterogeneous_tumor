function plot3Dparts(coords,nPart,rad)
%clear s1;
%clear sph;
[Xs,Ys,Zs] = sphere;

if (size(coords,1)~=3)
return
end


%figure
%plot3(coords(1,:),coords(2,:),coords(3,:))
%colorVec = hsv(100);

for step=1:nPart
sph=[coords(1,step), coords(2,step), coords(3,step),rad(step,1)];
s1=surf(Xs*sph(1,4)+sph(1,1),Ys*sph(1,4)+sph(1,2), ...
        Zs*sph(1,4)+sph(1,3));
hold on;
%axis( [0 5 0 5 0 5]);
end



end
