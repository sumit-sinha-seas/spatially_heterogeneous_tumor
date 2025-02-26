function plotscatter(coords,rad,nPart)

for step=1:nPart

scatter(coords(1,step),coords(2,step));
hold on;

end

end
