index=0;

dL = 4;
L=200;
x=L/dL;
count(x,1) = zeros; 
count1(x,1) = zeros;

ctrMx = 0;
ctrMy = 0;
ctrMz= 0;

for i=1:nPart
            ctrMx = ctrMx + coordsx(i,1)/nPart;
            ctrMy = ctrMy + coordsy(i,1)/nPart;
            ctrMz= ctrMz + coordsz(i,1)/nPart;
end



 
for k=1:part
 
    rctr = sqrt((coordsxl(k,1)-ctrMx)^2 + ...
    (coordsyl(k,1)-ctrMy)^2 + (coordszl(k,1)-ctrMz)^2);

  for dist = 1:dL:L;

    if dist == 1
    index = dist;
    else
    index = index+1;
    end
    
    if (rctr >= dist) && (rctr < (dist+dL))
        count(index,1) = count(index,1)+(1/(4*pi*dL*part*dist^2));
    end
  end
 end


% for index = 1:1:max(index);
% count(index,1)=count(index,1)/(4*pi*dL*part);
% end
figure
distvel = 1:dL:L;
plot(distvel,count,'b','LineWidth',2);
hold on

for k=1:part1
 
    rctr = sqrt((coordsxh(k,1)-ctrMx)^2 + ...
    (coordsyh(k,1)-ctrMy)^2 + (coordszh(k,1)-ctrMz)^2);

  for dist = 1:dL:L;

    if dist == 1
    index = dist;
    else
    index = index+1;
    end
    
    if (rctr >= dist) && (rctr < (dist+dL))
        count1(index,1) = count1(index,1)+(1/(4*pi*dL*part1*dist^2));
    end
  end
 end


% for index = 1:1:max(index);
% count1(index,1)=count1(index,1)/(4*pi*dL*part1);
% end


plot(distvel,count1,'r','LineWidth',2);
xlabel('r-r_{cm}(\mum)')
ylabel('\rho')
legend('Lo','Hi')
title('f=')