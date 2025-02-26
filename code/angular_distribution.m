load 'coords.mat'

centerM=zeros(3,1);
centerM(1,1)=mean(coords(1,:));
centerM(2,1)=mean(coords(2,:));
centerM(3,1)=mean(coords(3,:));
  
x_axis= coords(:,1)-centerM(:,1);
x_axis=x_axis/norm(x_axis);

y_axis=zeros(3,1);
y_axis(1,1)=x_axis(2,1)/sqrt(x_axis(1,1)^2+x_axis(2,1)^2);
y_axis(2,1)= -x_axis(1,1)*y_axis(1,1)/(x_axis(2,1));
y_axis(3,1)=0;

z_axis=cross(x_axis,y_axis);
z_axis= z_axis/norm(z_axis);

initial_cells= coords(:,1:numin);

radial=zeros(3,size(initial_cells,2));

for i=1:size(initial_cells,2)
    
    v1=initial_cells(:,i)-centerM(:,1);
    radial(1,i)=norm(v1);
    radial(2,i)=acos(dot(z_axis,v1)/radial(1,i));
    
    v2=acos(dot(x_axis,v1)/(radial(1,i)*sin(radial(2,i))));
    v3=asin(dot(y_axis,v1)/(radial(1,i)*sin(radial(2,i))));
    if(sin(v3)>=0)
      radial(3,i)=acos(dot(x_axis,v1)/(radial(1,i)*sin(radial(2,i))));
    end
%    if(sin(v2)>0)
%       radial(3,i)=acos(dot(x_axis,v1)/(radial(1,i)*sin(radial(2,i))));
%    end
    if( sin(v3)<0)
      radial(3,i)=2*pi-acos(dot(x_axis,v1)/(radial(1,i)*sin(radial(2,i))));
    end
    
    
%     if(cos(v2)>0 && sin(v2)<0)
%       radial(3,i)=2*pi-acos(dot(x_axis,v1)/(radial(1,i)*sin(radial(2,i))));
%     end
    
    
    
end
  figure(1)
  histogram(radial(1,:),20)
  h_gca=gca;
%  h=h_gca.Children;
%  h.FaceColor=[.98 .98 .98];
%  h.EdgeColor=[.94 .94 .94];


  figure(2)
binWidth = 0.1;
 binCtrs = 0:binWidth:pi;
 hist(radial(2,:),binCtrs);
% %xlim([0 pi]);
% hist(radial(2,:),binCtrs);
% xlim([0 pi]);
% 
% counts=hist(radial(2,:),binCtrs);
% N = length(radial(2,:));
% x1=binCtrs;
% x1=sin(x1);
% for i=1:size(counts,2)
%     prob(i)=counts(i)/(N*binWidth*sin(binCtrs(i)));
% end
% figure(2)
% 
% 
%  plot(binCtrs(2:30),prob(2:30),'b-o');
% % xlim([0 pi]);
% % ylim([0 1]);
% 
% % 
  figure(3)
  binWidth = 0.3;
 binCtrs = 0:binWidth:2*pi;
 
hist(radial(3,:),binCtrs);
% % % hist(radial(3,:),binCtrs)
% %   N = length(radial(3,:));
% %   prob=counts/(N*binWidth);
% %   plot(binCtrs,prob,'r-o')
% %scatter3(coords(1,:),coords(2,:),coords(3,:))
% 
% 
% 
