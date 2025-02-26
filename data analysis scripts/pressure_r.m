clear;
d1=dir;

%p_c=5*10^(-4);

no_of_folders=25;
total_no_of_particles=0;
data_length=20;
fraction_folders=zeros(data_length,no_of_folders);  
pressure_profile=zeros(1,2);
begin_count=8;
count_er=0;

for looper=begin_count:begin_count+(no_of_folders-1)
    looper
 cd (d1(looper).name);

%if exist('lifetime2.mat', 'file')
   
        load('lifetime1.txt');
        count_er=count_er+1;


        %fraction=zeros(data_length,1);
        frame=1200;
          %  for frame=start:data_length
                    time_inquired=frame*500;
                    data_required=lifetime1(find(lifetime1(:,6)==time_inquired),:);
                    nPart=size(data_required,1);
% coords=zeros(2,nPart);
                    coords=data_required(:,1:3)';
                    rad=data_required(:,11);
                    poisson=data_required(:,10);
                    modulus=data_required(:,9);
                    receptor=data_required(:,7);
                    ligand=data_required(:,8);

                    etab=0.0001;

% Initialize all forces to 0
                    forces = zeros(size(coords));
                     gamma3 = zeros(size(rad));
                     pressure = zeros(size(rad));
%bond_orientation=zeros(size(rad));
%bond_orientation_mag=zeros(size(rad));
                    coordination_number=zeros(size(rad));

% Get the number of particles
%nPart = size(coords,2);

                    f = 0.0001;

for part=1:nPart
%instead of looping over all pairs in the force calculation, find the distance between all
%pairs



dlist= sqrt(sum(bsxfun(@minus, coords(:,part), coords).^2, 1));

%index each set of nearest distances
[d, ind] = sort(rad(part,1)+rad'-dlist);
                
                %begin_index=find(d==2*rad(part,1));
                begin_index=find(d>0.0);
                
                %ind_closest = ind(2:(min(nPart,26))); %find the n nearest neighbors
                ind_closest=ind((begin_index):(end));
                
                coords_closest = coords(:,ind_closest);
                %rad_closest=rad(ind_closest,1);
                %poission_closest=poisson(partA,1)
                for partA=1:(size(ind_closest,2))
                
                dr =  coords(:,part) - coords_closest(:,partA);
                
                if norm(dr) > 0.0
                
                Rij = norm(dr);
                
                RijHat = dr/Rij;
                
                hij = (rad(ind_closest(partA),1) + rad(part,1) - Rij);
                
                Eij = ((1 - poisson(ind_closest(partA),1)^2)/(modulus(ind_closest(partA),1)) ...
                       + (1 - poisson(part,1)^2)/(modulus(part,1)));
                
                Rijf = (1/rad(part,1) + 1/rad(ind_closest(partA),1));
                
                %hij = max(0,h0);
                
                
                areaint = pi*(1/Rijf)*hij;
                
                
                invDr2 = (hij)^(3/2); % 1/r^2
                
                forceFact = (invDr2/(0.75*(Eij)*sqrt(Rijf)))-areaint*f...
                *0.5*(receptor(ind_closest(partA),1)*ligand(part,1)+ ...
                      receptor(part,1)*ligand(ind_closest(partA),1));
                
                % if (areaint > 0)
                pressure(part,1) = pressure(part,1)+ ...
                abs(forceFact/areaint);
                % end
                
                forces(:,part) = forces(:,part) + (forceFact*RijHat);
                
                end
                
                end
                
                
         
               
end

% figure
% h1=histogram(pressure,10, 'Normalization', 'probability','DisplayStyle','stairs')
% h1.BinWidth=1*10^(-4);
% hold on
% plot([5*10^(-4) 5*10^(-4)],ylim,'--','Linewidth',1)

%pressure_up_down=zeros(nPart,1);

                    centerM=zeros(3,1);
                    centerM(1,1)=mean(coords(1,:));
                    centerM(2,1)=mean(coords(2,:));
                    centerM(3,1)=mean(coords(3,:));
                    
                    %xpos=zeros(nPart,1);
                    %ypos=zeros(nPart,1);

                    %xpos(:,1)=coords(1,:)';
                    %ypos(:,1)=coords(2,:)';
                    %k = boundary(xpos,ypos);
                    %nBound=size(k,1);
                    
                    %deltar=0;
                    %for i=1:size(k,1)
                   %    deltar=deltar + sqrt((xpos(k(i),1)-centerM(1,1))^2+(ypos(k(i),1)-centerM(2,1))^2)/nBound;
                   % end
                    
                    %radius_tumor=deltar;
                    
                    for i=1:nPart
                        total_no_of_particles=total_no_of_particles+1;
                         r=norm(coords(:,i)-centerM);
                         % if(r<deltar)
                            pressure_profile(end+1,:)=[r pressure(i)];
                          %end
                        
                    end
                    
                    
                    
                    
                    
                    
                    

% figure
% 
% for i=1:nPart
% x=coords(1,i);
% y=coords(2,i);
% r=rad(i);
% %h = circle(x,y,r)
% hold on
% th = 0:pi/50:2*pi;
% xunit = r * cos(th) + x;
% yunit = r * sin(th) + y;
% plot(xunit, yunit);
% fill(xunit, yunit, pressure(i))
% %if(i<=numin)
% %    scatter(track(1,1:step,i),scatter(2,1:step,i));
%     
% %end
% end

%figure(2)
%h2=histogram(coordination_list,100, 'Normalization', 'probability')
%hold on
%figure

% for i=1:nPart
% x=coords(1,i);
% y=coords(2,i);
% r=rad(i);
% %h = circle(x,y,r)
% hold on
% th = 0:pi/50:2*pi;
% xunit = r * cos(th) + x;
% yunit = r * sin(th) + y;
% plot(xunit, yunit);
% fill(xunit, yunit, pressure_up_down(i))
% %if(i<=numin)
% %    scatter(track(1,1:step,i),scatter(2,1:step,i));
%     
% % %end
%  end

%fraction(frame-(start-1),1)=mean(pressure_up_down);
           % end
            
  %fraction_folders(:,count_er)=fraction;

%end


cd ..
end


% av_fraction=zeros(data_length,1);
% for i=1:data_length
%     av_fraction(i)=mean(fraction_folders(i,:));
% end


% gfac=0.6597;
% alpha=log(2)/(gfac*54000)-9.728*10^(-7);
% time=10000:10000:data_length*10^4;
% time=alpha*time;

%plot(time,av_fraction,'*-','LineWidth',1)

dl=max(pressure_profile(:,1))/10;
L=0:dl:max(pressure_profile(:,1));
p_r=zeros(2,size(L,2));
for i=1:size(pressure_profile,1)
index=round(pressure_profile(i,1)/dl)+1;
p_r(1,index)=p_r(1,index)+1;
p_r(2,index)=p_r(2,index)+pressure_profile(i,2);
end
for i=1:size(p_r,2)
p_r(2,i)=p_r(2,i)/p_r(1,i);
end
plot(L,p_r(2,:),'*','LineWidth',1)

                
                
               
                
  
