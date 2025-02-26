clear;
d1=dir;

%p_c=5*10^(-4);

no_of_folders=10;

initial_time=500000;
dt=500;
final_time=600000;
time=0:dt:final_time-initial_time;
%av_scat_every_cell=zeros(size(time,2),no_of_folders);
%scat_interior_cell=zeros(size(time,2),1);
%scat_boundary_cell=zeros(size(time,2),1);
xi_4_layer_1=zeros(size(time,2),no_of_folders);
xi_4_layer_2=zeros(size(time,2),no_of_folders);
xi_4_layer_3=zeros(size(time,2),no_of_folders);
xi_4_layer_4=zeros(size(time,2),no_of_folders);
xi_4_layer_5=zeros(size(time,2),no_of_folders);
xi_4_layer_6=zeros(size(time,2),no_of_folders);
xi_4_layer_7=zeros(size(time,2),no_of_folders);
xi_4_layer_8=zeros(size(time,2),no_of_folders);
xi_4_layer_9=zeros(size(time,2),no_of_folders);
xi_4_layer_10=zeros(size(time,2),no_of_folders);

begin_count=14;
interval=10;
t_end=size(time,2)-1;
length_scale=3.33;
count_er=0;


for looper=begin_count:begin_count+(no_of_folders-1)
    
    count_er=count_er+1;
    looper
 cd (d1(looper).name);

%clear

load('lifetime1.txt');



data_required_initial=lifetime1(find(lifetime1(:,6)==initial_time),:);
data_required_final=lifetime1(find(lifetime1(:,6)==final_time),:);
com=zeros(1,3);
com(1,1)=mean(data_required_final(:,1));
com(1,2)=mean(data_required_final(:,2));
com(1,3)=mean(data_required_final(:,2));


label_initial=data_required_initial(:,4);
label_final=data_required_final(:,4);

intersecting_label=intersect(label_initial,label_final);


scat_every_cell=zeros(size(time,2),size(intersecting_label,1));
%for tau=0:interval:size(time,2)-1
for tau=0:interval:t_end
    tau
    data_xi_4_layer1=zeros(1,1);
     data_xi_4_layer2=zeros(1,1);
      data_xi_4_layer3=zeros(1,1);
       data_xi_4_layer4=zeros(1,1);
        data_xi_4_layer5=zeros(1,1);
         data_xi_4_layer6=zeros(1,1);
          data_xi_4_layer7=zeros(1,1);
           data_xi_4_layer8=zeros(1,1);
            data_xi_4_layer9=zeros(1,1);
             data_xi_4_layer10=zeros(1,1);
             no_of_part=zeros(10,1);
      
       
            
for i=1:size(intersecting_label,1)
    
     track_particle=lifetime1(find(lifetime1(:,4)==intersecting_label(i,1)),:);
    start=find(track_particle(:,6)==initial_time);
    finish=find(track_particle(:,6)==final_time);
    track_particle=lifetime1(find(lifetime1(:,4)==intersecting_label(i,1)),1:3);
    req_track_particle=track_particle(start:finish,:);
    dis_com=norm(com-req_track_particle(end,:));
    if(dis_com<100)
    no_of_part(floor(dis_com/10)+1)= no_of_part(floor(dis_com/10)+1)+1;
    end
    
      step=1;
    count=0;
    coordsx=zeros(size(req_track_particle,1),1);
    coordsy=zeros(size(req_track_particle,1),1);
    coordsz=zeros(size(req_track_particle,1),1);
    for j=1:step:size(req_track_particle,1)
        count=count+1;
        coordsx(count,1)=req_track_particle(j,1);
        coordsy(count,1)=req_track_particle(j,2);
        coordsz(count,1)=req_track_particle(j,3);
    end
    coords_particle=[coordsx coordsy coordsz];
  %  i
        for k=1:size(coordsx,1)-tau
    
   
 %   MSDt(1:size(coordsx,1),1:no_cells_in_track)=zeros;
  %  count=0;
   % for j=1:1:1
        %count=count+1;
       % count
       % for tau=1:size(coordsx,1)
           % for k=1:size(coordsx,1)-tau
                     mag_delta=norm(coords_particle(k+tau,:)-coords_particle(k,:));
                   %  data_xi_4_k(end+1,1)= (heaviside((length_scale)-mag_delta));
                    % scat_every_cell(tau,i)=scat_every_cell(tau,i)+(1/(size(coordsx,1)-tau))*(heaviside((length_scale)-mag_delta));
           % end
       % end
        
        
   
         %if(dis_com<10)
          %    data_xi_4_layer1(end+1,1)=(heaviside((length_scale)-mag_delta));
        % end
         if( dis_com<20)
             data_xi_4_layer2(end+1,1)=(heaviside((length_scale)-mag_delta));
         end
         if(dis_com>20 && dis_com<30)
            data_xi_4_layer3(end+1,1)=(heaviside((length_scale)-mag_delta));
        
         end
         if(dis_com>30 && dis_com<40)
             data_xi_4_layer4(end+1,1)=(heaviside((length_scale)-mag_delta));
        
         end
     
         if(dis_com>40 && dis_com<50)
            data_xi_4_layer5(end+1,1)=(heaviside((length_scale)-mag_delta));
         end
        if(dis_com>50 && dis_com<60)
            data_xi_4_layer6(end+1,1)=(heaviside((length_scale)-mag_delta));
        end
        if(dis_com>60 && dis_com<70)
            data_xi_4_layer7(end+1,1)=(heaviside((length_scale)-mag_delta));
        end
        if(dis_com>70 && dis_com<80)
             data_xi_4_layer8(end+1,1)=(heaviside((length_scale)-mag_delta));
        end
        if(dis_com>80 && dis_com<100)
            data_xi_4_layer9(end+1,1)=(heaviside((length_scale)-mag_delta));
        end
      %  if(dis_com>90 && dis_com<100)
           %  data_xi_4_layer10(end+1,1)=(heaviside((length_scale)-mag_delta));
       % end
    
    
     
     end
     
    
    
    
    
end

% data_xi_4_layer1(1)=[];
     data_xi_4_layer2(1)=[];
      data_xi_4_layer3(1)=[];
       data_xi_4_layer4(1)=[];
        data_xi_4_layer5(1)=[];
         data_xi_4_layer6(1)=[];
          data_xi_4_layer7(1)=[];
           data_xi_4_layer8(1)=[];
            data_xi_4_layer9(1)=[];
          %   data_xi_4_layer10(1)=[];
             
           %  xi_4_layer_1(tau,1)=no_of_part(1)*var(data_xi_4_layer1(:,1));
             xi_4_layer_2(tau+1,count_er)=(no_of_part(1)+no_of_part(2))*var(data_xi_4_layer2(:,1));
             xi_4_layer_3(tau+1,count_er)=no_of_part(3)*var(data_xi_4_layer3(:,1));
             xi_4_layer_4(tau+1,count_er)=no_of_part(4)*var(data_xi_4_layer4(:,1));
             xi_4_layer_5(tau+1,count_er)=no_of_part(5)*var(data_xi_4_layer5(:,1));
             xi_4_layer_6(tau+1,count_er)=no_of_part(6)*var(data_xi_4_layer6(:,1));
             xi_4_layer_7(tau+1,count_er)=no_of_part(7)*var(data_xi_4_layer7(:,1));
             xi_4_layer_8(tau+1,count_er)=no_of_part(8)*var(data_xi_4_layer8(:,1));
             xi_4_layer_9(tau+1,count_er)=(no_of_part(9)+no_of_part(10))*var(data_xi_4_layer9(:,1));
           %  xi_4_layer_10(tau,1)=no_of_part(10)*var(data_xi_4_layer10(:,1));

       
end
  
   
%for i=1:size(time,2)
 %   av_scat_every_cell(i,count_er)=mean(scat_every_cell(i,:));
%end
 cd ..   

end


%figure(1)
%plot(time(1:50:end-1),xi_4_layer_1(1:50:end-1),'-','LineWidth',2)

plot(time(1:interval:t_end),mean(xi_4_layer_2(1:interval:t_end,:),2),'-','LineWidth',2)
set(gca,'xscale','log');
hold on
plot(time(1:interval:t_end),mean(xi_4_layer_3(1:interval:t_end,:),2),'-','LineWidth',2)
plot(time(1:interval:t_end),mean(xi_4_layer_4(1:interval:t_end,:),2),'-','LineWidth',2)
plot(time(1:interval:t_end),mean(xi_4_layer_5(1:interval:t_end,:),2),'-','LineWidth',2)
plot(time(1:interval:t_end),mean(xi_4_layer_6(1:interval:t_end,:),2),'-','LineWidth',2)
plot(time(1:interval:t_end),mean(xi_4_layer_7(1:interval:t_end,:),2),'-','LineWidth',2)
plot(time(1:interval:t_end),mean(xi_4_layer_8(1:interval:t_end,:),2),'-','LineWidth',2)
plot(time(1:interval:t_end),mean(xi_4_layer_9(1:interval:t_end,:),2),'-','LineWidth',2)