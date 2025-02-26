clear;
d1=dir;

%p_c=5*10^(-4);

no_of_folders=5;

initial_time=54000*10;
dt=500;
final_time=600000;
time=0:dt:final_time-initial_time;
av_vel_corr_every_cell=zeros(size(time,2)-1,no_of_folders);
vel_corr_interior_cell=zeros(size(time,2)-1,1);
vel_corr_boundary_cell=zeros(size(time,2)-1,1);
vel_corr_layer_1=zeros(size(time,2)-1,1);
vel_corr_layer_2=zeros(size(time,2)-1,1);
vel_corr_layer_3=zeros(size(time,2)-1,1);
vel_corr_layer_4=zeros(size(time,2)-1,1);
vel_corr_layer_5=zeros(size(time,2)-1,1);
vel_corr_layer_6=zeros(size(time,2)-1,1);
vel_corr_layer_7=zeros(size(time,2)-1,1);
vel_corr_layer_8=zeros(size(time,2)-1,1);
vel_corr_layer_9=zeros(size(time,2)-1,1);
vel_corr_layer_10=zeros(size(time,2)-1,1);
vel_corr_layer_extra=zeros(size(time,2)-1,1);
exponents=zeros(1,1);
diffusion=zeros(1,1);
distance=zeros(1,1);
begin_count=16;
count_er=0;

for looper=begin_count:begin_count+(no_of_folders-1)
    
    count_er=count_er+1;
    looper
 cd (d1(looper).name);

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

time=0:dt:final_time-initial_time;
vel_corr_every_cell=zeros(size(time,2)-1,size(intersecting_label,1));

for i=1:size(intersecting_label,1)
    
    
    track_particle=lifetime1(find(lifetime1(:,4)==intersecting_label(i,1)),:);
    start=find(track_particle(:,6)==initial_time);
    finish=find(track_particle(:,6)==final_time);
    track_particle=lifetime1(find(lifetime1(:,4)==intersecting_label(i,1)),1:3);
    req_track_particle=track_particle(start:finish,:);
    dis_com=norm(com-req_track_particle(end,:));
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
    vel=zeros(size(time,2)-1,3);
    for m=1:size(time,2)-1
        vel(m,:)=(req_track_particle(m+1,:)-req_track_particle(m,:))/dt;
    end

 %   MSDt(1:size(coordsx,1),1:no_cells_in_track)=zeros;
  %  count=0;
   % for j=1:1:1
        %count=count+1;
       % count
        for tau=0:size(vel,1)-1
            for k=1:size(vel,1)-tau
                     mag1=norm(vel(k+tau,:));
                     mag2=norm(vel(k,:));
                     vel_corr_every_cell(tau+1,i)= vel_corr_every_cell(tau+1,i)+(1/(size(vel,1)-tau))*dot(vel(k+tau,:),vel(k,:))/(mag1*mag2);
            end
        end
      %  vel_corr_every_cell(:,i)=vel_corr_every_cell(:,i)/vel_corr_every_cell(1,i);
        
        if(dis_com<30)
                vel_corr_interior_cell(:,end+1)=vel_corr_every_cell(:,i);
        end
        if(dis_com>60)
                vel_corr_boundary_cell(:,end+1)=vel_corr_every_cell(:,i);
        end
        if(dis_com<10)
                vel_corr_layer_1(:,end+1)=vel_corr_every_cell(:,i);
        end
         if(dis_com>10 && dis_com<20)
                vel_corr_layer_2(:,end+1)=vel_corr_every_cell(:,i);
         end
         if(dis_com>20 && dis_com<30)
                vel_corr_layer_3(:,end+1)=vel_corr_every_cell(:,i);
         end
         if(dis_com>30 && dis_com<40)
                vel_corr_layer_4(:,end+1)=vel_corr_every_cell(:,i);
         end
         if(dis_com>40 && dis_com<50)
                vel_corr_layer_5(:,end+1)=vel_corr_every_cell(:,i);
         end
         if(dis_com>50 && dis_com<60)
                vel_corr_layer_6(:,end+1)=vel_corr_every_cell(:,i);
         end
         if(dis_com>60 && dis_com<70)
                vel_corr_layer_7(:,end+1)=vel_corr_every_cell(:,i);
         end
         if(dis_com>70 && dis_com<80)
                vel_corr_layer_8(:,end+1)=vel_corr_every_cell(:,i);
         end
         if(dis_com>80 && dis_com<90)
                vel_corr_layer_9(:,end+1)=vel_corr_every_cell(:,i);
         end
         if(dis_com>90 && dis_com<100)
                vel_corr_layer_10(:,end+1)=vel_corr_every_cell(:,i);
         end
         if(dis_com>100)
                vel_corr_layer_extra(:,end+1)=vel_corr_every_cell(:,i);
         end
         
         
        
        
     
        
        
    
    
    
end

%av_msd_every_cell=zeros(size(time,2),1);

for i=1:size(time,2)-1
    av_vel_corr_every_cell(i,count_er)=mean(vel_corr_every_cell(i,:));
end

cd ..

end
av_av_vel_corr_every_cell=zeros(size(time,2)-1,1);
for i=1:size(time,2)-1
    av_av_vel_corr_every_cell(i,1)=mean(av_vel_corr_every_cell(i,:));
end
time=time(1:end-1);
figure(1)
loglog(time,av_av_vel_corr_every_cell,'--','LineWidth',1)

hold on
vel_corr_boundary_cell(:,1)=[];
vel_corr_interior_cell(:,1)=[];
av_vel_corr_boundary_cell=zeros(size(time,2)-1,1);
for i=1:size(time,2)
    av_vel_corr_boundary_cell(i,1)=mean(vel_corr_boundary_cell(i,:));
end
loglog(time,av_vel_corr_boundary_cell,'--','LineWidth',2)
av_vel_corr_interior_cell=zeros(size(time,2)-1,1);
for i=1:size(time,2)
av_vel_corr_interior_cell(i,1)=mean(vel_corr_interior_cell(i,:));
end
loglog(time,av_vel_corr_interior_cell,'--','LineWidth',2)
hold off





vel_corr_layer_1(:,1)=[];
vel_corr_layer_2(:,1)=[];
vel_corr_layer_3(:,1)=[];
vel_corr_layer_4(:,1)=[];
vel_corr_layer_5(:,1)=[];
vel_corr_layer_6(:,1)=[];
vel_corr_layer_7(:,1)=[];
vel_corr_layer_8(:,1)=[];
vel_corr_layer_9(:,1)=[];
vel_corr_layer_10(:,1)=[];
vel_corr_layer_extra(:,1)=[];
av_vel_corr_layer_1=zeros(size(time,2)-1,1);
av_vel_corr_layer_2=zeros(size(time,2)-1,1);
av_vel_corr_layer_3=zeros(size(time,2)-1,1);
av_vel_corr_layer_4=zeros(size(time,2)-1,1);
av_vel_corr_layer_5=zeros(size(time,2)-1,1);
av_vel_corr_layer_6=zeros(size(time,2)-1,1);
av_vel_corr_layer_7=zeros(size(time,2)-1,1);
av_vel_corr_layer_8=zeros(size(time,2)-1,1);
av_vel_corr_layer_9=zeros(size(time,2)-1,1);
av_vel_corr_layer_10=zeros(size(time,2)-1,1);
av_vel_corr_layer_extra=zeros(size(time,2)-1,1);


for i=1:size(time,2)
av_vel_corr_layer_1(i,1)=mean(vel_corr_layer_1(i,:));
av_vel_corr_layer_2(i,1)=mean(vel_corr_layer_2(i,:));
av_vel_corr_layer_3(i,1)=mean(vel_corr_layer_3(i,:));
av_vel_corr_layer_4(i,1)=mean(vel_corr_layer_4(i,:));
av_vel_corr_layer_5(i,1)=mean(vel_corr_layer_5(i,:));
av_vel_corr_layer_6(i,1)=mean(vel_corr_layer_6(i,:));
av_vel_corr_layer_7(i,1)=mean(vel_corr_layer_7(i,:));
av_vel_corr_layer_8(i,1)=mean(vel_corr_layer_8(i,:));
av_vel_corr_layer_9(i,1)=mean(vel_corr_layer_9(i,:));
av_vel_corr_layer_10(i,1)=mean(vel_corr_layer_10(i,:));
av_vel_corr_layer_extra(i,1)=mean(vel_corr_layer_extra(i,:));
end
figure(2)
loglog(time,av_vel_corr_layer_1,'-','LineWidth',2)
hold on
loglog(time,av_vel_corr_layer_2,'-','LineWidth',2)
loglog(time,av_vel_corr_layer_3,'-','LineWidth',2)
loglog(time,av_vel_corr_layer_4,'-','LineWidth',2)
loglog(time,av_vel_corr_layer_5,'-','LineWidth',2)
loglog(time,av_vel_corr_layer_6,'-','LineWidth',2)
loglog(time,av_vel_corr_layer_7,'-','LineWidth',2)
loglog(time,av_vel_corr_layer_8,'-','LineWidth',2)
loglog(time,av_vel_corr_layer_9,'-','LineWidth',2)
loglog(time,av_vel_corr_layer_10,'-','LineWidth',2)
loglog(time,av_vel_corr_layer_extra,'-','LineWidth',2)



