clear;
d1=dir;

%p_c=5*10^(-4);

no_of_folders=50;

initial_time=54000;
dt=500;
final_time=600000;
time=0:dt:final_time-initial_time;
av_scat_every_cell=zeros(size(time,2),no_of_folders);
scat_interior_cell=zeros(size(time,2),1);
scat_boundary_cell=zeros(size(time,2),1);
scat_layer_1=zeros(size(time,2),1);
scat_layer_2=zeros(size(time,2),1);
scat_layer_3=zeros(size(time,2),1);
scat_layer_4=zeros(size(time,2),1);
scat_layer_5=zeros(size(time,2),1);
scat_layer_6=zeros(size(time,2),1);
scat_layer_7=zeros(size(time,2),1);
scat_layer_8=zeros(size(time,2),1);
scat_layer_9=zeros(size(time,2),1);
scat_layer_10=zeros(size(time,2),1);
scat_layer_extra=zeros(size(time,2),1);
begin_count=16;
interval=1;
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

for i=1:size(intersecting_label,1)
  %  i
    
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
    coords_particle=[coordsx coordsy coordsz];
 %   MSDt(1:size(coordsx,1),1:no_cells_in_track)=zeros;
  %  count=0;
   % for j=1:1:1
        %count=count+1;
       % count
        for tau=1:size(coordsx,1)
            for k=1:size(coordsx,1)-tau
                     mag_delta=norm(coords_particle(k+tau,:)-coords_particle(k,:));
                     scat_every_cell(tau,i)=scat_every_cell(tau,i)+(1/(size(coordsx,1)-tau))*(heaviside((length_scale)-mag_delta));
            end
        end
        
    if(dis_com<30)
        scat_interior_cell(:,end+1)=scat_every_cell(:,i);
    end
    if(dis_com>60)
        scat_boundary_cell(:,end+1)=scat_every_cell(:,i);
    end
    if(dis_com<10)
        scat_layer_1(:,end+1)=scat_every_cell(:,i);
    end
    if(dis_com>10 && dis_com<20)
        scat_layer_2(:,end+1)=scat_every_cell(:,i);
    end
     if(dis_com>20 && dis_com<30)
        scat_layer_3(:,end+1)=scat_every_cell(:,i);
        
     end
      if(dis_com>30 && dis_com<40)
        scat_layer_4(:,end+1)=scat_every_cell(:,i);
        
      end
     
     if(dis_com>40 && dis_com<50)
        scat_layer_5(:,end+1)=scat_every_cell(:,i);
     end
     if(dis_com>50 && dis_com<60)
        scat_layer_6(:,end+1)=scat_every_cell(:,i);
     end
     if(dis_com>60 && dis_com<70)
        scat_layer_7(:,end+1)=scat_every_cell(:,i);
     end
     if(dis_com>70 && dis_com<80)
        scat_layer_8(:,end+1)=scat_every_cell(:,i);
     end
     if(dis_com>80 && dis_com<90)
        scat_layer_9(:,end+1)=scat_every_cell(:,i);
     end
     if(dis_com>90 && dis_com<100)
        scat_layer_10(:,end+1)=scat_every_cell(:,i);
     end
     if(dis_com>100)
        scat_layer_extra(:,end+1)=scat_every_cell(:,i);
    end
     
    
    
    
    
end
for i=1:size(time,2)
    av_scat_every_cell(i,count_er)=mean(scat_every_cell(i,:));
end
cd ..

end


av_av_scat_every_cell=zeros(size(time,2),1);
for i=1:size(time,2)
    av_av_scat_every_cell(i,1)=mean(av_scat_every_cell(i,:));
end
figure(1)
loglog(time,av_av_scat_every_cell,'-','LineWidth',2)
hold on
scat_boundary_cell(:,1)=[];
scat_interior_cell(:,1)=[];
av_scat_boundary_cell=zeros(size(time,2),1);
for i=1:size(time,2)
    av_scat_boundary_cell(i,1)=mean(scat_boundary_cell(i,:));
end
loglog(time,av_scat_boundary_cell,'-','LineWidth',2)
av_scat_interior_cell=zeros(size(time,2),1);
for i=1:size(time,2)
av_scat_interior_cell(i,1)=mean(scat_interior_cell(i,:));
end
loglog(time,av_scat_interior_cell,'-','LineWidth',2)

hold off

scat_layer_1(:,1)=[];
scat_layer_2(:,1)=[];
scat_layer_3(:,1)=[];
scat_layer_4(:,1)=[];
scat_layer_5(:,1)=[];
scat_layer_6(:,1)=[];
scat_layer_7(:,1)=[];
scat_layer_8(:,1)=[];
scat_layer_9(:,1)=[];
scat_layer_10(:,1)=[];
scat_layer_extra(:,1)=[];
av_scat_layer_1=zeros(size(time,2),1);
av_scat_layer_2=zeros(size(time,2),1);
av_scat_layer_3=zeros(size(time,2),1);
av_scat_layer_4=zeros(size(time,2),1);
av_scat_layer_5=zeros(size(time,2),1);
av_scat_layer_6=zeros(size(time,2),1);
av_scat_layer_7=zeros(size(time,2),1);
av_scat_layer_8=zeros(size(time,2),1);
av_scat_layer_9=zeros(size(time,2),1);
av_scat_layer_10=zeros(size(time,2),1);
av_scat_layer_extra=zeros(size(time,2),1);


for i=1:size(time,2)
av_scat_layer_1(i,1)=mean(scat_layer_1(i,:));
av_scat_layer_2(i,1)=mean(scat_layer_2(i,:));
av_scat_layer_3(i,1)=mean(scat_layer_3(i,:));
av_scat_layer_4(i,1)=mean(scat_layer_4(i,:));
av_scat_layer_5(i,1)=mean(scat_layer_5(i,:));
av_scat_layer_6(i,1)=mean(scat_layer_6(i,:));
av_scat_layer_7(i,1)=mean(scat_layer_7(i,:));
av_scat_layer_8(i,1)=mean(scat_layer_8(i,:));
av_scat_layer_9(i,1)=mean(scat_layer_9(i,:));
av_scat_layer_10(i,1)=mean(scat_layer_10(i,:));
av_scat_layer_extra(i,1)=mean(scat_layer_extra(i,:));
end
figure(2)
plot(time(1:end-1),av_scat_layer_1(1:end-1),'-','LineWidth',2)
hold on
plot(time(1:end-1),av_scat_layer_2(1:end-1),'-','LineWidth',2)
plot(time(1:end-1),av_scat_layer_3(1:end-1),'-','LineWidth',2)
plot(time(1:end-1),av_scat_layer_4(1:end-1),'-','LineWidth',2)
plot(time(1:end-1),av_scat_layer_5(1:end-1),'-','LineWidth',2)
plot(time(1:end-1),av_scat_layer_6(1:end-1),'-','LineWidth',2)
plot(time(1:end-1),av_scat_layer_7(1:end-1),'-','LineWidth',2)
plot(time(1:end-1),av_scat_layer_8(1:end-1),'-','LineWidth',2)
plot(time(1:end-1),av_scat_layer_9(1:end-1),'-','LineWidth',2)
plot(time(1:end-1),av_scat_layer_10(1:end-1),'-','LineWidth',2)
plot(time(1:end-1),av_scat_layer_extra(1:end-1),'-','LineWidth',2)