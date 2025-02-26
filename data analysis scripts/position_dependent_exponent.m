clear;
d1=dir;

%p_c=5*10^(-4);

no_of_folders=100;

initial_time=54000;
dt=500;
final_time=600000;
time=0:dt:final_time-initial_time;
av_msd_every_cell=zeros(size(time,2),no_of_folders);
msd_interior_cell=zeros(size(time,2),1);
msd_boundary_cell=zeros(size(time,2),1);
msd_layer_1=zeros(size(time,2),1);
msd_layer_2=zeros(size(time,2),1);
msd_layer_3=zeros(size(time,2),1);
msd_layer_4=zeros(size(time,2),1);
msd_layer_5=zeros(size(time,2),1);
msd_layer_6=zeros(size(time,2),1);
msd_layer_7=zeros(size(time,2),1);
msd_layer_8=zeros(size(time,2),1);
msd_layer_9=zeros(size(time,2),1);
msd_layer_10=zeros(size(time,2),1);
msd_layer_extra=zeros(size(time,2),1);
begin_count=20;
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


msd_every_cell=zeros(size(time,2),size(intersecting_label,1));

for i=1:size(intersecting_label,1)
  %  i
    
    track_particle=lifetime1(find(lifetime1(:,4)==intersecting_label(i,1)),:);
    start=find(track_particle(:,6)==initial_time);
    finish=find(track_particle(:,6)==final_time);
    track_particle=lifetime1(find(lifetime1(:,4)==intersecting_label(i,1)),1:3);
    req_track_particle=track_particle(start:finish,:);
    dis_com=norm(com-req_track_particle(end,:));
    
    for j=1:size(time,2)
        msd_every_cell(j,i)=(norm(req_track_particle(j,:)-req_track_particle(1,:)))^2;
    end
    if(dis_com<30)
        msd_interior_cell(:,end+1)=msd_every_cell(:,i);
    end
    if(dis_com>60)
        msd_boundary_cell(:,end+1)=msd_every_cell(:,i);
    end
    if(dis_com<25)
        msd_layer_1(:,end+1)=msd_every_cell(:,i);
    end
    if(dis_com>25 && dis_com<50)
        msd_layer_2(:,end+1)=msd_every_cell(:,i);
    end
     if(dis_com>50 && dis_com<75)
        msd_layer_3(:,end+1)=msd_every_cell(:,i);
        
     end
      if(dis_com>75 && dis_com<100)
        msd_layer_4(:,end+1)=msd_every_cell(:,i);
        
      end
     
%      if(dis_com>40 && dis_com<50)
%         msd_layer_5(:,end+1)=msd_every_cell(:,i);
%      end
%      if(dis_com>50 && dis_com<60)
%         msd_layer_6(:,end+1)=msd_every_cell(:,i);
%      end
%      if(dis_com>60 && dis_com<70)
%         msd_layer_7(:,end+1)=msd_every_cell(:,i);
%      end
%      if(dis_com>70 && dis_com<80)
%         msd_layer_8(:,end+1)=msd_every_cell(:,i);
%      end
%      if(dis_com>80 && dis_com<90)
%         msd_layer_9(:,end+1)=msd_every_cell(:,i);
%      end
%      if(dis_com>90 && dis_com<100)
%         msd_layer_10(:,end+1)=msd_every_cell(:,i);
%      end
     if(dis_com>100)
        msd_layer_5(:,end+1)=msd_every_cell(:,i);
    end
     
    
    
    
end

% figure
% for i=1:size(intersecting_label,1)
% loglog(time,msd_every_cell(:,i),'--','LineWidth',1)
% hold on
% end
% 

%figure
%av_msd_every_cell=zeros(size(time,2),1);

for i=1:size(time,2)
    av_msd_every_cell(i,count_er)=mean(msd_every_cell(i,:));
end
%loglog(time,av_msd_every_cell,'*-','LineWidth',2)


cd ..
end

av_av_msd_every_cell=zeros(size(time,2),1);
for i=1:size(time,2)
    av_av_msd_every_cell(i,1)=mean(av_msd_every_cell(i,:));
end
figure(1)
loglog(time,av_av_msd_every_cell,'-','LineWidth',2)
hold on
msd_boundary_cell(:,1)=[];
msd_interior_cell(:,1)=[];
av_msd_boundary_cell=zeros(size(time,2),1);
for i=1:size(time,2)
    av_msd_boundary_cell(i,1)=mean(msd_boundary_cell(i,:));
end
loglog(time,av_msd_boundary_cell,'-','LineWidth',2)
av_msd_interior_cell=zeros(size(time,2),1);
for i=1:size(time,2)
av_msd_interior_cell(i,1)=mean(msd_interior_cell(i,:));
end
loglog(time,av_msd_interior_cell,'-','LineWidth',2)

hold off


msd_layer_1(:,1)=[];
msd_layer_2(:,1)=[];
msd_layer_3(:,1)=[];
msd_layer_4(:,1)=[];
msd_layer_5(:,1)=[];
msd_layer_6(:,1)=[];
msd_layer_7(:,1)=[];
msd_layer_8(:,1)=[];
msd_layer_9(:,1)=[];
msd_layer_10(:,1)=[];
msd_layer_extra(:,1)=[];
av_msd_layer_1=zeros(size(time,2),1);
av_msd_layer_2=zeros(size(time,2),1);
av_msd_layer_3=zeros(size(time,2),1);
av_msd_layer_4=zeros(size(time,2),1);
av_msd_layer_5=zeros(size(time,2),1);
av_msd_layer_6=zeros(size(time,2),1);
av_msd_layer_7=zeros(size(time,2),1);
av_msd_layer_8=zeros(size(time,2),1);
av_msd_layer_9=zeros(size(time,2),1);
av_msd_layer_10=zeros(size(time,2),1);
av_msd_layer_extra=zeros(size(time,2),1);


for i=1:size(time,2)
av_msd_layer_1(i,1)=mean(msd_layer_1(i,:));
av_msd_layer_2(i,1)=mean(msd_layer_2(i,:));
av_msd_layer_3(i,1)=mean(msd_layer_3(i,:));
av_msd_layer_4(i,1)=mean(msd_layer_4(i,:));
av_msd_layer_5(i,1)=mean(msd_layer_5(i,:));
av_msd_layer_6(i,1)=mean(msd_layer_6(i,:));
av_msd_layer_7(i,1)=mean(msd_layer_7(i,:));
av_msd_layer_8(i,1)=mean(msd_layer_8(i,:));
av_msd_layer_9(i,1)=mean(msd_layer_9(i,:));
av_msd_layer_10(i,1)=mean(msd_layer_10(i,:));
av_msd_layer_extra(i,1)=mean(msd_layer_extra(i,:));
end
figure(2)
loglog(time,av_msd_layer_1,'-','LineWidth',2)
hold on
loglog(time,av_msd_layer_2,'-','LineWidth',2)
loglog(time,av_msd_layer_3,'-','LineWidth',2)
loglog(time,av_msd_layer_4,'-','LineWidth',2)
loglog(time,av_msd_layer_5,'-','LineWidth',2)
loglog(time,av_msd_layer_6,'-','LineWidth',2)
loglog(time,av_msd_layer_7,'-','LineWidth',2)
loglog(time,av_msd_layer_8,'-','LineWidth',2)
loglog(time,av_msd_layer_9,'-','LineWidth',2)
loglog(time,av_msd_layer_10,'-','LineWidth',2)
loglog(time,av_msd_layer_extra,'-','LineWidth',2)










