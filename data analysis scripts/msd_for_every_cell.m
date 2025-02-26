clear;
d1=dir;

%p_c=5*10^(-4);

no_of_folders=5;

initial_time=54000*10;
dt=500;
final_time=54000*11;
time=0:dt:final_time-initial_time;
av_msd_every_cell=zeros(size(time,2),5);

begin_count=7;
count_er=0;

for looper=begin_count:begin_count+(no_of_folders-1)
    
    count_er=count_er+1;
    looper
 cd (d1(looper).name);

%clear

load('lifetime1.txt');



data_required_initial=lifetime1(find(lifetime1(:,6)==initial_time),:);
data_required_final=lifetime1(find(lifetime1(:,6)==final_time),:);

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
    
    for j=1:size(time,2)
        msd_every_cell(j,i)=(norm(req_track_particle(j,:)-req_track_particle(1,:)))^2;
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
loglog(time,av_av_msd_every_cell,'-','LineWidth',2)
