clear;
d1=dir;

%p_c=5*10^(-4);

no_of_folders=5;

initial_time=54000;
dt=500;
final_time=600000;
time=0:dt:final_time-initial_time;
av_msd_every_cell=zeros(size(time,2),no_of_folders);
st_interior_cell=zeros(1,1);
st_boundary_cell=zeros(1,1);
st_layer_1=zeros(1,1);
st_layer_2=zeros(1,1);
st_layer_3=zeros(1,1);
st_layer_4=zeros(1,1);
st_layer_5=zeros(1,1);
st_layer_6=zeros(1,1);
st_layer_7=zeros(1,1);
st_layer_8=zeros(1,1);
st_layer_9=zeros(1,1);
st_layer_10=zeros(1,1);
st_layer_extra=zeros(1,1);
distribution_si=zeros(1,2);
begin_count=17;
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


st_every_cell=zeros(1,size(intersecting_label,1));
 path_length=zeros(1,size(intersecting_label,1));
displacement=zeros(1,size(intersecting_label,1));
for i=1:size(intersecting_label,1)
  %  i
    
    track_particle=lifetime1(find(lifetime1(:,4)==intersecting_label(i,1)),:);
    start=find(track_particle(:,6)==initial_time);
    finish=find(track_particle(:,6)==final_time);
    track_particle=lifetime1(find(lifetime1(:,4)==intersecting_label(i,1)),1:3);
    req_track_particle=track_particle(start:finish,:);
    dis_com=norm(com-req_track_particle(end,:));
    
    for j=1:size(time,2)-1
        del_x=norm(req_track_particle(j+1,:)-req_track_particle(j,:));
        
        path_length(1,i)=path_length(1,i)+ del_x;
    end
    
    displacement(1,i)=norm(req_track_particle(end,:)-req_track_particle(1,:));
    st_every_cell(1,i)= displacement(1,i)/path_length(1,i);
    
    distribution_si(end+1,:)=[dis_com st_every_cell(1,i)];

     if(dis_com<30)
        st_interior_cell(1,end+1)=st_every_cell(1,i);
    end
    if(dis_com>60)
        st_boundary_cell(1,end+1)=st_every_cell(1,i);
    end
    if(dis_com<10)
        st_layer_1(1,end+1)=st_every_cell(1,i);
    end
    if(dis_com>10 && dis_com<20)
        st_layer_2(1,end+1)=st_every_cell(1,i);
    end
     if(dis_com>20 && dis_com<30)
        st_layer_3(1,end+1)=st_every_cell(1,i);
        
     end
      if(dis_com>30 && dis_com<40)
        st_layer_4(1,end+1)=st_every_cell(1,i);
        
      end
     
     if(dis_com>40 && dis_com<50)
        st_layer_5(1,end+1)=st_every_cell(1,i);
     end
     if(dis_com>50 && dis_com<60)
        st_layer_6(1,end+1)=st_every_cell(1,i);
     end
     if(dis_com>60 && dis_com<70)
        st_layer_7(1,end+1)=st_every_cell(1,i);
     end
     if(dis_com>70 && dis_com<80)
        st_layer_8(1,end+1)=st_every_cell(1,i);
     end
     if(dis_com>80 && dis_com<90)
        st_layer_9(1,end+1)=st_every_cell(1,i);
     end
     if(dis_com>90 && dis_com<100)
        st_layer_10(1,end+1)=st_every_cell(1,i);
     end
     if(dis_com>100)
        st_layer_extra(1,end+1)=st_every_cell(1,i);
    end
    
    
    
end
cd ..

end


st_layer_1(1,1)=[];
st_layer_2(1,1)=[];
st_layer_3(1,1)=[];
st_layer_4(1,1)=[];
st_layer_5(1,1)=[];
st_layer_6(1,1)=[];
st_layer_7(1,1)=[];
st_layer_8(1,1)=[];
st_layer_9(1,1)=[];
st_layer_10(1,1)=[];
st_layer_extra(1,1)=[];


av_st_layer_1=mean(st_layer_1(1,:));
av_st_layer_2=mean(st_layer_2(1,:));
av_st_layer_3=mean(st_layer_3(1,:));
av_st_layer_4=mean(st_layer_4(1,:));
av_st_layer_5=mean(st_layer_5(1,:));
av_st_layer_6=mean(st_layer_6(1,:));
av_st_layer_7=mean(st_layer_7(1,:));
av_st_layer_8=mean(st_layer_8(1,:));
av_st_layer_9=mean(st_layer_9(1,:));
av_st_layer_10=mean(st_layer_10(1,:));
av_st_layer_extra=mean(st_layer_extra(1,:));