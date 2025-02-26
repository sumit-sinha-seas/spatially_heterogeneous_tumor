clear;
d1=dir;

%p_c=5*10^(-4);

no_of_folders=1;

initial_time=54000;
dt=500;
final_time=600000;
time=0:dt:final_time-initial_time;
av_msd_every_cell=zeros(size(time,2),no_of_folders);
vanhove_interior_cell=zeros(3,1);
vanhove_boundary_cell=zeros(3,1);
vanhove_layer_1=zeros(3,1);
vanhove_layer_2=zeros(3,1);
vanhove_layer_3=zeros(3,1);
vanhove_layer_4=zeros(3,1);
vanhove_layer_5=zeros(3,1);
vanhove_layer_6=zeros(3,1);
vanhove_layer_7=zeros(3,1);
vanhove_layer_8=zeros(3,1);
vanhove_layer_9=zeros(3,1);
vanhove_layer_10=zeros(3,1);
vanhove_layer_extra=zeros(3,1);
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


vanhove_every_cell=zeros(1,size(intersecting_label,1));
 path_length=zeros(1,size(intersecting_label,1));
displacement=zeros(1,size(intersecting_label,1));
del_t=1;
for i=1:size(intersecting_label,1)
  %  i
    
    track_particle=lifetime1(find(lifetime1(:,4)==intersecting_label(i,1)),:);
    start=find(track_particle(:,6)==initial_time);
    finish=find(track_particle(:,6)==final_time);
    track_particle=lifetime1(find(lifetime1(:,4)==intersecting_label(i,1)),1:3);
    req_track_particle=track_particle(start:finish,:);
    dis_com=norm(com-req_track_particle(end,:));
    
    for j=1:size(time,2)-del_t
        del_x=(req_track_particle(j+del_t,1)-req_track_particle(j,1));
        del_y=(req_track_particle(j+del_t,2)-req_track_particle(j,2));
        del_z=(req_track_particle(j+del_t,3)-req_track_particle(j,3));
        
       % path_length(1,i)=path_length(1,i)+ del_x;

    
             % displacement(1,i)=norm(req_track_particle(end,:)-req_track_particle(1,:));
             % vanhove_every_cell(1,i)= displacement(1,i)/path_length(1,i);

               if(dis_com<30)
                   vanhove_interior_cell(:,end+1)=[del_x;del_y;del_z];
               end
               if(dis_com>60)
                   vanhove_boundary_cell(:,end+1)=[del_x;del_y;del_z];
               end
               if(dis_com<10)
                    vanhove_layer_1(:,end+1)=[del_x;del_y;del_z];
               end
               if(dis_com>10 && dis_com<20)
                     vanhove_layer_2(:,end+1)=[del_x;del_y;del_z];
               end
               if(dis_com>20 && dis_com<30)
                      vanhove_layer_3(:,end+1)=[del_x;del_y;del_z];
        
               end
               if(dis_com>30 && dis_com<40)
                      vanhove_layer_4(:,end+1)=[del_x;del_y;del_z];
        
               end
     
               if(dis_com>40 && dis_com<50)
                    vanhove_layer_5(:,end+1)=[del_x;del_y;del_z];
               end
              if(dis_com>50 && dis_com<60)
                       vanhove_layer_6(:,end+1)=[del_x;del_y;del_z];
              end
                if(dis_com>60 && dis_com<70)
                       vanhove_layer_7(:,end+1)=[del_x;del_y;del_z];
                end
                 if(dis_com>70 && dis_com<80)
                       vanhove_layer_8(:,end+1)=[del_x;del_y;del_z];
                 end
              if(dis_com>80 && dis_com<90)
                     vanhove_layer_9(:,end+1)=[del_x;del_y;del_z];
              end
               if(dis_com>90 && dis_com<100)
                         vanhove_layer_10(:,end+1)=[del_x;del_y;del_z];
               end
                if(dis_com>100)
                         vanhove_layer_extra(:,end+1)=[del_x;del_y;del_z];
                end
    end
    
    
end
cd ..

end


vanhove_layer_1(:,1)=[];
vanhove_layer_2(:,1)=[];
vanhove_layer_3(:,1)=[];
vanhove_layer_4(:,1)=[];
vanhove_layer_5(:,1)=[];
vanhove_layer_6(:,1)=[];
vanhove_layer_7(:,1)=[];
vanhove_layer_8(:,1)=[];
vanhove_layer_9(:,1)=[];
vanhove_layer_10(:,1)=[];
vanhove_layer_extra(:,1)=[];
vanhove_interior_cell(:,1)=[];
vanhove_boundary_cell(:,1)=[];



% av_vanhove_layer_1=mean(vanhove_layer_1(1,:));
% av_vanhove_layer_2=mean(vanhove_layer_2(1,:));
% av_vanhove_layer_3=mean(vanhove_layer_3(1,:));
% av_vanhove_layer_4=mean(vanhove_layer_4(1,:));
% av_vanhove_layer_5=mean(vanhove_layer_5(1,:));
% av_vanhove_layer_6=mean(vanhove_layer_6(1,:));
% av_vanhove_layer_7=mean(vanhove_layer_7(1,:));
% av_vanhove_layer_8=mean(vanhove_layer_8(1,:));
% av_vanhove_layer_9=mean(vanhove_layer_9(1,:));
% av_vanhove_layer_10=mean(vanhove_layer_10(1,:));
% av_vanhove_layer_extra=mean(vanhove_layer_extra(1,:));