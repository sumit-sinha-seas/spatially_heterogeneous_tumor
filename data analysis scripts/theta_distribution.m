clear;
d=dir;
no_of_folders=1;
nsteps=40000;
%begin=5000;
cell_msd=zeros(nsteps,no_of_folders);
tracer_msd=zeros(nsteps,no_of_folders);
no_of_cells=zeros(nsteps,no_of_folders);
pressure=zeros(nsteps,no_of_folders);
pressure_mean=zeros(nsteps,no_of_folders-20);
straightness_tracer=zeros(1,1);
distance_cell=zeros(1,1);
displacement_cell=zeros(1,1);
theta_cell=zeros(1,1);
straightness_cell=zeros(1,1);
distance_tracer=zeros(1,1);
displacement_tracer=zeros(1,1);
theta_tracer=zeros(1,1);
%begin=1;

begin_count=26;
count_er=0;


t_o=1000;
interval=5400;
for i=begin_count:begin_count+(no_of_folders-1)
    i
    cd (d(i).name);
    load ('track.mat');
    %load ('fds_msdav.mat');
   % load ('no_cells.mat');
    %load ('pressure_history.mat');
    %if exist('mean_pressure.mat', 'file')
        
%         load('mean_pressure.mat');
%         count_er=count_er+1;
%         pressure_mean(:,count_er)=mean_pressure;
%     end
    
    no_of_cells=size(track,3);
    for part=1:no_of_cells
        part
        distance=0;
        
            for step=t_o:nsteps-2*interval  
          
                  vec1=track(:,step+interval,part)-track(:,step,part);
                  vec2=track(:,step+2*interval,part)-track(:,step+interval,part);
                  mag1=norm(vec1);
                  mag2=norm(vec2);
                  if((mag1>0) && (mag2>0))
                      theta_cell(end+1,1)= dot(vec1,vec2)/(mag1*mag2);
                      
                  end
                 %distance=distance+norm(track(:,step+1,100+part)-track(:,step,100+part));
                % msd_step=norm(track(:,step+1,100+part)-track(:,t_o,100+part))^2;
                 %cell_msd(step,i-(begin_count-1))=cell_msd(step,i-(begin_count-1))+msd_step/no_of_cells; 
         
            end
         
         
    end
   % for part=1:no_of_cells
         
%     for part=1:100
%         part
%          distance=0;
%         for step=t_o:nsteps-2*interval  
%              
%                    vec1=track(:,step+interval,part)-track(:,step,part);
%                   vec2=track(:,step+2*interval,part)-track(:,step+interval,part);
%                   mag1=norm(vec1);
%                   mag2=norm(vec2);
%                   if((mag1>0) && (mag2>0))
%                       theta_tracer(end+1,1)= dot(vec1,vec2)/(mag1*mag2);
%                       
%                   end
%                  
%                 % msd_step=norm(track(:,step+1,part)-track(:,t_o,part))^2;
%                %  tracer_msd(step,i-(begin_count-1))=tracer_msd(step,i-(begin_count-1))+msd_step/100; 
%          
%          end
%           
%     end
%     
%     no_of_cells(:,i-(begin_count-1))=no_cells; 
%     
%     for m=1:nsteps
%     pressure(m,i-(begin_count-1))=mean(pressure_history(1,m,:));
%     end

    
    cd ..
    
    
    
end
%tracer_msd(:,11)=[];
%tracer_msd(:,15)=[];
% av_cell_msd=zeros(nsteps,1);
% av_tracer_msd=zeros(nsteps,1);
% %av_no_of_cells=zeros(nsteps,1);
% %av_pressure=zeros(nsteps,1);
% %av_pressure_mean=zeros(nsteps,1);
% for i=1:nsteps
%     av_cell_msd(i,1)=mean(cell_msd(i,:));
%     av_tracer_msd(i,1)=mean(tracer_msd(i,:));
% %     av_no_of_cells(i,1)=mean(nonzeros(no_of_cells(i,:)));
% %     av_pressure(i,1)=mean(pressure(i,:));
% %     av_pressure_mean(i,1)=mean(pressure_mean(i,:));
% end
% 
% dt=10; 
% time=dt:dt:nsteps*dt;
% gfac=1;
% alpha=log(2)/(gfac*54000)-9.728*10^(-7);
% time1=alpha*time;
% alpha2=1/(gfac*54000);
% time2=alpha2*time;
% 
% %figure
% %loglog(time,av_c;endell_msd,'*-','LineWidth',3)
% loglog(time2(t_o:end),av_tracer_msd(t_o:end),'*-','LineWidth',3)
% hold on
% loglog(time2(t_o:end),av_cell_msd(t_o:end),'*-','LineWidth',3)

% figure(3)
% plot(time1,av_no_of_cells,'-','LineWidth',3)
% hold on
% figure(4)
% loglog(time,av_pressure,'-','LineWidth',0.5)
% hold on
% loglog(time,av_pressure_mean,'-','LineWidth',1)

    %if exist('mean_pressure.mat', 'file')
        
%         load('mean_pressure.mat');
%         count_er=count_er+1;
%         pressure_mean(:,count_er)=mean_pressure;
%     end
figure
 histogram(theta_cell,10, 'Normalization', 'probability') 
% hold on
 %histogram(theta_tracer,10, 'Normalization', 'probability') 
   