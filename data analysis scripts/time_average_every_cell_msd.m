clear;
d1=dir;

%p_c=5*10^(-4);

no_of_folders=5;

initial_time=54000*7;
dt=500;
final_time=54000*8;
time=0:dt:final_time-initial_time;
av_msd_every_cell=zeros(size(time,2),5);
exponents=zeros(1,1);
diffusion=zeros(1,1);
distance=zeros(1,1);
begin_count=8;
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
msd_every_cell=zeros(size(time,2),size(intersecting_label,1));

for i=1:size(intersecting_label,1)
    i
    
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

 %   MSDt(1:size(coordsx,1),1:no_cells_in_track)=zeros;
  %  count=0;
   % for j=1:1:1
        %count=count+1;
       % count
        for tau=1:size(coordsx,1)
            for k=1:size(coordsx,1)-tau
                     msd_every_cell(tau,i)=msd_every_cell(tau,i)+(1/(size(coordsx,1)-tau))*((coordsx(k+tau,1)-coordsx(k,1))^2+(coordsy(k+tau,1)-coordsy(k,1))^2 + (coordsz(k+tau,1)-coordsz(k,1))^2);
            end
        end
        
%         x1=time(2:20);
%         y1=msd_every_cell(2:20,i);
% 
%         [xData, yData] = prepareCurveData( x1', y1');
% 
% % Set up fittype and options.
%             ft = fittype( 'power1' );
%             opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
%             opts.Display = 'Off';
%             opts.StartPoint = [3.44222686723259 1.74147532981549];
% 
% % Fit model to data.
%             [fitresult, gof] = fit( xData, yData, ft, opts );
% 
% %fit_day_5(1,1,i)=fitresult.a;
%             diffusion(end+1)=fitresult.a;
%             exponents(end+1)=fitresult.b;
%             distance(end+1)=dis_com;
  %  end
%MSDtavg=mean(MSDt,2);
%time=10*step:10*step:size(coordsx,1)*10*step;
%time(end)
%size(time)
%loglog(time,MSDtavg,'--','LineWidth',1)
    
    
    
end

%av_msd_every_cell=zeros(size(time,2),1);

for i=1:size(time,2)
    av_msd_every_cell(i,count_er)=mean(msd_every_cell(i,:));
end

cd ..

end
av_av_msd_every_cell=zeros(size(time,2),1);
for i=1:size(time,2)
    av_av_msd_every_cell(i,1)=mean(av_msd_every_cell(i,:));
end
figure(1)
loglog(time/60,av_av_msd_every_cell,'--','LineWidth',1)
