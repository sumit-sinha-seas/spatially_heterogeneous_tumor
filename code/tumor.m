
 


load lebedev.txt;
b1=lebedev(:,1);
b2=lebedev(:,2);
b3=lebedev(:,3);
rng('shuffle');

nSteps = 40000;
Polar=zeros(nSteps,1);

cyclenum = 1;
scat_fun1=zeros(nSteps,1);  
scat_fun=zeros(nSteps,cyclenum);

radG=zeros(nSteps,cyclenum); 

den=zeros(nSteps,cyclenum);

df=zeros(nSteps,cyclenum);

ds_msdav=zeros(nSteps,1);

tds_msdav=zeros(nSteps,1);

radG_av=zeros(nSteps,1);

den_av=zeros(nSteps,1);

df_av=zeros(nSteps,1);

%numcells=0;

%scat_fun=zeros(nSteps,1);



liferows=1;
visdata=zeros(liferows,15);
visdata_row=zeros(liferows,15);
qt=zeros(nSteps,cyclenum);
tqt=zeros(nSteps,cyclenum);

 for cyclemsd=1:cyclenum
    
   numin=200;
   numin0=200;
   nPart = numin;        
   density = 0.001;
   mass = 1;          
   nDim = 3;           
    
    
   dt = 10.0;      
   dt2 = dt*dt;       
    
   
   sampleFreq = 100; 
   sampleCounter = 1; 
   printFreq = 1000;   
   plotFreq = 100;  
   radmitosis = 5.0; 
   samplecounter = 1;
   eta=0.005; 
   gfac=1;
   taumin = gfac*54000;
   p_c=1.7*10^(-6); %check for ik pressure definition 
   [coords L] = initCubicGrid(nPart,density);
   vels = zeros(nDim,nPart);
   rad = zeros(nPart,1);
   modulus = zeros(nPart,1);
   poisson = zeros(nPart,1);
   receptor = zeros(nPart,1);
   ligand = zeros(nPart,1);
   %lifetime = zeros(nPart,1);
   label= zeros(nPart,1);
     lifetime=zeros(nPart,1);
   countdr=0; 
   
    track(nDim,nSteps+1,numin0)=zeros;
    
    
    for part = 1:nPart
        rad(part,1) = randgaussrad(4.5,0.5);  
        modulus(part,1) = randgaussrad(10^-3,10^-4);
        poisson(part,1) = randgaussrad(0.5,0.02);
        receptor(part,1) = randgaussrad(1.0,0.02);
        ligand(part,1) = randgaussrad(1.0,0.02);
        lifetime(part,1) = 0;
        label(part,1)=part;
    end
    
    
%    for part=1:numin0/2
        
 %       rad(part,1)=4.5;
        
  %  end
  
    
 
    
    time = 0; 
    
    volrate = (2*pi*(radmitosis)^3)/(3*taumin);
    


    initial=coords(:,1:end);
    force_cell_tracer=zeros(3,nSteps,nPart); % added to record the force on cell and tps
    no_of_cells_ever_involved =  nPart;  
    
    for avini=1:numin
         track(:,1,avini)=coords(:,avini);
    end
    
    for step = 1:nSteps
       
        centerM = zeros(3,1);
        
 

        [forces,gamma3,pressure] = LJ_Force(coords,rad,poisson,modulus,...
            nPart,receptor,ligand);
         
        
        %[forces,gamma3,pressure] = LJ_Force(coords,rad,poisson,modulus,...
         %   nPart,receptor,ligand);


        
     
        
        
        gammat = 0;

        for part =1:nPart
            
            gammat = 6*pi*eta*rad(part,1) + gamma3(part,1);
            coords(:,part) = coords(:,part) + (dt*forces(:,part))/(gammat);
            vels(:,part) = forces(:,part)/(gammat);
            lifetime(part,1) = lifetime(part,1)+ dt; 
        end

        for looper=1:numin
            force_cell_tracer(:,step,looper)=forces(:,looper);
        end

         
      
        death = 9.728*10^(-7);
        deadpart=0;
	    deadnumin=0;
        
        for part=1:nPart
            
            if (rand <= death*dt) && (part > numin0/2)
                
                deadpart=deadpart+1;
                deadindex(deadpart,1) = part;
                
            elseif (rad(part,1) < radmitosis) && ...
                    (pressure(part,1) < p_c) && (part > numin0/2)
                
                grate = (volrate/(4*pi*rad(part,1)*rad(part,1)));

                rad(part,1) = rad(part,1) + dt*randgaussrad(grate,10^-5);
                
                
            elseif (rad(part,1) >= radmitosis) && ...
                    (pressure(part,1) < p_c) && (part > numin0/2)
                
                rad(end+1,1)=(2^(-1/3))*rad(part,1);
                
                rad(part,1) = (2^(-1/3))*rad(part,1); 
                modulus(end+1,1) = randgaussrad(10^-3,10^-4);
                poisson(end+1,1) = randgaussrad(0.5,0.02);
                 
    
                receptor(end+1,1) = randgaussrad(1.0,0.02);
                ligand(end+1,1) = randgaussrad(1.0,0.02);
                

                 lifetime(end+1,1) = 0; 
                no_of_cells_ever_involved= no_of_cells_ever_involved +1;
                label(end+1,1)= no_of_cells_ever_involved; 
                
              
                vels(:,end+1)= ([0,0,0]');
              
                
                
               
                a=0;
                b=1;
                r3=(b-a).*rand(1) + a;
                r4=pi*(b-a).*rand(1) + pi*a;
                r5=2*pi*(b-a).*rand(1) + 2*pi*a;
                psi=size(coords,2)+1;
                
                coords(1,psi) = coords(1,part)+radmitosis*(1-2^(-1/3))*sin(r4)*cos(r5);
                coords(2,psi) = coords(2,part)+radmitosis*(1-2^(-1/3))*sin(r4)*sin(r5);
                coords(3,psi) = coords(3,part)+radmitosis*(1-2^(-1/3))*cos(r4);
                
                coords(1,part)=coords(1,part)-radmitosis*(1-2^(-1/3))*sin(r4)*cos(r5);
                coords(2,part)=coords(2,part)-radmitosis*(1-2^(-1/3))*sin(r4)*sin(r5);
                coords(3,part) = coords(3,part)-radmitosis*(1-2^(-1/3))*cos(r4);
                
            end
            
            
            
        end

        if deadpart > 0 
            
                coords(:,deadindex(1:deadpart))=[];
                rad(deadindex(1:deadpart))=[];
                modulus(deadindex(1:deadpart)) = [];
                poisson(deadindex(1:deadpart)) = [];
                receptor(deadindex(1:deadpart)) =[];
                ligand(deadindex(1:deadpart)) = [];
                vels(:,deadindex(1:deadpart))= [];
                 lifetime(deadindex(1:deadpart))=[];
                label(deadindex(1:deadpart))= [];        
	
	        for part=1:deadpart
		    
		     if (deadindex(part,1) > numin0/2) && ...
                 (deadindex(part,1)<= numin) 
                        deadnumin=deadnumin+1;
                        deadindexn(deadnumin,1)=deadindex(part,1);
             end
            
            end
        
        end
        
        if deadnumin > 0
           
                   initial(:,deadindexn(1:deadnumin,1))=[];
                   force_cell_tracer(:,:,deadindexn(1:deadnumin,1))=[];
                   track(:,:,deadindexn(1:deadnumin,1))=[];            
        end

        numin=size(initial,2);

        
        nPart=size(coords,2);
        
       
        
        for part=1:nPart
            centerM(1,1) = centerM(1,1) + coords(1,part)/nPart;
            centerM(2,1) = centerM(2,1) + coords(2,part)/nPart;
            centerM(3,1) = centerM(3,1) + coords(3,part)/nPart;
        end
        
        for part=1:nPart
            radG(step,cyclemsd) = radG(step,cyclemsd) + ...
                ((norm(coords(:,part) - centerM(:,1)))^2)/nPart;  
        end
        
        den(step,cyclemsd) =nPart/((1/3)*4*pi*radG(step,cyclemsd)^(3/2));  
        
        radG_av(step) = radG_av(step) + (radG(step,cyclemsd)/cyclenum);
        
        den_av(step) = den_av(step) + (den(step,cyclemsd)/cyclenum);
        
        
        time = time + dt;
        
        if mod(step,printFreq) == 0
            step 
            cyclemsd
            
        end
        
%         if mod(step,plotFreq) == 0
%             
%             numcells(samplecounter,cyclemsd) = nPart;
%             samplecounter = samplecounter + 1;
% 
%             
%             
%         end
        
        
        for avini=((numin0/2)+1):numin
        ds_msd=(norm(coords(:,avini)-initial(:,avini)))^2;  
        ds_msdav(step) = ds_msdav(step) + (ds_msd)/(cyclenum*(numin-(numin0/2)));
        
%         ws=norm(coords(:,avini)-initial(:,avini));
% 
%                 if ws < 1.335 
%                    qt(step,cyclemsd)=qt(step,cyclemsd)+1/(numin-(numin0/2));
%                   
%                 else
%                    qt(step,cyclemsd)=qt(step,cyclemsd)+0;
%                 end
        
        df(step,cyclemsd)=(norm(coords(:,avini)-initial(:,avini)))^4;
        df_av(step)=df_av(step)+(df(step,cyclemsd))/(cyclenum*(numin-(numin0/2)));
        end
        
        for part=1:(numin0/2)
     tds_msd=(norm(coords(:,part)-initial(:,part)))^2;  
     tds_msdav(step) = tds_msdav(step) + (tds_msd)/(cyclenum*(numin0/2));
        
%             tws=norm(coords(:,part)-initial(:,part));
% 
%               if tws < 1.335
%                  tqt(step,cyclemsd)=tqt(step,cyclemsd)+1/(numin0/2);
%                    
%               else
%                  tqt(step,cyclemsd)=tqt(step,cyclemsd)+0;
%               end

       end
        
        
       
        save('fnumin.txt','numin','-ascii','-append');
        
        
    
    
        
        
     for avini=1:numin
         
        track(:,step+1,avini)=coords(:,avini);  
        
     end
        
                

	     if mod(step,10) == 0
        	%	countdr = countdr+1;
        	%	deltart(countdr,cyclemsd) = rtumor(coords,centerM);
              		
			    for part=1:nPart
      			visdata(1,1:3)=coords(1:3,part);
      			visdata(1,4)= label(part,1);
       			visdata(1,5)= lifetime(part,1);
      			visdata(1,6)= step*dt;
      			visdata(1,7)= receptor(part,1); 
                        visdata(1,8)= ligand(part,1);  
                        visdata(1,9)= modulus(part,1);
                        visdata(1,10)= poisson(part,1);
                        visdata(1,11)= rad(part,1);
                        visdata(1,12:14)=vels(1:3,part);
                        visdata(1,15)=cyclemsd;
                        visdata_row= visdata(1,:);
       			save('lifetime1.txt','visdata_row', '-ascii','-append'); 
    			end
          end	 

%             save('lifetimex.txt','coords(1,:)','-ascii','-append');     
%             save('lifetimey.txt','coords(2,:)','-ascii','-append');     
%             save('lifetimez.txt','coords(3,:)','-ascii','-append');     


    end

    
    
    
     for step = 3:nSteps
              for avini=1:numin0/2
                     displ(:,avini) = track(:,step,avini)-track(:,2,avini);
              end
              
                deltast=step-2;
              
      for avini=1:numin0/2
          
        
          
            
         
           for angles=1:26
                    angles_coords=[sin(b2(angles)*pi/180)*cos(b1(angles)*pi/180);...
                            sin(b2(angles)*pi/180)*sin(b1(angles)*pi/180);cos(b2(angles)*pi/180)];
         
                      dotp=dot(angles_coords,displ(:,avini));
             
                 
                     scat_fun(deltast,cyclemsd)=((exp(1i*dotp*2*pi/8)*b3(angles))...
                                                   /((numin0/2)))+scat_fun(deltast,cyclemsd);
             
                       scat_fun1(deltast)=((exp(1i*dotp*2*pi/8)*b3(angles))...
                                                 /(cyclenum*(numin0/2)))+scat_fun1(deltast);
                                             
              %         particle_scat=particle_scat+ (exp(1i*dotp*2*pi/8)*b3(angles));                      
            
             
           end
          
      end
      
    end
    
  
    
    
    
    % Simulation results
    % ===================
    radG_inst(1,1) = sqrt(radG(step,cyclemsd));
    rho = nPart/((4/3)*pi*radG_inst^3);
    
    L=2*radG_inst(1,1);
    dL = 1.0;
    
    
end
 
 

%   save('numcellf.mat','numcells');
%   save('fds_msd.txt','ds_msd','-ascii');
  save('fds_msdav.mat','ds_msdav');
  save('t_msd.mat','tds_msd');
  save('t_msdav.mat','tds_msdav');
  %save('fradG.txt','radG','-ascii');
  %save('fradG_av.txt','radG_av','-ascii');
 % save('fden.txt','den','-ascii');
  %save('fden_av.txt','den_av','-ascii');
  save('radius.mat','rad');
  save('scat_fun.mat','scat_fun');
  save('scat_fun1.mat','scat_fun1');
  %save('deltart.txt','deltart','-ascii'); 
  %save('qpara.txt','qt','-ascii');
  %save('tqpara.txt','tqt','-ascii');
  save('initial.mat','initial');
  save('coords.mat','coords');
  save('track.mat','track');
  save('modulus.mat','modulus');
  save('poisson.mat','poisson');
  save('receptor.mat','receptor');
  save('ligand.mat','ligand');
  save('forces.mat','forces');
  save('pressure.mat','pressure');
save('force_cell_tracer','force_cell_tracer');

                                             
   
