load 'coords.mat';
load 'track.mat';

%calculate the boundary of the tumor
%--------------------------------------
 deltar=0;
  
  centerM=zeros(3,1);
  centerM(1,1)=mean(coords(1,:));
  centerM(2,1)=mean(coords(2,:));
  centerM(3,1)=mean(coords(3,:));
  
  xpos(:,1)=coords(1,:);
  ypos(:,1)=coords(2,:);
  zpos(:,1)=coords(3,:);

 DT=delaunayTriangulation(xpos,ypos,zpos);

allF = [ DT(:,1) DT(:,2) DT(:,3); ...
    DT(:,1) DT(:,3) DT(:,4); ...
    DT(:,1) DT(:,4) DT(:,2); ...
    DT(:,2) DT(:,4) DT(:,3)];
 
 sortedF = sort(allF,2);

[u,m,n] = unique(sortedF,'rows');

counts = accumarray(n(:), 1);

 sorted_exteriorF = u(counts == 1,:);
 
 F = allF(ismember(sortedF,sorted_exteriorF,'rows'),:);
 
b = unique(F(:));

 for i=1:size(b,1)
   k=b(i,1);
    boundary(i,1)=xpos(k,1);
      boundary(i,2)=ypos(k,1);
        boundary(i,3)=zpos(k,1);
 end
 
 nBound=size(boundary,1);
 
 for k=1:nBound

       deltar= deltar + sqrt((boundary(k,1)-centerM(1,1))^2 + ...
       (boundary(k,2)-centerM(2,1))^2 + (boundary(k,3)-centerM(3,1))^2)/nBound;
   
 end
 %---------------------------------------------------------------------------
 %distance of initial cells from the center of mass
 %-----------------------------------------------------------------------
 for i=1:numin
     dist_init(i)=norm(coords(:,i)-centerM(:,1));
 end
 %for boundary cells
 %------------------------------------------------------
 index_b=find(dist_init>(deltar-2*mean(rad))); % indices of boundary cells
 
 no_of_boundary_cells=size(index_b,2)
 
 
 track_b=zeros(3,nSteps,size(index_b,2));
 for i=1:size(index_b,2)
     track_b(:,:,i)=track(:,2:end,index_b(i));
 end
 msd_b=zeros(nSteps,1);
 for i=1:nSteps
     for j=1:size(index_b,2)
         msd_b(i)=msd_b(i)+(norm(track_b(:,i,j)-track_b(:,1,j)))^2/(size(index_b,2));
     end
 end
 figure(1)
 loglog(msd_b,'r-o')
 %-----------------------------------------------------------
 %for interior cells
 %---------------------------------------------------------
     index_i=find(dist_init<0.4*deltar); % indices of boundary cells
     
      no_of_interior_cells=size(index_i,2)
 
 
 track_i=zeros(3,nSteps,size(index_i,2));
 for i=1:size(index_i,2)
     track_i(:,:,i)=track(:,2:end,index_i(i));
 end
 msd_i=zeros(nSteps,1);
 for i=1:nSteps
     for j=1:size(index_i,2)
         msd_i(i)=msd_i(i)+(norm(track_i(:,i,j)-track_i(:,1,j)))^2/(size(index_i,2));
     end
 end
 hold on
 loglog(msd_i,'b-o')
 %---------------------------------------------------------
 
 