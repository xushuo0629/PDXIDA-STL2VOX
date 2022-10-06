clc;clear;
addpath './utilize'
%%  parameter
name = 'pipe';
% path_in = '.\data_STL\pipe-in\realistic-bendable-pipe-in sheetmetal.STL';
path_in = '.\data_STL\pipe2\Tubo Dobrado 180.STL';
path_out = ['.\data_mat\' , name,'.mat'];
pixels = 512;
B = [0,1,0;1,1,1;0,1,0];
% B = [1,1,1;1,1,1;1,1,1];
% mu = 1;
%% read STL and plot
[vertice,face,normal] = stlRead(path_in);
%-----------Surface plot ------
% figure('color','w','position',[200,100,1200,500]);
figure('color','w');
license = 'Copyright (c) 2022 Xu Shuo, INET501, Tsinghua';
sgtitle(license,'FontSize',14,'FontName','Times New Roman','FontAngle','italic');
subplot(1,2,1);
T0 = triangulation(face,vertice(:,1),vertice(:,2),vertice(:,3));
trimesh( T0,'FaceColor','r','EdgeColor','b','FaceAlpha',0.5);
axis equal
title({name,'Triangular Surface'});
set(gca, 'Units', 'normalized', 'Position', [0.10 0.05 0.35 0.90])
%-----------Patch plot ------ 
subplot(1,2,2);
stlPlot(vertice,face,{name,'Patch Render'});
set(gca, 'Units', 'normalized', 'Position', [0.50 0.05 0.35 0.90])


%% grid boundary voxel
maxN=ceil(max(vertice)+2);
minN=ceil(min(vertice)-2);
dstep=max(maxN-minN)/(pixels-1);
%------ voxel size balance-------
middle = (maxN+ minN)/2;
scale = maxN-minN;
index = find( scale == max(scale));
xrange = linspace(middle(1)-dstep*pixels/2,middle(1)+dstep*pixels/2,pixels);
yrange = linspace(middle(2)-dstep*pixels/2,middle(2)+dstep*pixels/2,pixels);
zrange = linspace(middle(3)-dstep*pixels/2,middle(3)+dstep*pixels/2,pixels);
%------ voxel mesh   -------------
% xrange = minN(1):dstep:maxN(1);  
% yrange = minN(2):dstep:maxN(2);
% zrange = minN(3):dstep:maxN(3);
img_surf=surf2vol(vertice,face(:,1:3),xrange,yrange,zrange);
img_surf = single(img_surf);
% zslice=15;
% imagesc(squeeze(img(:,:,zslice))); % z=10
%% imdilate if necessary
img_surf = imdilate(img_surf,B);
%%  imfill and noise
img = imfill(img_surf,'holes');
img   ( img~= 0) = 1;
% img   ( img~= 0) = mu;
% noise = normrnd(mu*0,mu*0.02,size(img)) ;
% img = img  + noise;
img = uint8(img);
%% show voxel and save
ph = img;
volumeViewer(ph);
% ph = permute(ph,[1 3 2]);
save(path_out,'ph')
% save('out/gunACPint.mat','img');

