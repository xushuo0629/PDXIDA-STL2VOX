path = pwd;

% Add dependencies
addpath( [ path,'/ICP/'] );
addpath( [ path,'/stlTools/'] );
addpath( [ path,'/Example_Data/'] );

[Ver,Tri,Nor,~] = stlRead( 'Example.stl' );
[V,M,T] = computeSymmetry(Ver,Tri);

red = [1 0 0];
green = [0 1 0];
figure,hold on,...
trimesh( Tri,Ver(:,1),Ver(:,2),Ver(:,3),'FaceColor',red,'EdgeColor',(red-0.2).*((red-0.2)>0) ),...
trimesh( Tri,V(:,1),V(:,2),V(:,3),'FaceColor',green,'EdgeColor',(green-0.2).*((green-0.2)>0)),hold off;

stlWrite( [ path '\Example_aligned.stl'],Tri,V );