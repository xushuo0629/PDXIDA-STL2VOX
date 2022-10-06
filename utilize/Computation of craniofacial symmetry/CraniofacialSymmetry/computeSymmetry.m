function [V,M,T] = computeSymmetry(Ver,Tri)
% -------------------------------------------------------------------------
%
% This function provides an implementation of the algorithm propose by
% "M Pinheiro, X Ma, MJ Fagan, GT McIntyre, P Lin, G Sivamurthy, PA Mossey 
% A 3D cephalometric protocol for the accurate quantification of the cranio-
% facial symmetry and facial growth", published in May 2019 in the Journal 
% of Biological Engineering.
%
% Input parameters:
%   Ver: a [ m x 3 ] matrix contraining the vertices of the 3D model
%   Tri: a [ n x 3 ] trinagulation of the surface
%
% Output parameters:
%   V: a [ m x 3 ] matrix contraining the vertices after transformation
%   M: a [ 3 x 3 ] rotation matrix from the initial to the final position
%   T: a [ 1 x 3 ] translation vector
% 
% The implementation uses by default the pcregrigid/pcregistericp function(s)
% and therefore it depends on the MATLAB Computer Vision System Toolbox. 
% An additional implementation using the ICP method provided by Martin Kjer 
% and Jakob Wilm from the Technical University of Denmark is provided.  
% The source code for the ICP implementation can be found in:
%
% http://mathworks.com/matlabcentral/fileexchange/27804-iterative-closest-point
%
% This code also uses the stlTools provided by Pau Micó for STL input and output.
% This toolbox can be found in:
%
% http://mathworks.com/matlabcentral/fileexchange/51200-stltools
%
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% Pre-alignment of the data
CM  = mean( Ver,1 );

% Remove any existing displacement
Ver = [ Ver(:,1)-CM(1,1), Ver(:,2)-CM(1,2), Ver(:,3)-CM(1,3) ] ;

% Computation of the PCA
[COEFF,~,~] = pca( Ver );

% Re-order the eigenvectors and check for unwanted reflections on the PCA scores
M1 = PCAcoefftoMatrix( COEFF );

V01 = zeros( size(Ver) );
for i = 1:size( V01,1 )
    V01(i,:) = ( M1*Ver(i,:)' ); 
end
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% Test the 3 symmetry planes
rmse = zeros(1,3);

% Symmetry on the Ox
V02 = [ -V01(:,1), V01(:,2), V01(:,3) ];
% Symmetry on the Oy
V03 = [ V01(:,1), -V01(:,2), V01(:,3) ];    

% Symmetry on the Oz
V04 = [ V01(:,1), V01(:,2), -V01(:,3) ];

bMCVST = 0;
try
    % Using MATLAB Computer Vision System Toolbox
    pcRef = pointCloud( V01(1:2:end,:) );
    pcMov = pointCloud( V02(1:2:end,:) );
    [~,~,rmse(1,1)] = pcregrigid( pcMov,pcRef,'Metric','pointToPoint','Extrapolate', true,'MaxIterations',10 );

    pcRef = pointCloud( V01(1:2:end,:) );
    pcMov = pointCloud( V03(1:2:end,:) );
    [~,~,rmse(1,2)] = pcregrigid( pcMov,pcRef,'Metric','pointToPoint','Extrapolate', true,'MaxIterations',10 );

    pcRef = pointCloud( V01(1:2:end,:) );
    pcMov = pointCloud( V04(1:2:end,:) );
    [~,~,rmse(1,3)] = pcregrigid( pcMov,pcRef,'Metric','pointToPoint','Extrapolate', true,'MaxIterations',10 );
    bMCVST = 1;
catch
    % Using MathWorks file exchange ICP method
    p = unique( V01(1:2:end,:),'rows' ); 
    q = unique( V02(1:2:end,:),'rows' ); 
    [~,~,err,~] = icp( p',q','Minimize','point' );
    rmse(1,1)   = err(end,1);
    
    q = unique( V03(1:2:end,:),'rows' ); 
    [~,~,err,~] = icp( p',q','Minimize','point' );
    rmse(1,2)   = err(end,1);
    
    q = unique( V04(1:2:end,:),'rows' ); 
    [~,~,err,~] = icp( p',q','Minimize','point' );
    rmse(1,3)   = err(end,1);
end
clear V02 V03 V04 
    
% Symmetric dataset
[~,ind] = min(rmse);
switch ind
    case 1, V02 = [ -V01(:,1), V01(:,2), V01(:,3) ];
    case 2, V02 = [ V01(:,1), -V01(:,2), V01(:,3) ];
    case 3, V02 = [ V01(:,1), V01(:,2), -V01(:,3) ];
end
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% Compute new normals and compute normals in the points
Nor01 = computePointNormals( V01,Tri ); 
Nor02 = computePointNormals( V02,Tri );

ntri  = size(Tri,1);
N1    = zeros( size(V01) );
N2    = zeros( size(V02) );
for j = 1:ntri
    N1( Tri(j,:),: ) = N1( Tri(j,:),: ) + repmat( Nor01(j,:),3,1 );
    N2( Tri(j,:),: ) = N2( Tri(j,:),: ) + repmat( Nor02(j,:),3,1 );
end
[counts,~] = hist( Tri(:), size(V01,1) );
N1 = ( N1./counts' );
N1 = N1./sqrt( sum( N1.^2, 2 ) );
N2 = ( N2./counts' );
N2 = N2./sqrt( sum( N2.^2, 2 ) ); 
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% Align the initial and mirrored meshes
if bMCVST == 1
    pcRef = pointCloud( V01,'Normal',N1 );
    pcMov = pointCloud( V02,'Normal',N2 );    
    [tform,~,rmse] = pcregrigid( pcMov,pcRef,'Metric','pointToPlane','Extrapolate', true );
    
    M2 = tform.T(1:3,1:3);
    T2 = tform.T(4,1:3);  
else
    p = unique( V01,'rows' ); 
    q = unique( V02,'rows' ); 
    [r,t,rmse,~] = icp(p',q','Extrapolation',true,'Minimize','plane','Normals',N1' ); 
    M2   = inv(r);
    T2   = t';
    rmse = rmse(end,1);
end
disp( [ 'RMSE : ' num2str( rmse ) ] );

for j = 1:size(V02,1),V02(j,:) = ( M2\V02(j,:)' )'; end
V02 = [ V02(:,1)+T2(1,1), V02(:,2)+T2(1,2), V02(:,3)+T2(1,3) ];   
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% Compute the joint PCA of the two models combined
[COEFF,~,~] = pca( [ V01; V02 ] );

% Check for unwanted reflections on the eigenvectors of the PCA
M3 = PCAcoefftoMatrix( COEFF );

% Compute the final model alignment
for i = 1:size( V01,1 ),V01(i,:) = ( M3*V01(i,:)' ); end
V = V01;
M = M3*M1;
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% Double-check computed transformation (Optinal)
V03 = zeros( size(Ver) );
for i = 1:size( V03,1 ), V03(i,:) = ( M*Ver(i,:)' ); end

disp('Max. diff. between iterative and combined transformation:');
disp( max( sqrt( sum( ( V01 - V03 ).^2,2 ) ) ) );

% Double-check for additional translations (Optinal)
Nor01 = computePointNormals( V01,Tri );  
Nor02 = computePointNormals( Ver,Tri );

ntri  = size(Tri,1);
N1    = zeros( size(V01) );
N2    = zeros( size(V03) );
for j = 1:ntri
    N1( Tri(j,:),: ) = N1( Tri(j,:),: ) + repmat( Nor01(j,:),3,1 );
    N2( Tri(j,:),: ) = N2( Tri(j,:),: ) + repmat( Nor02(j,:),3,1 );
end
[counts,~] = hist( Tri(:), size(V01,1) );
N1 = ( N1./counts' );
N1 = N1./sqrt( sum( N1.^2, 2 ) );
N2 = ( N2./counts' );
N2 = N2./sqrt( sum( N2.^2, 2 ) ); 

if bMCVST == 1
    pcRef = pointCloud( V01,'Normal',N1 );
    pcMov = pointCloud( Ver,'Normal',N2 );    
    [tform,~,rmse] = pcregrigid( pcMov,pcRef,'Metric','pointToPlane','Extrapolate', true );
    tM = inv( tform.T( 1:3,1:3) );
    T  = tform.T(4,1:3);
else
    % NOTE: ignore possible warnings - we are double checking the solution
    [r,t,rmse,~] = icp(V01',Ver','Extrapolation',true,'Minimize','plane','Normals',N1' );
    tM   = r;
    T    = t';
    rmse = rmse(end,1);
end
disp( [ 'RMSE : ' num2str( rmse ) ] );

if max( max( M - tM ) ) > 1e-10
    disp( 'Initial transformation:' ); disp(M);
    M = tM;
    disp( 'Final transformation:'); disp(tM);
end
% -------------------------------------------------------------------------
end

function pn = computePointNormals(v,t)
% Compute new normals and compute normals in the points
facets = single(v');
facets = reshape(facets(:,t'), 3, 3, []);
v1     = squeeze(facets(:,2,:) - facets(:,1,:));
v2     = squeeze(facets(:,3,:) - facets(:,1,:));
nor    = v1([2 3 1],:).*v2([3 1 2],:) - v2([2 3 1],:).*v1([3 1 2],:);
pn     = bsxfun( @times,nor,1./sqrt( sum( nor.*nor, 1) ) )';

end

function rm = PCAcoefftoMatrix(coeff)
% Convert the PCA matrix into a rotation matrix
[~,ind] = max( abs(coeff), [], 1 );
if size(unique(ind),2) < 3
    if     isempty( find( ind==1,1 ) ),ind(1,3) = 1;
    elseif isempty( find( ind==2,1 ) ),ind(1,3) = 2;
    else,  isempty( find( ind==3,1 ) ),ind(1,3) = 3;
    end
end
rm = coeff( :,ind );

if det(rm) > 0,rm = inv( rm );
else, error('Matrix not invertible.');
end

end