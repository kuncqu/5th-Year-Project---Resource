%PATCH_RASTERIZE Rasterization of a 2D triangulation
%   PATCH_RASTERIZE rasterize a 2D triangulation into a image. The function 
%   returns an image rasterizing the triangulation. This image is similar 
%   to the one used for OpenGL's glSelectBuffer for fast "picking".
% 
%   I = PATCH_RASTERIZE(p) creates an image of size=max(p.vertices). If the 
%   f-th triangle p.faces(f,:) contains the pixels at position (i,j) then
%   I(i,j)=f.
%
%   Once I has been computed, retrieving a the index of the triangle 
%   containing the point (x,y) can be achieved by f = I( round(x), round(y) )
%
%   The function run without parameters executes Example 1. 
%   
%   Example 1:
%      clc, clear, close all;
%      % generate some data
%      P = gallery('uniformdata',[100 2],0);
%      DT = delaunayTriangulation(P);
%      p.faces = DT.ConnectivityList;
%      p.vertices = DT.Points * 300;
%      % rasterize it
%      I = patch_rasterize(p);
%      figure, hold on;
%      imagesc(I); axis image;
%      triplot(p.faces, p.vertices(:,1), p.vertices(:,2), 'color', 'white');
%      
%   See also PATCH, TRIANGULATION, 

%   Copyright 2013 Andrea.Tagliasacchi@epfl.ch
%   $Revision: 1.0 $  $Date: 2013/05/14 16:27 $
function [ I ] = patch_rasterize( p )

% with no arguments run the demo
run_demo = (nargin==0 || (nargin==1 && ischar(p) && strcmp(p,'demo')));

% build input for demo
if run_demo
    clc, clear, close all;
    run_demo = true;
    disp('patch_rasterize demo')
    P = gallery('uniformdata',[100 2],0);
    DT = delaunayTriangulation(P);
    p.faces = DT.ConnectivityList;
    p.vertices = DT.Points * 300;
end

% aliases (MATLAB is smart)
vertices = p.vertices;
faces = p.faces;

% Build triangulation data structure (for barycentric)
TR = triangulation(faces,vertices);

% Top-right vertex gives image size
I_size = fliplr( ceil( max(vertices) ) );

% outside of triangulation you have NAN
I = nan( I_size );

% for every face
for ti=1:size(faces,1)
    face = faces(ti,:);
    
    % get this triangle's bounding box   
    vs = vertices(face,:);
    pmin = floor( min(vs) );
    pmax = ceil( max(vs) );
 
    % implicit coordinates for triangle
    [XX,YY] = meshgrid( pmin(1):pmax(1), pmin(2):pmax(2) );
    XXYY = [XX(:),YY(:)];
    TI = ti*ones(size(XXYY,1),1);
    
    % get barycentric & check interior
    bary = TR.cartesianToBarycentric(TI,XXYY);
    inside = ( all( bary>=0, 2 ) & all( bary<=1, 2 ) );
    
    % set face index inside the image
    idxs = sub2ind(I_size',XXYY(inside,2),XXYY(inside,1));
    I(idxs) = ti;
end


if run_demo
	figure, hold on;
    imagesc(I); % set(h,'ydir','normal');
    triplot(p.faces, ...
            p.vertices(:,1), p.vertices(:,2), ...
            'color', 'white');
    axis image
    I = [];
end

end

