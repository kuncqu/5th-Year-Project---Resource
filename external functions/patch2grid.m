function [Z,M] = patch2grid(varargin)
% PATCH2GRID  Convert patch data to gridded data.
%   Z = patch2grid(P,xgv,ygv) converts the patch data from PATCH object (or
%   handle) P into gridded data Z. Each element in Z is calculated by
%   taking the Cdata of the patch which occupies the largest area on a
%   given grid point. xgv and ygv are grid vectors and the size of Z is
%   numel(ygv)-by-numel(xgv). Note that the grid point is rectangular in
%   shape and its dimensions are (xgv(2)-xgv(1))-by-(ygv(2)-ygv(1)).
% 
%   Z = patch2grid(X,Y,C,xgv,ygv) calculates gridded data from X, Y, and C
%   patch data. X and Y are m-by-n matrices with m vertices and n polygons.
%   C is the color data of the patch.
% 
%   Z = patch2grid(fvc,xgv,ygv) calculates gridded data using structure
%   fvc, which contains the patch properties fields vertices, faces, and
%   facevertexcdata.
% 
%   Z = patch2grid(...,'weighted') uses an area weighted method to compute
%   Z using any of the input arguments in the previous syntaxes. The area
%   weighted method computes each element in Z by weighing the Cdata value
%   of each patch by the total area that they occupy on a given patch and
%   then calculating the weighted average.
% 
%   [Z, M] = patch2grid(...) returns the gridded data Z and the cell array
%   M of size numel(ygv)-by-numel(xgv) which lists the fractional area
%   occupied by each patch on a specific grid point. Fractional area is
%   defined as the intersecting area of a patch with a grid point divided
%   by the total area of the grid point. The numeric array at M{i,j} is of
%   size 2-by-p where p is the number of patches inside the grid point at
%   the index coordinate (i,j), the 1st row is the fractional area value,
%   and the 2nd row is the corresponding patch index.
% 
%   Example 1:
%      X = [2 2 0 2 5;...
%           2 8 2 4 5;...
%           8 8 2 4 8];
%      Y = [4 4 4 2 0;...
%           8 4 6 2 2;...
%           4 0 4 0 0];
%      C = [15 30 25 2 60]';
%      
%      patch(X,Y,C');
%      axis off
% 
%      xgv=0:8;
%      ygv=xgv;
%      [Z,M]=patch2grid(X,Y,C,xgv,ygv);
% 
%      figure, pcolor(xgv,ygv,Z)
%      axis off
% 
%   Example 2:
%      vertices = [0 0; 0 .5; 4/9 1; 1 1; 1 .5; 5/9 0];
%      faces = [1 2 3 6; 3 4 5 6];
%      facevertexcdata = [0;1];
%      fvc = struct('vertices',vertices,...
%                   'faces',faces,...
%                   'facevertexcdata',facevertexcdata);
% 
%      patch(fvc)
%      shading faceted
%      colormap summer
%      axis off
%      colorbar
% 
%      xgv=linspace(0,1,10);
%      ygv=xgv;
%      Z=patch2grid(fvc,xgv,ygv,'weighted');
% 
%      figure, pcolor(xgv,ygv,Z)
%      colormap summer
%      axis off
%      colorbar
% 
%   Requires the functions POLY2CW and POLYBOOL from the Mapping Toolbox™.
% 
%   Copyright 2015 Ivan Kostylev
%   Version: 1.0

%Determine the grid vectors
endarg = ischar(varargin{nargin}); %1 if char else 0
if nargin==3+endarg
    xgv=varargin{2};
    ygv=varargin{3};
elseif nargin==5+endarg
    xgv=varargin{4};
    ygv=varargin{5};
else
    error('Not enough input arguments.')
end

if numel(xgv)==1||numel(ygv)==1 %Error check
    error('The grid vectors xgv and ygv both need to have at least 2 elements.')
end

%Determine the type of patch data input
if isa(varargin{1},'matlab.graphics.primitive.Patch') %IF patch class object
    fvc.vertices = varargin{1}.Vertices;
    fvc.faces = varargin{1}.Faces;
    C = varargin{1}.FaceVertexCData;
elseif ishandle(varargin{1}) %IF patch handle (for older version support)
    fvc.vertices = get(varargin{1},'Vertices');
    fvc.faces = get(varargin{1},'Faces');
    C = get(varargin{1},'FaceVertexCData');
elseif isa(varargin{1},'numeric')&&prod(single(size(varargin{1})==size(varargin{2}))) %IF X and Y array
    fvc = {varargin{1},varargin{2}}; %Note: fvc doesn't represent Face-Vertex here anymore
    C = varargin{3}; %C array
elseif isstruct(varargin{1}) %IF structure
    fvc = varargin{1};
    C = fvc.facevertexcdata;
else
    error('Unsupported input')
end

%Extract parameters
x0=xgv(1);
y0=ygv(1);
dX=xgv(2)-xgv(1);
dY=ygv(2)-ygv(1);
sizeX=numel(xgv);
sizeY=numel(ygv);

%Determine area contribution of each patch on each grid rectangle
M = patchareaongrid(fvc,x0,y0,dX,dY,sizeX,sizeY);

%Create 2D grid with C values
if endarg&&strcmp(varargin{nargin},'weighted')
    Z = cellfun(@(x) WeightedAverageCValue(x,C),M);
else
    Z = cellfun(@(x) MaxAreaCValue(x,C),M);
end
end

%-----------------------------------------------------------------

function M = patchareaongrid(fvc,x0,y0,dX,dY,sizeX,sizeY)
%Determine Input Type
if iscell(fvc)
    getxy = @(x,y) cell2xy(x,y) ; %Function to extract polygon vertices from cell array
    NoP = size(fvc{1},2); %Number of polygons
else %structure
    getxy = @(x,y) fv2xy(x,y); %Function to extract polygon vertices from structure
    NoP = size(fvc.faces,1); %Number of polygons
end

A2=dX*dY; %area of a single grid rectangle
M=cell(sizeY,sizeX); %preallocate
for p=1:NoP %Loop over all the polygons
    [x,y]=getxy(fvc,p); %extract x and y vectors
    %upper and lower bounds on affected grid points by patch
    jmin=floor((min(x)-x0)/dX)+1;
    jmax=floor((max(x)-x0)/dX)+1;
    kmin=floor((min(y)-y0)/dY)+1;
    kmax=floor((max(y)-y0)/dY)+1;
    for j=jmin:jmax %for all horizontal grid points
        if (j>=1&&j<=sizeX) %If within bounds of the sampling grid
            x2=([j j j+1 j+1 j]-1)*dX+x0;
            for k=kmin:kmax %for all vertical grid points
                if (k>=1&&k<=sizeY) %If within bounds of the sampling grid
                    y2=([k k+1 k+1 k k]-1)*dY+y0;
                    [xb, yb] = polybool('intersection', x, y, x2, y2); %Polygon defining the intersection between the grid point and the patch
                    A1=polyarea(xb, yb); %Area of the intersecting polygon
                    if A1~=0 %A1 returns 0 if xb,yb are empty so this value isn't recorded
                        M{k,j}=[M{k,j} [A1/A2;p]]; %A1/A2 was originally used
                    end
                end
            end
        end
    end
end
end

%-----------------------------------------------------------------

function [x,y]=fv2xy(fvc,p)
% Convert patch face and vertices to x and y polygonal region
xy=fvc.vertices(fvc.faces(p,:),:);
x=xy(:,1)';
y=xy(:,2)';
if ~ispolycw(x, y) %Make CW instead of CCW for polybool
    [x, y] = poly2cw(x, y);
end
end

%-----------------------------------------------------------------

function [x,y]=cell2xy(fvc,p)
% Split XY cell array to x and y variables
x=fvc{1}(:,p)';
y=fvc{2}(:,p)';
end

%-----------------------------------------------------------------

function y=MaxAreaCValue(x,C)
% Take the C value of the patch that occupies the most area on the grid
if ~isempty(x)&&sum(x(1,:))>=0.5 %unoccupied area should be less than 50% of total area
    x(:,isnan(x(1,:)))=[]; %Remove NaNs
    if ~isempty(x)
        mx=max(x(1,:));
        ind=x(2,x(1,:)==mx);
        y=C(ind(1)); %(1) ensures it takes only the 1st index if there are multiple with the same maximum value
    else
        y=NaN;
    end
else
    y=NaN;
end
end

%-----------------------------------------------------------------

function y=WeightedAverageCValue(x,C)
% Calculate the patch area weighted average of C
if ~isempty(x)&&sum(x(1,:))>=0.5 %unoccupied area should be less than 50% of total area
    x(:,isnan(x(1,:)))=[]; %Remove NaNs
    if ~isempty(x)
        if size(C,1)==1, C=C'; end
        y=x(1,:)*C(x(2,:))./sum(x(1,:)); %Area Weighted Average
    else
        y=NaN;
    end
else
    y=NaN;
end
end