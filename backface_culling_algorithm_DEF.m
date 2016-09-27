clear all
close all
clc
addpath('C:\Users\GBenedetti\Google Drive\TESI LM\tesi\tps\FOSTRAD\aerothermal-master\aerothermal-master\altmany-export_fig-2bad961');

                                            %% BACKFACE CULLING ALGORITHM %%

                                        %%%% Gianluca Benedetti, May 2016 %%%%

tic
%%%%% Read in the STL file for the Geometry %%%%%%%%%
%V= vertices %N=trinorms, F=triangles

[F,V,N] = stlread('C:\Users\GBenedetti\Google Drive\TESI LM\tesi\tps\FOSTRAD\aerothermal-master\aerothermal-master\CAD MODELS\GOCE_simplified_FOSTRAD_veryfine.stl');
%[V, F, N] = import_stl_fast('C:\Users\GBenedetti\Google Drive\TESI LM\tesi\tps\FOSTRAD\aerothermal-master\aerothermal-master\CAD MODELS\GOCE_right_meas.stl',1);
b = 5;            %Characteristic Length of Object

%%%%% Truncating for unwanted rows %%%%%%%
if length(N) > length(F) 
    N = N(1:end-1,:);    
end

%%%%%% Calculating Elemental Areas %%%%%%%%%
A = Area(F,V);

%Use triangulation to create an in-memory representation of any 2-D or 3-D triangulation data that is in matrix format, 
%such as the matrix output from the delaunay function or other software tools. 
%When your data is represented using triangulation, you can perform topological and geometric queries, 
%which you can use to develop geometric algorithms. For example, you can find the triangles or tetrahedra attached to a vertex, 
%those that share an edge, their circumcenters, and other features.
tr = triangulation(F,V);

P = incenter(tr);
CG = COG(P);

r=zeros(length(N),3);

for i = 1:length(N)
r(i,:)= P(i,:) - CG;   %Moment Arm
i=i+1;
end

%%%%%% Cleaning up duplicates %%%%%%%%
[ N, F, A ] = cleanup( N, F, A );

%%%%% Modified Newtonian Theory Pressure Calculations %%%%%%%%

alpha = 10^-10;     % Angle of attack, degrees
beta = -5;       % Yaw angle, degrees
Vinf = 300*2;     % Free Stream Velocity, m/s
Vinf = [Vinf * cosd(alpha) * cosd(beta), - Vinf * sind(beta), Vinf * sind(alpha) * cosd(beta)]; % %Free Stream Velocity Vector
%%

directionX = [1 0 0];
directionY = [0 1 0];
directionZ = [0 0 1];

figure(1)
p=patch('Faces',F,'Vertices',V,'FaceVertexCdata',F,'FaceColor','flat');
rotate(p,directionZ,beta)
rotate(p,directionY,alpha)
axis off
view([-90,0,0])


% F_old=F;
% V_old=V;



% %% selezione triangoli WINDWARD- LEEWARD
% 
% LEEWARD=0;
% WINDWARD=0;
% % W=zeros(length(N),1);
% 
% a=zeros(1,length(N));
% % W=zeros(length(N),1);
% 
% for i=1:length(N)
% a(i)=dot((Vinf/norm(Vinf)),N(i,:));
% 
% if a(i)<=0
%     LEEWARD=LEEWARD+1;
%     W(i)=1;
% else
%     WINDWARD=WINDWARD+1;
%     W(i)=0;
% end
% 
% end
% W=W';
% 
% patch %windward(RED)-leeward(BLUE) triangles
% cmap = hsv(length(F));
% 
% figure(2)
% 
% subplot(2,2,1)
% p=patch('Faces',F,'Vertices',V,'FaceVertexCdata',W,'FaceColor','flat');
% hold on
% rotate(p,directionZ,beta)
% rotate(p,directionY,alpha)
% axis equal
% axis off
% view([-90,0,0])
% hold off
% 
% subplot(2,2,2)
% p=patch('Faces',F,'Vertices',V,'FaceVertexCdata',W,'FaceColor','flat');
% hold on
% rotate(p,directionZ,beta)
% rotate(p,directionY,alpha)
% axis equal
% axis off
% view([0,90,0])
% hold off
% 
% subplot(2,2,3)
% p=patch('Faces',F,'Vertices',V,'FaceVertexCdata',W,'FaceColor','flat');
% hold on
% rotate(p,directionZ,beta)
% rotate(p,directionY,alpha)
% axis equal
% axis off
% view([0,-90,0])
% hold off
% 
% subplot(2,2,4)
% p=patch('Faces',F,'Vertices',V,'FaceVertexCdata',W,'FaceColor','flat');
% hold on
% rotate(p,directionZ,beta)
% rotate(p,directionY,alpha)
% axis equal
% axis off
% view([90,0,0])
% hold off
% suptitle('STS Orbiter, AoA=30 deg, windward(RED)-leeward(BLUE) ')
% print(figure(1),'windwar_leeward,','-dmeta')
% 
% for j=1:length(N)
% if W(j)==0;
% N(j,:)= 0;
% F(j,:)= 0;
% A(j,:)= 0;
% end
% end
% 
% N( ~any(N,2), : ) = [];
% F( ~any(F,2), : ) = [];
% A( ~any(A,2), : ) = [];
% W( ~any(W,2), : ) = [];
% 

%% Back-face culling 

for i=1:1:length(F)
s = randomRGB();     %generazione random colore in EXA
G(i,:)=hex2rgb(s);   %generazione matrice per assegnazioni colori random per ogni triangolo
end


% G=unique(randi([0 255],length(F),3),'rows'); %generazione random colore
% 
% G=uint8(G);         %conversione G in uint8

myColorMap = G;

myColorMap1=im2double(myColorMap);  %conversione G in double [0 255] ----> [0 1]

colornumber(:,1:3)=[myColorMap(:,1).*256^2 myColorMap(:,2).*256 myColorMap(:,3)];

ID_COL_OR=colornumber(:,1)+colornumber(:,2)+colornumber(:,3); %ID original colors

%generate patch

figure(3)
p=patch('Faces',F,'Vertices',V,'FaceVertexCdata',myColorMap./255,'FaceColor','flat','LineStyle','none');
hold on
rotate(p,directionZ,beta)
rotate(p,directionY,alpha)
axis off
view([-90,0,0])
hold off
print('Front_Shuttle','-dpng','-r0') %Salva immagine %for a lower
% accurate analysis, low computational time
%imageData = export_fig('Front_Shuttle','-png'); %higher accurate analysis, high comp. time

CDATA_OR=get(p,'CData'); %Patch color data
CDATA_OR=im2uint8(CDATA_OR);

%patch color channels
CDATA_OR_R=CDATA_OR(:,:,1);
CDATA_OR_B=CDATA_OR(:,:,2);
CDATA_OR_G=CDATA_OR(:,:,3);

CDATA_OR=[CDATA_OR_R CDATA_OR_B CDATA_OR_G ]; %Patch color data [0 255]

CDATA_OR32=uint32(CDATA_OR);

CDATA_OR_NUMBER(:,1:3)=[CDATA_OR32(:,1).*256^2 CDATA_OR32(:,2).*256 CDATA_OR32(:,3)]; %univoco 

ID_COL_DATA_OR=(CDATA_OR_NUMBER(:,1)+CDATA_OR_NUMBER(:,2)+CDATA_OR_NUMBER(:,3));

%read image from folder

AA= imread('C:\Users\GBenedetti\Google Drive\TESI LM\tesi\tps\FOSTRAD\aerothermal-master\aerothermal-master\Front_Shuttle.png');   
figure(4)
O=image(AA,'CDataMapping','direct'); 
B1=im2double(AA); 

CDATA_FIN=get(O,'CData');       %Image color data
CDATA_FIN32=uint32(CDATA_FIN);  %Image color data [0 255]

clear i
clear k

CDATA_FIN_NUMBER=zeros(size(AA,1),size(AA,2));

    for i=1:size(CDATA_FIN32,1)%rassegna righe px
        for k=1:size(CDATA_FIN32,2) %rassegna colonnne px
            CDATA_FIN_R=CDATA_FIN32(i,k,1)*256^2;
            CDATA_FIN_G=CDATA_FIN32(i,k,2)*256;
            CDATA_FIN_B=CDATA_FIN32(i,k,3);    
            CDATA_FIN_NUMBER(i,k)=CDATA_FIN_R+CDATA_FIN_B+CDATA_FIN_G;
        end
    end
    
%eliminare elementi uguali CDATA_FIN_NUMBER

a = unique(CDATA_FIN_NUMBER);
out = [a,histc(CDATA_FIN_NUMBER(:),a)]; % mi dà i numeri asssociati ai colori che compaiono nell'immagine e la loro frequenza (n pixel con quel colore)
CDATA_FIN_NUMBER_UNICI=out(:,1);
CDATA_FIN_NUMBER_UNICI(end, :) = [];    % rimozione colore bianco (sfondo)

clear i
clear y

% creazione matrici per confronto liste
MATR_COLOR=kron(double(ID_COL_DATA_OR),ones(1,length(CDATA_FIN_NUMBER_UNICI)));


CDATA_FIN_NUMBER_UNICI_tr=CDATA_FIN_NUMBER_UNICI';
for i=1:length(F)   
   MATR_COLOR1(i,:)=CDATA_FIN_NUMBER_UNICI_tr;
end



%comparazione lista colori iniziali con colori matrice

check=MATR_COLOR-double(MATR_COLOR1);
Ch=find(~check);  %cerca elementi di check uguali a zero. Quelli elementi saranno i colori 
                    %utilizzati nei triangolini frontali della mesh di partenza

[row,col]=find(~check);      %row contiene i numeri dei triangolini in vista!!!!       

VIEW_TRIANGLES=row; %triangolini in vista da usare per il calcolo
                    

%% Selezione triangolini necessari all'analisi

num_tri=linspace(1,length(F),length(F))';
Fmod=[F num_tri];

INTER_TRIANGLES=intersect(Fmod(:,4),row);

% NON_INT_TRIANGLES=setdiff(Fmod(:,4),row);
NON_INT_TRIANGLES=setxor(Fmod(:,4),row);

%re-define F,N,A
F_new=F;
N_new=N;
A_new=A;

F_new([NON_INT_TRIANGLES(:)], :) = [];
N_new([NON_INT_TRIANGLES(:)], :) = [];
A_new([NON_INT_TRIANGLES(:)], :) = [];

% F_new([INTER_TRIANGLES(:)], :) = [];
% N_new([INTER_TRIANGLES(:)], :) = [];
% A_new([INTER_TRIANGLES(:)], :) = [];

%creazione patch con triangoli da utilizzare per l'analisi
figure(5)
p=patch('Faces',F_new,'Vertices',V,'FaceVertexCdata',F_new,'FaceColor','flat');
rotate(p,directionZ,beta)
rotate(p,directionY,alpha)
axis    off
view([-90,0,0])



figure(6)
subplot(1,2,2)
p=patch('Faces',F_new,'Vertices',V,'FaceVertexCdata',F_new,'FaceColor','flat');
rotate(p,directionZ,beta)
rotate(p,directionY,alpha)
axis    equal
view([-90,-60,30])

subplot(1,2,1)
p=patch('Faces',F_old,'Vertices',V_old,'FaceVertexCdata',F_old,'FaceColor','flat');
rotate(p,directionZ,beta)
rotate(p,directionY,alpha)
axis equal
view([-90,-60,30])


toc


