clear;close all;clc
rm=3393.5;
aj=0.4;
ej=0.4;
Az=aj:aj:360;
El=ej:ej:180;
az=Az/180*pi;
el=El/180*pi;
 maxD = 110;
 n=maxD^2+2*maxD;
[AZ,EL]=meshgrid(az,el);
al=120;
h=rm+al*ones(length(el),length(az));
phi=reshape(AZ,1,[]);
theta=reshape(EL,1,[]);
r=reshape(h,1,[])/3393.5;
[A_r_gsm, A_theta_gsm, A_phi_gsm] = SH(r, theta, phi,maxD);
load cof_test.mat;
x0=x0(1:n);
br=A_r_gsm*x0;
bt=A_theta_gsm*x0;
bp=A_phi_gsm*x0;


br=reshape(br,[],length(az));
bt=reshape(bt,[],length(az));
bp=reshape(bp,[],length(az));
br=flipud(br);
bt=flipud(bt);
bp=flipud(bp);
ball=(br.^2+bt.^2+bp.^2).^0.5;

lat=linspace(-90,90,352);
lon=linspace(-180,180,720);
[Lon,Lat]=meshgrid(lon,lat);

%%     
   ax=figure('Color','White');
 
   load geoid
   load mola128_dx
         load('BlRe.rgb')
      load('hotres.rgb')
      hotres=hotres/256;
      hotres(1:4,:)=[0.9,0.9,0.9;0.9,0.9,0.9;0.9,0.9,0.9;0.9,0.9,0.9];
   BlRe=BlRe/256;
   BlRe(48:49,:)=[0.9,0.9,0.9;0.9,0.9,0.9];
      %%
   
   pos=[1 1 20 30];
set(gcf,'unit','centimeters','position',pos)

%%
ln=1000;


ax1=axes('position',[0.11,0.67,0.18,0.18])  ;
axesm('ortho','maplatlim',[60 90]);

gridm on;
framem on;
axis off

%mlabel('equator')
%mlabel;
%mlabel('south')
%plabel('west')
%plabel('fontweight','bold')


geoidrefvec=[1/aj,90,0];  
B=geoshow(ball, geoidrefvec, 'DisplayType', 'texturemap'); hold on;
text(0.47,0.37,'120^o','Rotation',-60)
text(0.45,-0.35,'60^o','Rotation',60)
text(-0.55,-0.2,'300^o','Rotation',-60)
text(-0.57,0.22,'240^o','Rotation',60)
h1=hot(256); 
h1(239,:)=[0.9 0.9 0.9];
colormap(ax1,hotres)
caxis([0 ln]);


%%

ax1=axes('position',[0.3,0.66,0.4,0.2])  ;
axesm eckert4; 
framem;
gridm;
axis off

geoidrefvec=[1/aj,90,180];  
B=geoshow(ball, geoidrefvec, 'DisplayType', 'texturemap'); hold on; 
[c,h] = contourm(Lat,Lon,mola128_dx,-40000:20000:40000,'k','LineWidth',0.3);; hold on;

h1=hot(256); 
h1(239,:)=[0.9 0.9 0.9];
colormap(ax1,hotres)
caxis([0 ln]);


% text(-3.3,-0.5,'Latitude','FontSize',10,'rotation',90)

hcb = colorbar('location','southoutside','position',[0.25 0.85 0.5 0.015]);
hcb.Label.String = 'nT';
hcb.Label.Position = [1030, 0.1, 0]
hcb.FontSize=10;
%hcb = colorbar('position',[0.1 0.6 0.8 0.03]);

%set(get(hcb,'Xlabel'),'String','std (nT)')

text(-0.1,1.7,'|B|','FontSize',12)
text(2.1,1.2,'90^o')
text(2.75,0,'0^o')
text(2,-1.2,'-90^o')

text(-2.54,1.2,'90^o')
text(-3.0,0,'0^o')
text(-2.54,-1.2,'-90^o')

text(-1.7,-1.5,'0^o')
text(-0.1,-1.5,'180^o')
text(1.5,-1.5,'360^o')





   %%
ax1=axes('position',[0.71,0.67,0.18,0.18])  ;
axesm('ortho','maplatlim',[-90 -60]);

gridm on;
framem on;
axis off

%mlabel('equator')
%mlabel;
%mlabel('north')
%plabel;
%plabel('meridian');
%plabel('fontweight','bold')


geoidrefvec=[1/aj,90,0];  
B=geoshow(ball, geoidrefvec, 'DisplayType', 'texturemap'); hold on;
h1=hot(256); 
h1(239,:)=[0.9 0.9 0.9];
text(0.47,0.37,'60^o','Rotation',-60)
text(0.45,-0.35,'120^o','Rotation',60)
text(-0.55,-0.2,'240^o','Rotation',-60)
text(-0.57,0.22,'300^o','Rotation',60)
colormap(ax1,hotres)
caxis([0 ln]);


%%
   
L1=500;

ax2=axes('position',[0.11,0.51,0.18,0.18])  ;
axesm('ortho','maplatlim',[60 90]);

gridm on;
framem on;
axis off

geoidrefvec=[1/aj,90,0];  
B=geoshow(br, geoidrefvec, 'DisplayType', 'texturemap'); hold on;
text(0.47,0.37,'120^o','Rotation',-60)
text(0.45,-0.35,'60^o','Rotation',60)
text(-0.55,-0.2,'300^o','Rotation',-60)
text(-0.57,0.22,'240^o','Rotation',60)
colormap(ax2,BlRe)
caxis([-L1 L1]);

ax2=axes('position',[0.3,0.50,0.4,0.2])  ;
axesm eckert4; 
framem;
gridm;
axis off

geoidrefvec=[1/aj,90,180];  
B=geoshow(br, geoidrefvec, 'DisplayType', 'texturemap'); hold on;
[c,h] = contourm(Lat,Lon,mola128_dx,-40000:20000:40000,'k','LineWidth',0.3);; hold on;
colormap(ax2,BlRe)
caxis([-L1 L1]);

text(-0.1,1.7,'B_r','FontSize',12)
% text(-3.3,-0.5,'Latitude','FontSize',10,'rotation',90)

text(2.1,1.2,'90^o')
text(2.75,0,'0^o')
text(2,-1.2,'-90^o')

text(-2.54,1.2,'90^o')
text(-3.0,0,'0^o')
text(-2.54,-1.2,'-90^o')

text(-1.7,-1.5,'0^o')
text(-0.1,-1.5,'180^o')
text(1.5,-1.5,'360^o')

% 
% ax2=axes('position',[0.3,0.50,0.4,0.2])  ;
% axesm eckert4; 
% %framem;
% %gridm;
% axis off



ax2=axes('position',[0.71,0.51,0.18,0.18])  ;
axesm('ortho','maplatlim',[-90 -60]);

gridm on;
framem on;
axis off
geoidrefvec=[1/aj,90,0];  
B=geoshow(br, geoidrefvec, 'DisplayType', 'texturemap'); hold on;
text(0.47,0.37,'60^o','Rotation',-60)
text(0.45,-0.35,'120^o','Rotation',60)
text(-0.55,-0.2,'240^o','Rotation',-60)
text(-0.57,0.22,'300^o','Rotation',60)
colormap(ax2,BlRe)
caxis([-L1 L1]);

%%


ax2=axes('position',[0.11,0.35,0.18,0.18])  ;
axesm('ortho','maplatlim',[60 90]);

gridm on;
framem on;
axis off

geoidrefvec=[1/aj,90,0];  
B=geoshow(bt, geoidrefvec, 'DisplayType', 'texturemap'); hold on;
text(0.47,0.37,'120^o','Rotation',-60)
text(0.45,-0.35,'60^o','Rotation',60)
text(-0.55,-0.2,'300^o','Rotation',-60)
text(-0.57,0.22,'240^o','Rotation',60)
colormap(ax2,BlRe)
caxis([-L1 L1]);

ax2=axes('position',[0.3,0.34,0.4,0.2])  ;
axesm eckert4; 
framem;
gridm;
axis off

geoidrefvec=[1/aj,90,180];  
B=geoshow(bt, geoidrefvec, 'DisplayType', 'texturemap'); hold on;
[c,h] = contourm(Lat,Lon,mola128_dx,-40000:20000:40000,'k','LineWidth',0.3);; hold on;
colormap(ax2,BlRe)
caxis([-L1 L1]);

text(-0.1,1.7,'B_{\theta}','FontSize',12)
% text(-3.3,-0.5,'Latitude','FontSize',10,'rotation',90)

text(2.1,1.2,'90^o')
text(2.75,0,'0^o')
text(2,-1.2,'-90^o')

text(-2.54,1.2,'90^o')
text(-3.0,0,'0^o')
text(-2.54,-1.2,'-90^o')

text(-1.7,-1.5,'0^o')
text(-0.1,-1.5,'180^o')
text(1.5,-1.5,'360^o')

% 
% ax2=axes('position',[0.3,0.34,0.4,0.2])  ;
% axesm eckert4; 
% %framem;
% %gridm;
% axis off



ax2=axes('position',[0.71,0.35,0.18,0.18])  ;
axesm('ortho','maplatlim',[-90 -60]);

gridm on;
framem on;
axis off
geoidrefvec=[1/aj,90,0];  
B=geoshow(bt, geoidrefvec, 'DisplayType', 'texturemap'); hold on;
text(0.47,0.37,'60^o','Rotation',-60)
text(0.45,-0.35,'120^o','Rotation',60)
text(-0.55,-0.2,'240^o','Rotation',-60)
text(-0.57,0.22,'300^o','Rotation',60)
colormap(ax2,BlRe)
caxis([-L1 L1]);
%%

ax2=axes('position',[0.11,0.19,0.18,0.18])  ;
axesm('ortho','maplatlim',[60 90]);

gridm on;
framem on;
axis off

geoidrefvec=[1/aj,90,0];  
B=geoshow(bp, geoidrefvec, 'DisplayType', 'texturemap'); hold on;
text(0.47,0.37,'120^o','Rotation',-60)
text(0.45,-0.35,'60^o','Rotation',60)
text(-0.55,-0.2,'300^o','Rotation',-60)
text(-0.57,0.22,'240^o','Rotation',60)
colormap(ax2,BlRe)
caxis([-L1 L1]);

ax2=axes('position',[0.3,0.18,0.4,0.2])  ;
axesm eckert4; 
framem;
gridm;
axis off

geoidrefvec=[1/aj,90,180];  
B=geoshow(bp, geoidrefvec, 'DisplayType', 'texturemap'); hold on;
[c,h] = contourm(Lat,Lon,mola128_dx,-40000:20000:40000,'k','LineWidth',0.3);; hold on;
colormap(ax2,BlRe)
caxis([-L1 L1]);

text(-0.1,1.7,'B_{\phi}','FontSize',12)
% text(-3.3,-0.5,'Latitude','FontSize',10,'rotation',90)

text(2.1,1.2,'90^o')
text(2.75,0,'0^o')
text(2,-1.2,'-90^o')

text(-2.54,1.2,'90^o')
text(-3.0,0,'0^o')
text(-2.54,-1.2,'-90^o')

text(-1.7,-1.5,'0^o')
text(-0.1,-1.5,'180^o')
text(1.5,-1.5,'360^o')

hcb = colorbar('location','southoutside','position',[0.25 0.18 0.5 0.015]);
% text(0.75,0.18,'nT','FontSize',10)
% set(get(hcb,'Xlabel'),'String','nT','FontSize',10)
 hcb.Ticks=[-500 -400 -300 -200 -100 0 100 200 300 400 500];
%   hcb.TickLabels={'-500'; '-400'; '-300'; '-200'; '-100'; '0'; '100'; '200'; '300'; '400'; '500'};
hcb.FontAngle = 'normal';
hcb.Label.String = 'nT';
hcb.Label.Position = [530, 1.1, 0]
hcb.FontSize=10;
% ax2=axes('position',[0.3,0.18,0.4,0.2])  ;
% axesm eckert4; 
% %framem;
% %gridm;
% axis off
% 
% geoidrefvec=[2,90,180]; 


ax2=axes('position',[0.71,0.19,0.18,0.18])  ;
axesm('ortho','maplatlim',[-90 -60]);

gridm on;
framem on;
axis off
geoidrefvec=[1/aj,90,0];  
B=geoshow(bp, geoidrefvec, 'DisplayType', 'texturemap'); hold on;
text(0.47,0.37,'60^o','Rotation',-60)
text(0.45,-0.35,'120^o','Rotation',60)
text(-0.55,-0.2,'240^o','Rotation',-60)
text(-0.57,0.22,'300^o','Rotation',60)
colormap(ax2,BlRe)
clim([-L1 L1]);

%%
