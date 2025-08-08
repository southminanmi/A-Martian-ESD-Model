clear;close all;clc
[r,theta,phi]=generatemesh(3373.5/3393.5,15000);
load('C:\Users\41047\mfinal_10.txt');
load mola128_dx
br=mfinal_10(1:15000);
bt=mfinal_10(15001:30000);
bp=mfinal_10(30001:45000);
m=sqrt(br.^2+bt.^2+bp.^2)*3393.5^3/3.81e7;
br=br*3393.5^3/3.81e7;
F = scatteredInterpolant(theta,phi,br) ;
F1=scatteredInterpolant(theta,phi,m) ;
F.Method = 'natural';
aj=0.4;
rm=3393.5;
ej=0.4;
Az=aj:aj:360;
El=ej:ej:180;
az=Az/180*pi-pi;
el=El/180*pi;
[xq,yq] = meshgrid(el,az);
pbr=F({el,az});
pbr=circshift(pbr,450,2);
pbr=flipud(pbr);
pm=F1({el,az});
pm=circshift(pm,450,2);
pm=flipud(pm);
lat=linspace(-90,90,352);
lon=linspace(-180,180,720);
[Lon,Lat]=meshgrid(lon,lat);
      %%
  load('hotres.rgb')
  hotres=hotres/256;
  hotres(1:4,:)=[0.9,0.9,0.9;0.9,0.9,0.9;0.9,0.9,0.9;0.9,0.9,0.9];
  load('BlRe.rgb')
   BlRe=BlRe/256;
   BlRe(48:49,:)=[0.9,0.9,0.9;0.9,0.9,0.9];
a=figure('Color','White');
   pos=[1 1 20 20];
set(gcf,'unit','centimeters','position',pos)
%subplot(2,1,1)
ax1=axes('position',[0.13,0.67,0.4,0.2])  ; 
axesm eckert4; 
framem on;
gridm on;
axis off;
geoidrefvec=[1/aj,90,180];  
[c,h] = contourm(Lat,Lon,mola128_dx,-40000:20000:40000,'k'); hold on;
% text(2.5,0.9,'45^o')
% text(2.75,0,'0^o')
% text(2.4,-0.9,'-45^o')

text(-2.9,0.9,'45^o')
text(-3.1,0,'0^o')
text(-3.0,-0.9,'-45^o')

text(-1.7,-1.5,'0^o')
text(-0.3,-1.5,'180^o')
text(1.5,-1.5,'360^o')
%load mycolormap4.mat
%colormap(mycolormap4)
cbar = colorbar('location','southoutside','position',[0.23 0.635 0.2 0.02]);
cbar.FontSize=10;
ax=gca;
colormap(ax,BlRe)
%colormap(gca,mycolormap4)
cbar.Label.String = 'M_r(A/m)';%30/48=0.625
cbar.Label.Position = [20.5, 1.1, 0]
cbar.Ticks=[-15 0 15];
cbar.TickLabels=[{'-15' ;'0' ;'15' }];
% 调整 colorbar 的宽度
clim([-15 15]);
ax1=axes('position',[0.53,0.67,0.4,0.2])  ;
axesm eckert4; 
framem;
gridm;
axis off
i=0:255;
i=i';
cmap=[i/256,i/256,i/256];
axis off

% geoidrefvec=[1024/360,90,180]; 
% topo=flipud([test(:,513:1024),test(:,1:512)]);
% A=geoshow(topo, geoidrefvec, 'DisplayType', 'texturemap'); hold on;
% colormap(ax1,cmap)
% clim([60,200])
hold on;
% ax1=axes('position',[0.3,0.40,0.4,0.2])  ;
% axesm eckert4; 
% framem;
% gridm;
% axis off
geoidrefvec=[1/aj,90,180];  
B=geoshow(pm, geoidrefvec, 'DisplayType', 'texturemap'); hold on;
[c,h] = contourm(Lat,Lon,mola128_dx,-40000:20000:40000,'k'); hold on;
text(2.5,0.9,'45^o')
text(2.75,0,'0^o')
text(2.4,-0.9,'-45^o')

text(-3.1,0.9,'45^o')
text(-3.1,0,'0^o')
text(-3.2,-0.9,'-45^o')

text(-1.7,-1.5,'0^o')
text(-0.3,-1.5,'180^o')
text(1.5,-1.5,'360^o')
b = colorbar('location','southoutside','position',[0.63 0.635 0.2 0.02]);
b.FontSize=10;
h1=hotres;
% h1=hot(100); 
% h1(99,:)=[0.9 0.9 0.9];
% colormap(gca,flipud(h1(1:2:100,:) ))
colormap(gca,h1(1:2:254,:))%24/127*2=0.189*2=0.38
b.Label.String = 'M(A/m)';
b.Label.Position = [28, 1.1, 0]
b.Ticks=[0 8 16 24];
b.TickLabels=[{'0' ;'8' ;'16';'24' }];
clim([0 24]);