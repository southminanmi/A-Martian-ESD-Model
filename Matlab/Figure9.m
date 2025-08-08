clear;close all;clc
misfit=load('C:\Users\41047\weight_IMF_dawn.txt');
misfit1=load('C:\Users\41047\weight_IMF_dusk.txt');
load mola128_dx;
br=flipud(reshape(misfit(:,1),[400,200])')*21;
bt=flipud(reshape(misfit(:,2),[400,200])')*21;
bp=flipud(reshape(misfit(:,3),[400,200])')*21;
br1=flipud(reshape(misfit1(:,1),[400,200])')*21;clim([7 21]);
bt1=flipud(reshape(misfit1(:,2),[400,200])')*21;
bp1=flipud(reshape(misfit1(:,3),[400,200])')*21;
aj=0.9;
rm=3393.5;
ej=0.9;
Az=aj:aj:360;
El=ej:ej:180;
az=Az/180*pi-pi;
el=El/180*pi;
[xq,yq] = meshgrid(el,az);
      %%
 a=figure('Color','White');%,'position',[424,100,560,520]);  
   pos=[1 1 20 20];
set(gcf,'unit','centimeters','position',pos)

ax1=axes('position',[0.13,0.67,0.4,0.2])  ; 
axesm eckert4; 
framem;
gridm;
axis off
title('Residual in Dawn','FontSize',12,'FontWeight','bold')


% text(2.5,0.9,'45^o','FontSize',9)
% text(2.75,0,'0^o','FontSize',9)
% text(2.4,-0.9,'-45^o','FontSize',9)

text(-2.85,0.9,'45^o','FontSize',9)
text(-3.0,0,'0^o','FontSize',9)
text(-2.9,-0.9,'-45^o','FontSize',9)

text(-1.7,-1.5,'0^o','FontSize',9)
text(-0.3,-1.5,'180^o','FontSize',9)
text(1.5,-1.5,'360^o','FontSize',9)
geoidrefvec=[1/aj,90,180];  
geoshow(br, geoidrefvec, 'DisplayType', 'texturemap');hold on;
geoidrefvec=[2,90,180]; 
contourm(mola128_dx,geoidrefvec,-40000:20000:40000,'k'); 
%colormap(mycolormap4)
text(-0.15,-2,'B_r','FontSize',11);
clim([6 14]);

ax1=axes('position',[0.13,0.42,0.4,0.2])  ; 
axesm eckert4; 
framem;
gridm;
axis off
geoidrefvec=[1/aj,90,180];  
B=geoshow(bt, geoidrefvec, 'DisplayType', 'texturemap'); hold on;
geoidrefvec=[2,90,180]; 
contourm(mola128_dx,geoidrefvec,-40000:20000:40000,'k'); 
text(-0.15,-2,'B_\theta','FontSize',11);
% text(2.5,0.9,'45^o','FontSize',9)
% text(2.75,0,'0^o','FontSize',9)
% text(2.4,-0.9,'-45^o','FontSize',9)

text(-2.85,0.9,'45^o','FontSize',9)
text(-3.0,0,'0^o','FontSize',9)
text(-2.9,-0.9,'-45^o','FontSize',9)

text(-1.7,-1.5,'0^o','FontSize',9)
text(-0.3,-1.5,'180^o','FontSize',9)
text(1.5,-1.5,'360^o','FontSize',9)
clim([6 14]);

ax1=axes('position',[0.13,0.17,0.4,0.2])  ; 
axesm eckert4; 
framem;
gridm;
axis off
geoidrefvec=[1/aj,90,180];  
B=geoshow(bp, geoidrefvec, 'DisplayType', 'texturemap'); hold on;
geoidrefvec=[2,90,180]; 
contourm(mola128_dx,geoidrefvec,-40000:20000:40000,'k'); 
%colormap(mycolormap4)
text(-0.15,-2,'B_\phi','FontSize',11);
% text(2.5,0.9,'45^o','FontSize',9)
% text(2.75,0,'0^o','FontSize',9)
% text(2.4,-0.9,'-45^o','FontSize',9)

text(-2.85,0.9,'45^o','FontSize',9)
text(-3.0,0,'0^o','FontSize',9)
text(-2.9,-0.9,'-45^o','FontSize',9)


text(-1.7,-1.5,'0^o','FontSize',9)
text(-0.3,-1.5,'180^o','FontSize',9)
text(1.5,-1.5,'360^o','FontSize',9)
clim([6 14]);

ax1=axes('position',[0.53,0.67,0.4,0.2])  ;

axesm eckert4; 
framem;
gridm;
axis off
title('Residual in Dusk','FontSize',12,'FontWeight','bold')
geoidrefvec=[1/aj,90,180];  
geoshow(br1, geoidrefvec, 'DisplayType', 'texturemap');hold on;
geoidrefvec=[2,90,180]; 
contourm(mola128_dx,geoidrefvec,-40000:20000:40000,'k'); 
%colormap(mycolormap4)
text(-0.15,-2,'B_r','FontSize',11);
text(2.5,0.9,'45^o','FontSize',9)
text(2.75,0,'0^o','FontSize',9)
text(2.4,-0.9,'-45^o','FontSize',9)

text(-3.1,0.9,'45^o','FontSize',9)
text(-3.1,0,'0^o','FontSize',9)
text(-3.2,-0.9,'-45^o','FontSize',9)

text(-1.7,-1.5,'0^o','FontSize',9)
text(-0.3,-1.5,'180^o','FontSize',9)
text(1.5,-1.5,'360^o','FontSize',9)
clim([6 14]);
ax1=axes('position',[0.53,0.42,0.4,0.2])  ;
axesm eckert4; 
framem;
gridm;
axis off
geoidrefvec=[1/aj,90,180];  
B=geoshow(bt1, geoidrefvec, 'DisplayType', 'texturemap'); hold on;
geoidrefvec=[2,90,180]; 
contourm(mola128_dx,geoidrefvec,-40000:20000:40000,'k'); 
text(-0.15,-2,'B_\theta','FontSize',11);
text(2.5,0.9,'45^o','FontSize',9)
text(2.75,0,'0^o','FontSize',9)
text(2.4,-0.9,'-45^o','FontSize',9)

text(-3.1,0.9,'45^o','FontSize',9)
text(-3.1,0,'0^o','FontSize',9)
text(-3.2,-0.9,'-45^o','FontSize',9)

text(-1.7,-1.5,'0^o','FontSize',9)
text(-0.3,-1.5,'180^o','FontSize',9)
text(1.5,-1.5,'360^o','FontSize',9)
clim([6 14]);
ax1=axes('position',[0.53,0.17,0.4,0.2])  ;
axesm eckert4; 
framem;
gridm;
axis off
geoidrefvec=[1/aj,90,180];  
B=geoshow(bp1, geoidrefvec, 'DisplayType', 'texturemap'); hold on;
geoidrefvec=[2,90,180]; 
contourm(mola128_dx,geoidrefvec,-40000:20000:40000,'k'); 
%colormap(mycolormap4)
text(-0.15,-2,'B_\phi','FontSize',11);
text(2.5,0.9,'45^o','FontSize',9)
text(2.75,0,'0^o','FontSize',9)
text(2.4,-0.9,'-45^o','FontSize',9)

text(-3.1,0.9,'45^o','FontSize',9)
text(-3.1,0,'0^o','FontSize',9)
text(-3.2,-0.9,'-45^o','FontSize',9)

text(-1.7,-1.5,'0^o','FontSize',9)
text(-0.3,-1.5,'180^o','FontSize',9)
text(1.5,-1.5,'360^o','FontSize',9)
clim([6 14]);
hcb = colorbar('location','southoutside','position',[0.27 0.1 0.5 0.02]);
cm=hot;
cmap = crameri('lajolla');
colormap(cmap);
%colormap(flip(cm(49:240,:)))
hcb.Label.String = 'nT';
hcb.FontSize=10;
hcb.Label.Position = [14.2, 1.1, 0];
%export_fig misfit_MSO_IMF.png -r600