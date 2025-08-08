% 示例数据
clear;close all;clc
%path='D:\datap\';
load D:\datapro\1997AB.txt
load D:\datapro\2006MGS.txt
load D:\datapro\2014MAVEN.txt
load D:\datapro\2015MAVEN.txt
load D:\datapro\2016MAVEN.txt
load D:\datapro\2017MAVEN.txt
load D:\datapro\2018MAVEN.txt
load D:\datapro\2019MAVEN.txt
load D:\datapro\2020MAVEN.txt
load D:\datapro\2021MAVEN.txt
load D:\datapro\2022MAVEN.txt
MAVEN=[X2014MAVEN;X2015MAVEN;X2016MAVEN;X2017MAVEN;X2018MAVEN;X2019MAVEN;X2020MAVEN;X2021MAVEN;X2022MAVEN];
MGS=[X1997AB;X2006MGS];
clear X2014MAVEN X2015MAVEN X2016MAVEN X2017MAVEN X2018MAVEN X2019MAVEN X2020MAVEN X2021MAVEN X2022MAVEN
load D:\datapro\tianwen.txt

j1=zeros(60,181);
j2=zeros(60,181);
j3=zeros(60,181);
r1=MGS(:,1)-3393.5;
r2=MAVEN(:,1)-3393.5;
r3=tianwen(:,1)-3393.5;
MGS_r=zeros(61,1);
MAVEN_r=zeros(61,1);
tianwen_r=zeros(61,1);
for g=1:length(r1)
    rj=round(r1(g)/10);
    ej=180-round(MGS(g,2)/pi*180);
    MGS_r(rj)=MGS_r(rj)+1;
    j1(rj,ej)=j1(rj,ej)+1;
end
for g=1:length(r2)
    rj=round(r2(g)/10);
    ej=180-round(MAVEN(g,2)/pi*180);
    MAVEN_r(rj)=MAVEN_r(rj)+1;
    j2(rj,ej)=j2(rj,ej)+1;
end
for g=1:length(r3)
    rj=round(r3(g)/10);
    ej=180-round(tianwen(g,2)/pi*180);
    tianwen_r(rj)=tianwen_r(rj)+1;
    j3(rj,ej)=j3(rj,ej)+1;
end
jj1=j1+j2+j3;
j1(j1==0)=NaN;
j2(j2==0)=NaN;
j3(j3==0)=NaN;
jj1(jj1==0)=NaN;
%pos=[1 1 25 20];
%set(gcf,'unit','centimeters','position',pos)
%    axes('position',[0.08,0.7,0.32,0.20])
%%
a=figure('Color','White','position',[424,100,700,500]);
subplot(1,2,2)
x = 0:1:60;


% 将数据合并到一个矩阵中
counts = [MAVEN_r'; MGS_r'; tianwen_r'];

% 创建堆叠条形图
b=bar(x, counts, 'stacked', 'EdgeColor','None');
% for k = 1:length(b) 
%     b(k).BarWidth = 1; % 设置条形宽度为1 
%     b(k).XOffset = -0.15 + (k-1) * 0.15; % 调整每个条形的水平位置 
% end
xlim([0,60])
% 添加标签和标题
xlabel('Altitude(km)','fontsize',11);
ylabel('Counts','fontsize',11);
ax = gca; ax.YAxisLocation = 'right';
 set(gca,'Xtick',[0 10 20 30 40 50 60],'Xticklabel', {'0' '100' '200' '300' '400' '500' '600'} )
% 设置图例

hLegend=legend([b(2) b(1) b(3)], {'MGS', 'MAVEN', 'Tianwen-1'}, 'Location', 'northeast');
 set(hLegend, 'FontSize', 9);
 set(hLegend, 'Position', [0.605, 0.847, 0.1, 0.05]);
 set(hLegend, 'EdgeColor', 'none');
 text(-10,85000,'(b)','fontsize',13,'FontWeight','bold')
% 显示图形
grid on;


xj=91;
yj=181;
j4=zeros(xj,yj);
j5=zeros(xj,yj);
j6=zeros(xj,yj);

for g=1:length(r1) 
    wj=91-round(MGS(g,2)/pi*180 /(180/xj) );
    jj=round(MGS(g,3)/pi*180 /(360/yj) )+round(yj/2);
    j4(wj,jj)=j4(wj,jj)+1;
end
for g=1:length(r2) 
    wj=91-round(MAVEN(g,2)/pi*180 /(180/xj) );
    jj=round(MAVEN(g,3)/pi*180 /(360/yj) )+round(yj/2);
    j5(wj,jj)=j5(wj,jj)+1;
end
for g=1:length(r3) 
    wj=91-round(tianwen(g,2)/pi*180 /(180/xj) );
    jj=round(tianwen(g,3)/pi*180 /(360/yj) )+round(yj/2);
    j6(wj,jj)=j6(wj,jj)+1;
end
jj2=[j6(2:xj,2:yj);j5(2:xj,2:yj);j4(2:xj,2:yj)];
j4(j4==0)=NaN;
j5(j5==0)=NaN;
j6(j6==0)=NaN;
jj2(jj2==0)=NaN; 
subplot(1,2,1)
   pcolor(jj2)
   shading flat;
yline([90,180]);
   
    set(gca,'Ytick',[1 45 90 135 180 225 270],'Yticklabel', {'-90' '0' '∓90' '0' '∓90' '0' '90'} )
   set(gca,'Xtick',[1 45 90 135 180],'Xticklabel', {'0' '90' '180' '270' '360'} )
   
%    h_2=  colorbar;
%      h_2.Label.String={'counts'};
%      h_2.Ticks=[0 20 40 60];
       clim([0 250]);
%      h_2.TickLabels=[{'0' ;'20';'40';'60' }];
 ylabel('Latitude (°)','fontsize',11);
 xlabel('Longitude (\circ)','fontsize',11); 
 text(-36,270/16+270,'(a)','fontsize',13,'FontWeight','bold')
 text(150,277,'MGS')
 text(142,177,'MAVEN')
 text(132,87,'Tianwen-1')
 set(gca,'tickdir','out')
set(gca,'xminortick','on');
set(gca,'yminortick','on');

% set(gca,'ticklength',[0.05 0.025]);
  % text(280,15,'MGS','fontsize',10)
    grid on;


%  cmap=hot(256);
%  cmap=flip(cmap(64:192,:));
%  colormap(cmap);

  cmap=load('MPL_Oranges.rgb');
%  cmap=flip(cmap(64:192,:));
 colormap(cmap(33:128,:));

    h_2 = colorbar('Position', [0.495, 0.25, 0.02, 0.5]); 
    h_2.Label.String = 'Counts';     
    h_2.FontSize = 10;
    h_2.Label.Rotation = 0;
    h_2.Label.Units = 'normalized'; 
    h_2.Label.Position = [0.5, 1.1, 0];
     h_2.Ticks=[0 50 100 150 200 250];
    h_2.TickLabels={'0' ;'50';'100';'150';'200';'250'};

   % h_2.TickLabelLocation='left';
