clear;close all;clc
load cof
load g_110_mm_q
gao_g=g;
load h_110_mm_q
gao_h=h;
load cof_final.mat
load cof_m110.txt
maxD=160;
Rn= zeros(maxD,1);
for i =1:maxD
    Rn(i)=sum(x0(i^2:i^2+2*i).^2)*(i+1)*(3395.5/(3395.5+0))^(2*i+4);
end
maxD=134;
Rn1= zeros(maxD,1);
for i =1:maxD
        Rn1(i)=sum(cof(i^2:i^2+2*i).^2)*(i+1)*(3395.5/(3395.5+0))^(2*i+4);
end
maxD=110;
Rn2= zeros(maxD,1);
for i=1:maxD
    g2=gao_g.^2;
    h2=gao_h.^2;
    Rn2(i)=(i+1)*(3395.5/(3395.5+0))^(2*i+4)*(sum(g2(:,i))+sum(h2(:,i)));
end
maxD=110;
Rn3= zeros(maxD,1);
for i=1:maxD
    Rn3(i)=sum(cof_m110(i^2:i^2+2*i).^2)*(i+1)*(3395.5/(3395.5+0))^(2*i+4);
end
figure('Color','White','position',[424,100,800,350])
subplot(1,2,1)
plot(Rn3,'g','Linewidth',1.5);hold on;
plot(Rn1,'b','Linewidth',1.5);hold on;
plot(Rn2,'black','Linewidth',1.5);hold on;
%plot(Rn(1:110),'c','Linewidth',1.5);hold on;
plot(Rn(1:110),'r','Linewidth',1.5);hold on;
plot(111:3:160,Rn(111:3:160)','Color','r','LineStyle','none','Linewidth',0.3,'Marker','.')
%plot(111:3:160,Rn(111:3:160)','Color','c','LineStyle','none','Linewidth',0.3,'Marker','.')
ylim([1e0,1e5]);
xlabel("SH Degree")
ylabel("{R_n}(nT^2)")
title('Surface','FontSize',15,'FontWeight','bold')
set(gca,'yscale','log')
set(gca,'tickdir','out')
 % set(gca,'FontSize',10)
set(gca,'xminortick','on');
set(gca,'yminortick','on');    
grid on;
text(-35,1e5,'(a)','fontsize',15,'FontWeight','bold')
set(gca,'Xtick',[20 40 60 80 100 120 140 160] )
legend({'M110','L134','G110','This study'},'Location','south','FontSize',11)
maxD=160;
Rn= zeros(maxD,1);
for i =1:maxD
    Rn(i)=sum(x0(i^2:i^2+2*i).^2)*(i+1)*(3395.5/(3395.5+120))^(2*i+4);
end
maxD=134;
Rn1= zeros(maxD,1);
for i =1:maxD
        Rn1(i)=sum(cof(i^2:i^2+2*i).^2)*(i+1)*(3395.5/(3395.5+120))^(2*i+4);
end
maxD=110;
Rn2= zeros(maxD,1);
for i=1:maxD
    g2=gao_g.^2;
    h2=gao_h.^2;
    Rn2(i)=(i+1)*(3395.5/(3395.5+120))^(2*i+4)*(sum(g2(:,i))+sum(h2(:,i)));
end
Rn3= zeros(maxD,1);
for i=1:maxD
    Rn3(i)=sum(cof_m110(i^2:i^2+2*i).^2)*(i+1)*(3395.5/(3395.5+120))^(2*i+4);
end
subplot(1,2,2)
plot(Rn3,'g','Linewidth',1.5);hold on;
plot(Rn1,'b','Linewidth',1.5);hold on;
plot(Rn2,'black','Linewidth',1.5);hold on;
%plot(Rn(1:110),'c','Linewidth',1.5);hold on;
plot(Rn(1:110),'r','Linewidth',1.5);hold on;
plot(111:3:160,Rn(111:3:160)','Color','r','LineStyle','none','Linewidth',0.3,'Marker','.')
%plot(111:3:160,Rn(111:3:160)','Color','c','LineStyle','none','Linewidth',0.3,'Marker','.')
ylim([1e-2,1e3]);
xlabel("SH Degree")
ylabel("{R_n}(nT^2)")
title('120 km Altitude','FontSize',15,'FontWeight','bold')
set(gca,'yscale','log')
set(gca,'tickdir','out')
 % set(gca,'FontSize',10)
set(gca,'xminortick','on');
set(gca,'yminortick','on');    
grid on;
text(-35,1e3,'(b)','fontsize',15,'FontWeight','bold')
legend({'M110','L134','G110','This study'},'Location','southwest','FontSize',11)
set(gca,'Xtick',[20 40 60 80 100 120 140 160] )
%export_fig spectrum.png -r600