clear all;close all;clc

rr0 = load('C:\Users\41047\res.txt');
load C:\Users\41047\r_l2.txt
rm=3393.5;




  % wsph= qd1_res.Bsph_res_m14';
  % wpc= qd1_res.Bpc_res_m14';
   
%   wsph= qd1_res.Bsph_res_m14';
%   wpc= qd1_res.Bpc_res_m14';
   

   jie=229038;
   wsphmav=rr0(:,jie:end-3060);
 
   wsphmgs=rr0(:,1:jie);



   wsphmav_l2=r_l2(:,jie:end-3060);
   
   wsphmgs_l2=r_l2(:,1:jie);
   
    
   %%
  %  Bt=(wpc(1,:).^2+wpc(2,:).^2+wpc(3,:).^2).^0.5;
  %  stot=std(Bt)
   %%
   
   figure('Color','White');
   subplot(2,3,1)
   s1=std(wsphmav(1,:));
   edg=linspace(-3*s1,3*s1,50);
  y = normpdf(edg, 0 ,s1); 
  h= histogram(wsphmav(1,:),edg,'Normalization','pdf'); hold on;
  plot(edg, y, 'LineWidth', 1.5, 'Color', [0.5, 0, 0.5]);
 binCenters = (edg(1:end-1) + edg(2:end))/2;
   h.EdgeColor=[0 0 1] 
  h.FaceColor=[1 1 1]
 % 绘制折线图
 s11=std(wsphmav_l2(1,:));
   edg1=linspace(-3*s11,3*s11,50);
 counts= histcounts(wsphmav_l2(1,:),edg1,'Normalization','pdf');
 plot(binCenters, counts,'LineWidth', 1.5, 'Color', [1 0 0]);
  
   xlim([-3*s1,3*s1]);
   ylim([0,0.1]);
   
   
   xticks([-2*s1:s1:2*s1]);
   yticks([0:0.01:0.1]);
   
   xticklabels({'-2\sigma','-\sigma','0','\sigma','2\sigma'});
   yticklabels({'0','1','2','3','4','5','6','7','8','9','10'});
   
   title('MAVEN B_r misfit','FontWeight','bold')
   ylabel('Relative counts (%)')
%[mu,sigma]=normfit(wsphmav(1,:));
set(gca,'xminortick','on')
set(gca,'yminortick','on')
set(gca,'tickdir','out')
grid on;
   
   subplot(2,3,2)
   s2=std(wsphmav(2,:));
   edg=linspace(-3*s2,3*s2,50);
  y = normpdf(edg, 0 ,s2); 
  h= histogram(wsphmav(2,:),edg,'Normalization','pdf'); hold on;
  plot(edg, y, 'LineWidth', 1.5, 'Color', [0.5, 0, 0.5]);
 binCenters = (edg(1:end-1) + edg(2:end))/2;
   h.EdgeColor=[0 0 1] 
  h.FaceColor=[1 1 1]
 % 绘制折线图
 s22=std(wsphmav_l2(2,:));
   edg1=linspace(-3*s22,3*s22,50);
 counts= histcounts(wsphmav_l2(2,:),edg1,'Normalization','pdf');
 plot(binCenters, counts,'LineWidth', 1.5, 'Color', [1 0 0]);
  
   xlim([-3*s2,3*s2]);
   ylim([0,0.1]);
   
   
   xticks([-2*s2:s2:2*s2]);
    yticks([0:0.01:0.1]);
   
   xticklabels({'-2\sigma','-\sigma','0','\sigma','2\sigma'});
   yticklabels({'0','1','2','3','4','5','6','7','8','9','10'});

   title('MAVEN B_\theta misfit','FontWeight','bold')
   ylabel('Relative counts (%)')
   set(gca,'xminortick','on')
set(gca,'yminortick','on')
set(gca,'tickdir','out')
grid on;
   
   subplot(2,3,3)
       s3=std(wsphmav(3,:));
   edg=linspace(-3*s3,3*s3,50);
  y = normpdf(edg, 0 ,s3); 
  h= histogram(wsphmav(3,:),edg,'Normalization','pdf'); hold on;
  plot(edg, y, 'LineWidth', 1.5, 'Color', [0.5, 0, 0.5]);
 binCenters = (edg(1:end-1) + edg(2:end))/2;
   h.EdgeColor=[0 0 1] 
  h.FaceColor=[1 1 1]
 % 绘制折线图
 s33=std(wsphmav_l2(3,:));
   edg1=linspace(-3*s33,3*s33,50);
 counts= histcounts(wsphmav_l2(3,:),edg1,'Normalization','pdf');
 plot(binCenters, counts,'LineWidth', 1.5, 'Color', [1 0 0]);
 hLegend=legend('Huber weight','Norm','Least Squares');
 set(hLegend, 'FontSize', 7);
set(hLegend, 'Location', 'northeast');
set(hLegend, 'Position', [0.79, 0.85, 0.1, 0.05]);

   xlim([-3*s3,3*s3]);
   ylim([0,0.1]);
  
   xticks([-2*s3:s3:2*s3]);
    yticks([0:0.01:0.1]);
   
   xticklabels({'-2\sigma','-\sigma','0','\sigma','2\sigma'});
   yticklabels({'0','1','2','3','4','5','6','7','8','9','10'});
   
   title('MAVEN B_\phi misfit','FontWeight','bold')
   ylabel('Relative counts (%)')
set(gca,'xminortick','on')
set(gca,'yminortick','on')
set(gca,'tickdir','out')
grid on;
   
   
   %%
      
   subplot(2,3,4)
   s1=std(wsphmgs(1,:));
   edg=linspace(-3*s1,3*s1,50);
  y = normpdf(edg, 0 ,s1); 
  h= histogram(wsphmgs(1,:),edg,'Normalization','pdf'); hold on;
  plot(edg, y, 'LineWidth', 1.5, 'Color', [0.5, 0, 0.5]);
 binCenters = (edg(1:end-1) + edg(2:end))/2;
   h.EdgeColor=[0 0 1] 
  h.FaceColor=[1 1 1]
 % 绘制折线图
 s11=std(wsphmgs_l2(1,:));
   edg1=linspace(-3*s11,3*s11,50);
 counts= histcounts(wsphmgs_l2(1,:),edg1,'Normalization','pdf');
 plot(binCenters, counts,'LineWidth', 1.5, 'Color', [1 0 0]);

  
   xlim([-3*s1,3*s1]);
   ylim([0,0.1]);
   
   
   xticks([-2*s1:s1:2*s1]);
   yticks([0:0.01:0.1]);
   
   xticklabels({'-2\sigma','-\sigma','0','\sigma','2\sigma'});
   yticklabels({'0','1','2','3','4','5','6','7','8','9','10'});
   
   title('MGS B_r misfit','FontWeight','bold')
   ylabel('Relative counts (%)')
   set(gca,'xminortick','on')
set(gca,'yminortick','on')
set(gca,'tickdir','out')
grid on;
%[mu,sigma]=normfit(wsphmgs(1,:));

   
   subplot(2,3,5)
  s2=std(wsphmgs(2,:));
   edg=linspace(-3*s2,3*s2,50);
  y = normpdf(edg, 0 ,s2); 
  h= histogram(wsphmgs(2,:),edg,'Normalization','pdf'); hold on;
  plot(edg, y, 'LineWidth', 1.5, 'Color', [0.5, 0, 0.5]);
 binCenters = (edg(1:end-1) + edg(2:end))/2;
   h.EdgeColor=[0 0 1] 
  h.FaceColor=[1 1 1]
 % 绘制折线图
 s22=std(wsphmgs_l2(2,:));
   edg1=linspace(-3*s22,3*s22,50);
 counts= histcounts(wsphmgs_l2(2,:),edg1,'Normalization','pdf');
 plot(binCenters, counts,'LineWidth', 1.5, 'Color', [1 0 0]);

   xlim([-3*s2,3*s2]);
   ylim([0,0.1]);
   
   xticks([-2*s2:s2:2*s2]);
   yticks([0:0.01:0.1]);
   
   xticklabels({'-2\sigma','-\sigma','0','\sigma','2\sigma'});
   yticklabels({'0','1','2','3','4','5','6','7','8','9','10'});
   
   title('MGS B_\theta misfit','FontWeight','bold')
   ylabel('Relative counts (%)')
   set(gca,'xminortick','on')
set(gca,'yminortick','on')
set(gca,'tickdir','out')
grid on;
   
   subplot(2,3,6)
     s3=std(wsphmgs(3,:));
   edg=linspace(-3*s3,3*s3,50);
  y = normpdf(edg, 0 ,s3); 
  h= histogram(wsphmgs(3,:),edg,'Normalization','pdf'); hold on;
  plot(edg, y, 'LineWidth', 1.5, 'Color', [0.5, 0, 0.5]);
 binCenters = (edg(1:end-1) + edg(2:end))/2;
   h.EdgeColor=[0 0 1] 
  h.FaceColor=[1 1 1]
 % 绘制折线图
 s33=std(wsphmgs_l2(3,:));
   edg1=linspace(-3*s33,3*s33,50);
 counts= histcounts(wsphmgs_l2(3,:),edg1,'Normalization','pdf');
 plot(binCenters, counts,'LineWidth', 1.5, 'Color', [1 0 0]);

   xlim([-3*s3,3*s3]);
   ylim([0,0.1]);
   
   
   xticks([-2*s3:s3:2*s3]);
   yticks([0:0.01:0.1]);
   
   xticklabels({'-2\sigma','-\sigma','0','\sigma','2\sigma'});
   yticklabels({'0','1','2','3','4','5','6','7','8','9','10'});
   
   title('MGS B_\phi misfit','FontWeight','bold')
   ylabel('Relative counts (%)')
   
set(gca,'xminortick','on')
set(gca,'yminortick','on')
set(gca,'tickdir','out')

  grid on; 
   
   %%
  pos=[0 0 25 15];
 set(gcf,'PaperPositionMode','manual','PaperUnits','centimeters',...
         'PaperPosition',pos,'units','centimeters','position',pos);
    









