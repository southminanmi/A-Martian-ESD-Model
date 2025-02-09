clear;close all;clc
maxD = 160;
 n=maxD^2+2*maxD;
 i=0;
  load('C:\Users\41047\r.txt');
load('C:\Users\41047\theta.txt');
load('C:\Users\41047\phi.txt');
step=10000;
count=1:step:length(r);
G=zeros(n,n);
d=zeros(n,1);
load('C:\Users\41047\br.txt');
load('C:\Users\41047\bt.txt');
load('C:\Users\41047\bp.txt');
 b=[br,bt,bp];
%[A_r_gsm, A_theta_gsm, A_phi_gsm]=synth_values(r, theta, phi, cof)
for i=1:length(count)
    if i~=length(count)
    [A_r_gsm, A_theta_gsm, A_phi_gsm] = design_SHA(r((i-1)*step+1:i*step), theta((i-1)*step+1:i*step),phi((i-1)*step+1:i*step),maxD,'int');
    X=[A_r_gsm; A_theta_gsm; A_phi_gsm];
    d=d+X'*reshape(b((i-1)*step+1:i*step,:),[],1);
    G=G+X'*X;
    else
    [A_r_gsm, A_theta_gsm, A_phi_gsm] = design_SHA(r((i-1)*step+1:length(r)), theta((i-1)*step+1:length(r)),phi((i-1)*step+1:length(r)),maxD,'int');
    X=[A_r_gsm; A_theta_gsm; A_phi_gsm];
    d=d+X'*reshape(b((i-1)*step+1:length(r),:),[],1);
    G=G+X'*X;  
    end
end
x0 = zeros(n,1);
% x0=G\d;
r0 = G*x0 - d;
p0 = -r0;
iter_max = 200;
n=0;
for i =1:iter_max
    alpha = (r0'*r0)/((p0'*G)*p0);
    x = x0 + alpha*p0;
    r = r0 + alpha*(G*p0);
    beta = (r'*r)/(r0'*r0);
    p = -r + beta*p0;
    if norm(x-x0,2)/norm(x0,2)<=1e-3
       break
    end
    n=n+1
    x0 = x;
    r0 = r;
    p0 = p;
    m0=x0;
end
%export_fig spectrum.png -r600