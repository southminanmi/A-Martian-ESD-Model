# -*- coding: utf-8 -*-
"""
Created on Tue Sep  5 15:21:51 2023

@author: 41047
"""

import pandas as pd
import numpy as np
import os
import cupy as cp
os.environ['KMP_DUPLICATE_LIB_OK']='TRUE'
path = "D:\\datapro"
files= os.listdir(path) 
df = pd.DataFrame(columns=['r','theta','phi','br','btheta','bphi','s_theta','s_phi'])
for file in files:
    file = os.path.join(path,file)
    data = pd.read_csv(file, sep='\s+',header=None,names=['r','theta','phi','br','btheta','bphi','s_theta','s_phi'])
    df=pd.concat([df,data],axis=0)
b=df[['br','btheta','bphi']]
phi=cp.array(df['phi'],dtype='float32')
theta=cp.array(df['theta'],dtype='float32')
r=cp.array(df['r'],dtype='float32')/3393.5#Be careful the r
del df
r=r.reshape(-1,1)
phi=phi.reshape(-1,1)
theta=theta.reshape(-1,1)
# In[mesh1]
Ntheta=114#In fact,Ntheta=Ntheta+1
def generatemesh(r,Ntheta):
    points=cp.array([[0,0,r],[np.pi,0,r]])
    theta=cp.linspace(0,np.pi,Ntheta,endpoint=False,dtype='float32')
    n=0
    for lat in theta[1:]:
        Nphi=round(1/2+np.sin(cp.asnumpy(lat))*3**0.5*Ntheta)
        phi=cp.linspace(0+2*np.pi*n/Ntheta,2*np.pi+2*np.pi*n/Ntheta,Nphi,endpoint=False,dtype='float32')
        temp=cp.concatenate((lat*cp.ones((Nphi,1)),phi.reshape((Nphi,1)),r*cp.ones((Nphi,1))),axis=1)        
        points=cp.concatenate((points,temp),axis=0)
        n+=1
    return points
mesh=generatemesh(3373.5/3393.5, Ntheta)
mesh=mesh.astype(cp.float32)
# In[mesh]
N=15000
def generatemesh(r,N):
    cofphi = (np.sqrt(5) - 1) / 2
    n = cp.arange(0, N)
    z = ((2*n + 1) / N - 1)
    x = (cp.sqrt(1 - z**2)) * cp.cos(2 * np.pi * (n + 1) *cofphi)
    y = (cp.sqrt(1 - z**2)) * cp.sin(2 * np.pi * (n + 1) * cofphi)
    r=r*cp.ones((N,1))
    theta=cp.arccos(z)
    phi=cp.arctan2(y,x)
    points=cp.concatenate((theta.reshape((N,1)),phi.reshape((N,1)),r),axis=1)
    return points
mesh=generatemesh(3373.5/3393.5, N)
mesh=mesh.astype(cp.float32)
# In[dipole_to_magneticfield]
def A1(theta,thetad,phi,phid):
    result=cp.cos(theta)*cp.cos(thetad)+cp.sin(theta)*cp.sin(thetad)*cp.cos(phi-phid)
    return result
def l(theta,thetad,phi,phid,r,rd):
    result=(rd**2+r**2-2*rd*r*(cp.cos(theta)*cp.cos(thetad)+cp.sin(theta)*cp.sin(thetad)*cp.cos(phi-phid)))**0.5
    return result
def D1(theta,thetad,phi,phid,r,rd):
    result=r-rd*A1(theta,thetad,phi,phid)
    return result
    return result
def F1(theta,thetad,phi,phid,r,rd):
    result=r*A1(theta,thetad,phi,phid)-rd
    return result 
def c11(phi,theta,r,mesh):
    c11=((3*D1(theta,mesh[:,0],phi,mesh[:,1],r,mesh[:,2])*F1(theta,mesh[:,0],phi,mesh[:,1],r,mesh[:,2]))/l(theta,mesh[:,0],phi,mesh[:,1],r,mesh[:,2])**2-A1(theta,mesh[:,0],phi,mesh[:,1]))
    length=1/l(theta,mesh[:,0],phi,mesh[:,1],r,mesh[:,2])**3
    length[length<1/0.5**3]=0
    c11=c11*length
    return c11  
def generatemagneticbr(phi,theta,r,mesh,M):
    Br=cp.dot(c11(phi,theta,r,mesh),M)
    return Br
# In[inversion]
import cupy as cp
n=mesh.shape[0]
b=cp.array(b,dtype='float32')
lamda=1e7
G=cp.zeros([n,n])#+lamda*cp.diag(mesh.shape[0]*[0]+2*mesh.shape[0]*[1])
d=cp.zeros([n,1])
step=5000
count=np.arange(0,len(r),step);
for i in count:
    if i==count[-1]:
        C11=c11(phi[i:],theta[i:],r[i:],mesh)
        d=d+cp.dot(cp.transpose(C11),b[i:,0]).reshape(-1,1)
        G=G+cp.dot(cp.transpose(C11),C11)
        del C11
    else:
        C11=c11(phi[i:i+step],theta[i:i+step],r[i:i+step],mesh)
        d=d+cp.dot(cp.transpose(C11),b[i:i+step,0]).reshape(-1,1)
        G=G+cp.dot(cp.transpose(C11),C11)
        del C11
G=G
d=d
x0=cp.zeros([n,1])
r0 = cp.dot(G,x0) - d;
p0 = -r0;
iter_max = 20;
misfit=[];
n=0
for i in range(iter_max):
    alpha = cp.dot(cp.transpose(r0),r0)/cp.dot(cp.dot(cp.transpose(p0),G),p0);
    x = x0 + alpha*p0;
    r = r0 + alpha*cp.dot(G,p0);
    beta = cp.dot(cp.transpose(r),r)/cp.dot(cp.transpose(r0),r0);
    p = -r + beta*p0;
    misfit.append(cp.linalg.norm(x=r0,ord=2))
    if cp.asnumpy(cp.linalg.norm(x=r-r0,ord=2)/(1+cp.linalg.norm(x=r0,ord=2)))<=0.01:
        break
    n=n+1
    print(n)
    x0 = x;
    r0 = r;
    p0 = p;
m0=x0
# In[b]
theta=cp.linspace(0.01,np.pi-0.01,200)
phi=cp.linspace(0,2*np.pi,400,dtype='float32')
[phi,theta]=np.meshgrid(phi,theta)
phi=phi.reshape(-1,1)
theta=theta.reshape(-1,1)
r=(cp.ones(phi.shape[0]*phi.shape[1])*3493.5)/3393.5
r=r.astype(cp.float16)
r=r.reshape(-1,1)
step=10000
count=np.arange(0,len(r),step);
for i in count:
    if i==count[0]:
        br=generatemagneticbr(phi[i:i+step].reshape(step,1),theta[i:i+step].reshape(step,1),r[i:i+step].reshape(step,1),mesh,m0)
    elif i==count[-1]:
        btemp=generatemagneticbr(phi[i:i+step].reshape(step,1),theta[i:i+step].reshape(step,1),r[i:i+step].reshape(step,1),mesh,m0)
        br=cp.concatenate((br,btemp),axis=0)
    else:
        btemp=generatemagneticbr(phi[i:i+step].reshape(step,1),theta[i:i+step].reshape(step,1),r[i:i+step].reshape(step,1),mesh,m0)
        br=cp.concatenate((br,btemp),axis=0)
# In[b]
theta=(cp.random.rand(100000,1))*np.pi
phi=(cp.random.rand(100000,1))*np.pi*2
phi=phi.reshape(-1,1)
theta=theta.reshape(-1,1)
r=(cp.random.rand(100000,1)*310+3483.5)/3393.5
r=r.astype(cp.float16)
r=r.reshape(-1,1)
step=10000
count=np.arange(0,len(r),step);
for i in count:
    if i==count[0]:
        br=generatemagneticbr(phi[i:i+step].reshape(step,1),theta[i:i+step].reshape(step,1),r[i:i+step].reshape(step,1),mesh,m0)
    elif i==count[-1]:
        btemp=generatemagneticbr(phi[i:i+step].reshape(step,1),theta[i:i+step].reshape(step,1),r[i:i+step].reshape(step,1),mesh,m0)
        br=cp.concatenate((br,btemp),axis=0)
    else:
        btemp=generatemagneticbr(phi[i:i+step].reshape(step,1),theta[i:i+step].reshape(step,1),r[i:i+step].reshape(step,1),mesh,m0)
        br=cp.concatenate((br,btemp),axis=0)
'''
np.savetxt('r.txt',cp.asnumpy(r))
np.savetxt('theta.txt',cp.asnumpy(theta))
np.savetxt('phi.txt',cp.asnumpy(phi))
'''
# In[b1]
N=100000
def generatemesh(r,N):
    cofphi = (np.sqrt(5) - 1) / 2
    n = cp.arange(0, N)
    z = ((2*n + 1) / N - 1)
    x = (cp.sqrt(1 - z**2)) * cp.cos(2 * np.pi * (n + 1) *cofphi)
    y = (cp.sqrt(1 - z**2)) * cp.sin(2 * np.pi * (n + 1) * cofphi)
    r=r*cp.ones((N,1))
    theta=cp.arccos(z)
    phi=cp.arctan2(y,x)
    points=cp.concatenate((theta.reshape((N,1)),phi.reshape((N,1)),r),axis=1)
    return points
mesh1=generatemesh(3543.5/3393.5, N)
mesh1=mesh1.astype(cp.float32)
theta=mesh1[:,0]
phi=mesh1[:,1]
r=mesh1[:,2]
phi=phi.reshape(-1,1)
theta=theta.reshape(-1,1)
r=r.astype(cp.float16)
r=r.reshape(-1,1)
step=20000
count=np.arange(0,len(r),step);
for i in count:
    if i==count[0]:
        br=generatemagneticbr(phi[i:i+step].reshape(step,1),theta[i:i+step].reshape(step,1),r[i:i+step].reshape(step,1),mesh,m0)
    elif i==count[-1]:
        btemp=generatemagneticbr(phi[i:i+step].reshape(step,1),theta[i:i+step].reshape(step,1),r[i:i+step].reshape(step,1),mesh,m0)
        br=cp.concatenate((br,btemp),axis=0)
    else:
        btemp=generatemagneticbr(phi[i:i+step].reshape(step,1),theta[i:i+step].reshape(step,1),r[i:i+step].reshape(step,1),mesh,m0)
        br=cp.concatenate((br,btemp),axis=0)
np.savetxt('br.txt',cp.asnumpy(br))
# In[misfit]
c=cp.random.choice(900000,10000)+203831#+208866
#c=cp.arange(210000,220000)
phit=phi[c]
thetat=theta[c]
rt=r[c]
b=cp.array(b,dtype='float32')
brt=b[c,0]
step=5000
count=np.arange(0,len(rt),step);
for i in count:
    if i==count[0]:
        br=generatemagneticbr(phit[i:i+step].reshape(step,1),thetat[i:i+step].reshape(step,1),rt[i:i+step].reshape(step,1),mesh,m0)
    elif i==count[-1]:
        btemp=generatemagneticbr(phit[i:i+step].reshape(step,1),thetat[i:i+step].reshape(step,1),rt[i:i+step].reshape(step,1),mesh,m0)
        br=cp.concatenate((br,btemp),axis=0)
    else:
        btemp=generatemagneticbr(phit[i:i+step].reshape(step,1),thetat[i:i+step].reshape(step,1),rt[i:i+step].reshape(step,1),mesh,m0)
        br=cp.concatenate((br,btemp),axis=0)
std_br=cp.std(brt.reshape(-1,1)-br)
print("std_br=",std_br)
