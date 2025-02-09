# -*- coding: utf-8 -*-
"""
Created on Wed Nov  1 13:40:50 2023

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
def B1(theta,thetad,phi,phid):
    result=cp.cos(theta)*cp.sin(thetad)-cp.sin(theta)*cp.cos(thetad)*cp.cos(phi-phid)
    return result
def C1(theta,thetad,phi,phid):
    result=cp.sin(theta)*cp.sin(phi-phid)
    return result
def A2(theta,thetad,phi,phid):
    result=-cp.sin(theta)*cp.cos(thetad)+cp.cos(theta)*cp.sin(thetad)*cp.cos(phi-phid)
    return result
def B2(theta,thetad,phi,phid):
    result=-cp.sin(theta)*cp.sin(thetad)-cp.cos(theta)*cp.cos(thetad)*cp.cos(phi-phid)
    return result
def C2(theta,thetad,phi,phid):
    result=cp.cos(theta)*cp.sin(phi-phid)
    return result
def A3(theta,thetad,phi,phid):
    result=-cp.sin(thetad)*cp.sin(phi-phid)
    return result
def B3(theta,thetad,phi,phid):
    result=cp.cos(thetad)*cp.sin(phi-phid)
    return result
def C3(theta,thetad,phi,phid):
    result=cp.cos(phi-phid)
    return result
def l(theta,thetad,phi,phid,r,rd):
    result=(rd**2+r**2-2*rd*r*(cp.cos(theta)*cp.cos(thetad)+cp.sin(theta)*cp.sin(thetad)*cp.cos(phi-phid)))**0.5
    return result
def D1(theta,thetad,phi,phid,r,rd):
    result=r-rd*A1(theta,thetad,phi,phid)
    return result
def D2(theta,thetad,phi,phid,r,rd):
    result=-rd*A2(theta,thetad,phi,phid)
    return result
def D3(theta,thetad,phi,phid,r,rd):
    result=-rd*A3(theta,thetad,phi,phid)
    return result
def F1(theta,thetad,phi,phid,r,rd):
    result=r*A1(theta,thetad,phi,phid)-rd
    return result
def F2(theta,thetad,phi,phid,r,rd):
    result=-r*B1(theta,thetad,phi,phid)
    return result 
def F3(theta,thetad,phi,phid,r,rd):
    result=r*C1(theta,thetad,phi,phid)
    return result  
def c11(phi,theta,r,mesh):
    c11=((3*D1(theta,mesh[:,0],phi,mesh[:,1],r,mesh[:,2])*F1(theta,mesh[:,0],phi,mesh[:,1],r,mesh[:,2]))/l(theta,mesh[:,0],phi,mesh[:,1],r,mesh[:,2])**2-A1(theta,mesh[:,0],phi,mesh[:,1]))
    length=1/l(theta,mesh[:,0],phi,mesh[:,1],r,mesh[:,2])**3
    length[length<1/0.5**3]=0
    c11=c11*length
    return c11
def c12(phi,theta,r,mesh):
    c12=((3*D1(theta,mesh[:,0],phi,mesh[:,1],r,mesh[:,2])*F2(theta,mesh[:,0],phi,mesh[:,1],r,mesh[:,2]))/l(theta,mesh[:,0],phi,mesh[:,1],r,mesh[:,2])**2+B1(theta,mesh[:,0],phi,mesh[:,1]))
    length=1/l(theta,mesh[:,0],phi,mesh[:,1],r,mesh[:,2])**3
    length[length<1/0.5**3]=0
    c12=c12*length
    return c12
def c13(phi,theta,r,mesh):
    c13=((3*D1(theta,mesh[:,0],phi,mesh[:,1],r,mesh[:,2])*F3(theta,mesh[:,0],phi,mesh[:,1],r,mesh[:,2]))/l(theta,mesh[:,0],phi,mesh[:,1],r,mesh[:,2])**2-C1(theta,mesh[:,0],phi,mesh[:,1]))
    length=1/l(theta,mesh[:,0],phi,mesh[:,1],r,mesh[:,2])**3
    length[length<1/0.5**3]=0
    c13=c13*length
    return c13
def c21(phi,theta,r,mesh):
    c21=((3*D2(theta,mesh[:,0],phi,mesh[:,1],r,mesh[:,2])*F1(theta,mesh[:,0],phi,mesh[:,1],r,mesh[:,2]))/l(theta,mesh[:,0],phi,mesh[:,1],r,mesh[:,2])**2-A2(theta,mesh[:,0],phi,mesh[:,1]))
    length=1/l(theta,mesh[:,0],phi,mesh[:,1],r,mesh[:,2])**3
    length[length<1/0.5**3]=0
    c21=c21*length
    return c21
def c22(phi,theta,r,mesh):
    c22=((3*D2(theta,mesh[:,0],phi,mesh[:,1],r,mesh[:,2])*F2(theta,mesh[:,0],phi,mesh[:,1],r,mesh[:,2]))/l(theta,mesh[:,0],phi,mesh[:,1],r,mesh[:,2])**2+B2(theta,mesh[:,0],phi,mesh[:,1]))
    length=1/l(theta,mesh[:,0],phi,mesh[:,1],r,mesh[:,2])**3
    length[length<1/0.5**3]=0
    c22=c22*length
    return c22
def c23(phi,theta,r,mesh):
    c23=((3*D2(theta,mesh[:,0],phi,mesh[:,1],r,mesh[:,2])*F3(theta,mesh[:,0],phi,mesh[:,1],r,mesh[:,2]))/l(theta,mesh[:,0],phi,mesh[:,1],r,mesh[:,2])**2-C2(theta,mesh[:,0],phi,mesh[:,1]))
    length=1/l(theta,mesh[:,0],phi,mesh[:,1],r,mesh[:,2])**3
    length[length<1/0.5**3]=0
    c23=c23*length
    return c23
def c31(phi,theta,r,mesh):
    c31=((3*D3(theta,mesh[:,0],phi,mesh[:,1],r,mesh[:,2])*F1(theta,mesh[:,0],phi,mesh[:,1],r,mesh[:,2]))/l(theta,mesh[:,0],phi,mesh[:,1],r,mesh[:,2])**2-A3(theta,mesh[:,0],phi,mesh[:,1]))
    length=1/l(theta,mesh[:,0],phi,mesh[:,1],r,mesh[:,2])**3
    length[length<1/0.5**3]=0
    c31=c31*length
    return c31
def c32(phi,theta,r,mesh):
    c32=((3*D3(theta,mesh[:,0],phi,mesh[:,1],r,mesh[:,2])*F2(theta,mesh[:,0],phi,mesh[:,1],r,mesh[:,2]))/l(theta,mesh[:,0],phi,mesh[:,1],r,mesh[:,2])**2+B3(theta,mesh[:,0],phi,mesh[:,1]))
    length=1/l(theta,mesh[:,0],phi,mesh[:,1],r,mesh[:,2])**3
    length[length<1/0.5**3]=0
    c32=c32*length
    return c32
def c33(phi,theta,r,mesh):
    c33=((3*D3(theta,mesh[:,0],phi,mesh[:,1],r,mesh[:,2])*F3(theta,mesh[:,0],phi,mesh[:,1],r,mesh[:,2]))/l(theta,mesh[:,0],phi,mesh[:,1],r,mesh[:,2])**2-C3(theta,mesh[:,0],phi,mesh[:,1]))
    length=1/l(theta,mesh[:,0],phi,mesh[:,1],r,mesh[:,2])**3
    length[length<1/0.5**3]=0
    c33=c33*length
    return c33         
def generatemagnetic(phi,theta,r,mesh,M):
    Br=cp.dot(c11(phi,theta,r,mesh),M[:,0])
    +cp.dot(c12(phi,theta,r,mesh),M[:,1])
    +cp.dot(c13(phi,theta,r,mesh),M[:,2])
    Btheta=cp.dot(c21(phi,theta,r,mesh),M[:,0])
    +cp.dot(c22(phi,theta,r,mesh),M[:,1])
    +cp.dot(c23(phi,theta,r,mesh),M[:,2])
    Bphi=cp.dot(c31(phi,theta,r,mesh),M[:,0])
    +cp.dot(c32(phi,theta,r,mesh),M[:,1])
    +cp.dot(c33(phi,theta,r,mesh),M[:,2])
    return [Br,Btheta,Bphi]
# In[inversion]
import cupy as cp
n=mesh.shape[0]*mesh.shape[1]
b=cp.array(b,dtype='float32')
# In[Huber]
from tqdm import tqdm 
x0=cp.zeros([mesh.shape[0]*3,1])
#x0=cp.array(np.loadtxt('C:/Users/41047/mfinal_10.txt')).reshape(-1,1)
#x0=cp.array(np.loadtxt('C:/Users/41047/m_inverse.txt')).reshape(-1,1)
x0[0:mesh.shape[0]]=cp.array(np.loadtxt('C:/Users/41047/mr_final.txt')).reshape(-1,1)
G=cp.zeros([n,n])#+lamda*cp.diag(mesh.shape[0]*[0]+2*mesh.shape[0]*[1])
d=cp.zeros([n,1])
step=5000
nr=len(r)
delta_MGS=6*1.32
delta_MAVEN=7*1.32
delta_tianwen=10*1.32
#delta=cp.array(203831*[9]+200564*[13]+30498*[9])
delta=cp.array((203831+25206)*[delta_MGS]+869868*[delta_MAVEN]+3060*[delta_tianwen])
count=np.arange(0,len(r),step);
for i in tqdm(count):
    if i==count[-1]:
        C11=c11(phi[i:],theta[i:],r[i:],mesh)
        C12=c12(phi[i:],theta[i:],r[i:],mesh)
        C13=c13(phi[i:],theta[i:],r[i:],mesh)
        X1=cp.column_stack((C11,C12,C13))
        del C11,C12,C13
        d=d+cp.dot(cp.transpose(X1),b[i:,0]).reshape(-1,1)
        G=G+cp.dot(cp.transpose(X1),X1)
        del X1
        C21=c21(phi[i:],theta[i:],r[i:],mesh)
        C22=c22(phi[i:],theta[i:],r[i:],mesh)
        C23=c23(phi[i:],theta[i:],r[i:],mesh)
        X2=cp.column_stack((C21,C22,C23))
        del C21,C22,C23
        d=d+cp.dot(cp.transpose(X2),b[i:,1]).reshape(-1,1)
        G=G+cp.dot(cp.transpose(X2),X2)
        del X2
        C31=c31(phi[i:],theta[i:],r[i:],mesh)
        C32=c32(phi[i:],theta[i:],r[i:],mesh)
        C33=c33(phi[i:],theta[i:],r[i:],mesh)
        X3=cp.column_stack((C31,C32,C33))
        del C31,C32,C33
        d=d+cp.dot(cp.transpose(X3),b[i:,2]).reshape(-1,1)
        G=G+cp.dot(cp.transpose(X3),X3)
        del X3
    else:
        C11=c11(phi[i:i+step],theta[i:i+step],r[i:i+step],mesh)
        C12=c12(phi[i:i+step],theta[i:i+step],r[i:i+step],mesh)
        C13=c13(phi[i:i+step],theta[i:i+step],r[i:i+step],mesh)
        X1=cp.column_stack((C11,C12,C13))
        del C11,C12,C13
        d=d+cp.dot(cp.transpose(X1),b[i:i+step,0]).reshape(-1,1)
        G=G+cp.dot(cp.transpose(X1),X1)
        del X1
        C21=c21(phi[i:i+step],theta[i:i+step],r[i:i+step],mesh)
        C22=c22(phi[i:i+step],theta[i:i+step],r[i:i+step],mesh)
        C23=c23(phi[i:i+step],theta[i:i+step],r[i:i+step],mesh)
        X2=cp.column_stack((C21,C22,C23))
        del C21,C22,C23
        d=d+cp.dot(cp.transpose(X2),b[i:i+step,1]).reshape(-1,1)
        G=G+cp.dot(cp.transpose(X2),X2)
        del X2
        C31=c31(phi[i:i+step],theta[i:i+step],r[i:i+step],mesh)
        C32=c32(phi[i:i+step],theta[i:i+step],r[i:i+step],mesh)
        C33=c33(phi[i:i+step],theta[i:i+step],r[i:i+step],mesh)
        X3=cp.column_stack((C31,C32,C33))
        del C31,C32,C33
        d=d+cp.dot(cp.transpose(X3),b[i:i+step,2]).reshape(-1,1)
        G=G+cp.dot(cp.transpose(X3),X3)
r0 = cp.dot(G,x0) - d;
iter_max2 = 1;
misfit=[];
number=0
p0 = -r0;
b=cp.concatenate((b[:,0].reshape(-1,1),b[:,1].reshape(-1,1),b[:,2].reshape(-1,1)),axis=0).reshape(-1,1)
for i in tqdm(count):
    if i==count[0]:
        C11=c11(phi[i:i+step],theta[i:i+step],r[i:i+step],mesh)
        C12=c12(phi[i:i+step],theta[i:i+step],r[i:i+step],mesh)
        C13=c13(phi[i:i+step],theta[i:i+step],r[i:i+step],mesh)
        X1=cp.column_stack((C11,C12,C13))
        del C11,C12,C13
        Gm1=cp.dot(X1,x0).reshape(-1,1)
        del X1
        C21=c21(phi[i:i+step],theta[i:i+step],r[i:i+step],mesh)
        C22=c22(phi[i:i+step],theta[i:i+step],r[i:i+step],mesh)
        C23=c23(phi[i:i+step],theta[i:i+step],r[i:i+step],mesh)
        X2=cp.column_stack((C21,C22,C23))
        del C21,C22,C23
        Gm2=cp.dot(X2,x0).reshape(-1,1)
        del X2
        C31=c31(phi[i:i+step],theta[i:i+step],r[i:i+step],mesh)
        C32=c32(phi[i:i+step],theta[i:i+step],r[i:i+step],mesh)
        C33=c33(phi[i:i+step],theta[i:i+step],r[i:i+step],mesh)
        X3=cp.column_stack((C31,C32,C33))
        del C31,C32,C33
        Gm3=cp.dot(X3,x0).reshape(-1,1)
        Gm_temp=cp.column_stack((Gm1,Gm2,Gm3))
        Gm=Gm_temp
    elif i==count[-1]:
        C11=c11(phi[i:],theta[i:],r[i:],mesh)
        C12=c12(phi[i:],theta[i:],r[i:],mesh)
        C13=c13(phi[i:],theta[i:],r[i:],mesh)
        X1=cp.column_stack((C11,C12,C13))
        del C11,C12,C13
        Gm1=cp.dot(X1,x0).reshape(-1,1)
        del X1
        C21=c21(phi[i:],theta[i:],r[i:],mesh)
        C22=c22(phi[i:],theta[i:],r[i:],mesh)
        C23=c23(phi[i:],theta[i:],r[i:],mesh)
        X2=cp.column_stack((C21,C22,C23))
        del C21,C22,C23
        Gm2=cp.dot(X2,x0).reshape(-1,1)
        del X2
        C31=c31(phi[i:],theta[i:],r[i:],mesh)
        C32=c32(phi[i:],theta[i:],r[i:],mesh)
        C33=c33(phi[i:],theta[i:],r[i:],mesh)
        X3=cp.column_stack((C31,C32,C33))
        del C31,C32,C33
        Gm3=cp.dot(X3,x0).reshape(-1,1)
        del X3
        Gm_temp=cp.column_stack((Gm1,Gm2,Gm3))
        Gm=cp.vstack((Gm,Gm_temp))
    else:
        C11=c11(phi[i:i+step],theta[i:i+step],r[i:i+step],mesh)
        C12=c12(phi[i:i+step],theta[i:i+step],r[i:i+step],mesh)
        C13=c13(phi[i:i+step],theta[i:i+step],r[i:i+step],mesh)
        X1=cp.column_stack((C11,C12,C13))
        del C11,C12,C13
        Gm1=cp.dot(X1,x0).reshape(-1,1)
        del X1
        C21=c21(phi[i:i+step],theta[i:i+step],r[i:i+step],mesh)
        C22=c22(phi[i:i+step],theta[i:i+step],r[i:i+step],mesh)
        C23=c23(phi[i:i+step],theta[i:i+step],r[i:i+step],mesh)
        X2=cp.column_stack((C21,C22,C23))
        del C21,C22,C23
        Gm2=cp.dot(X2,x0).reshape(-1,1)
        del X2
        C31=c31(phi[i:i+step],theta[i:i+step],r[i:i+step],mesh)
        C32=c32(phi[i:i+step],theta[i:i+step],r[i:i+step],mesh)
        C33=c33(phi[i:i+step],theta[i:i+step],r[i:i+step],mesh)
        X3=cp.column_stack((C31,C32,C33))
        del C31,C32,C33
        Gm3=cp.dot(X3,x0).reshape(-1,1)
        Gm_temp=cp.column_stack((Gm1,Gm2,Gm3))
        Gm=cp.vstack((Gm,Gm_temp))
Gm=Gm.T.reshape(-1,1)
rr0=Gm - b;
#delta=15
for i in range(iter_max2):
    w=abs(rr0)
    w=w.reshape(3,nr).T
    w[0:229038][w[0:229038]<delta_MGS]=delta_MGS
    w[229038:229038+869868][w[229038:229038+869868]<delta_MAVEN]=delta_MAVEN
    w[229038+869868:][w[229038+869868:]<delta_tianwen]=delta_tianwen
    rr0=rr0.reshape(3,nr).T
    del G
    G=cp.zeros([n,n])
    r0=cp.zeros([n,1])
    for i in tqdm(count):
        if i==count[-1]:
            C11=c11(phi[i:],theta[i:],r[i:],mesh)
            C12=c12(phi[i:],theta[i:],r[i:],mesh)
            C13=c13(phi[i:],theta[i:],r[i:],mesh)
            X1=cp.column_stack((C11,C12,C13))
            del C11,C12,C13
            r0=r0+cp.dot(cp.dot(X1.T,cp.diag((delta[i:]/w[i:,0].flatten()))),rr0[i:,0].reshape(-1,1))
            G=G+cp.dot(cp.dot(X1.T,cp.diag((delta[i:]/w[i:,0].flatten()))),X1)
            del X1
            C21=c21(phi[i:],theta[i:],r[i:],mesh)
            C22=c22(phi[i:],theta[i:],r[i:],mesh)
            C23=c23(phi[i:],theta[i:],r[i:],mesh)
            X2=cp.column_stack((C21,C22,C23))
            del C21,C22,C23
            r0=r0+cp.dot(cp.dot(X2.T,cp.diag((delta[i:]/w[i:,1].flatten()))),rr0[i:,1].reshape(-1,1))
            G=G+cp.dot(cp.dot(X2.T,cp.diag((delta[i:]/w[i:,1].flatten()))),X2)
            del X2
            C31=c31(phi[i:],theta[i:],r[i:],mesh)
            C32=c32(phi[i:],theta[i:],r[i:],mesh)
            C33=c33(phi[i:],theta[i:],r[i:],mesh)
            X3=cp.column_stack((C31,C32,C33))
            del C31,C32,C33
            r0=r0+cp.dot(cp.dot(X3.T,cp.diag((delta[i:]/w[i:,2].flatten()))),rr0[i:,2].reshape(-1,1))
            G=G+cp.dot(cp.dot(X3.T,cp.diag((delta[i:]/w[i:,2].flatten()))),X3)
            del X3
        else:
            C11=c11(phi[i:i+step],theta[i:i+step],r[i:i+step],mesh)
            C12=c12(phi[i:i+step],theta[i:i+step],r[i:i+step],mesh)
            C13=c13(phi[i:i+step],theta[i:i+step],r[i:i+step],mesh)
            X1=cp.column_stack((C11,C12,C13))
            del C11,C12,C13
            r0=r0+cp.dot(cp.dot(X1.T,cp.diag((delta[i:i+step]/w[i:i+step,0].flatten()))),rr0[i:i+step,0].reshape(-1,1))
            G=G+cp.dot(cp.dot(X1.T,cp.diag((delta[i:i+step]/w[i:i+step,0].flatten()))),X1)
            del X1
            C21=c21(phi[i:i+step],theta[i:i+step],r[i:i+step],mesh)
            C22=c22(phi[i:i+step],theta[i:i+step],r[i:i+step],mesh)
            C23=c23(phi[i:i+step],theta[i:i+step],r[i:i+step],mesh)
            X2=cp.column_stack((C21,C22,C23))
            del C21,C22,C23
            r0=r0+cp.dot(cp.dot(X2.T,cp.diag((delta[i:i+step]/w[i:i+step,1].flatten()))),rr0[i:i+step,1].reshape(-1,1))
            G=G+cp.dot(cp.dot(X2.T,cp.diag((delta[i:i+step]/w[i:i+step,1].flatten()))),X2)
            del X2
            C31=c31(phi[i:i+step],theta[i:i+step],r[i:i+step],mesh)
            C32=c32(phi[i:i+step],theta[i:i+step],r[i:i+step],mesh)
            C33=c33(phi[i:i+step],theta[i:i+step],r[i:i+step],mesh)
            X3=cp.column_stack((C31,C32,C33))
            del C31,C32,C33
            r0=r0+cp.dot(cp.dot(X3.T,cp.diag((delta[i:i+step]/w[i:i+step,2].flatten()))),rr0[i:i+step,2].reshape(-1,1))
            G=G+cp.dot(cp.dot(X3.T,cp.diag((delta[i:i+step]/w[i:i+step,2].flatten()))),X3)
    alpha = cp.dot(r0.T,r0)/cp.dot(cp.dot(p0.T,G),p0);
    x = x0 + alpha*p0;
    r1 = r0 + alpha*cp.dot(G,p0);
    beta = cp.dot(r1.T,r1)/cp.dot(r0.T,r0);
    p = -r1 + beta*p0;
    misfit.append(cp.linalg.norm(x=r0,ord=2))
    if cp.asnumpy(cp.linalg.norm(x=r1-r0,ord=2)/(cp.linalg.norm(x=r0,ord=2)))<=1e-2:
        break
    number=number+1
    print(number)
    x0 = x;
    for i in tqdm(count):
        if i==count[0]:
            C11=c11(phi[i:i+step],theta[i:i+step],r[i:i+step],mesh)
            C12=c12(phi[i:i+step],theta[i:i+step],r[i:i+step],mesh)
            C13=c13(phi[i:i+step],theta[i:i+step],r[i:i+step],mesh)
            X1=cp.column_stack((C11,C12,C13))
            del C11,C12,C13
            Gm1=cp.dot(X1,x0).reshape(-1,1)
            del X1
            C21=c21(phi[i:i+step],theta[i:i+step],r[i:i+step],mesh)
            C22=c22(phi[i:i+step],theta[i:i+step],r[i:i+step],mesh)
            C23=c23(phi[i:i+step],theta[i:i+step],r[i:i+step],mesh)
            X2=cp.column_stack((C21,C22,C23))
            del C21,C22,C23
            Gm2=cp.dot(X2,x0).reshape(-1,1)
            del X2
            C31=c31(phi[i:i+step],theta[i:i+step],r[i:i+step],mesh)
            C32=c32(phi[i:i+step],theta[i:i+step],r[i:i+step],mesh)
            C33=c33(phi[i:i+step],theta[i:i+step],r[i:i+step],mesh)
            X3=cp.column_stack((C31,C32,C33))
            del C31,C32,C33
            Gm3=cp.dot(X3,x0).reshape(-1,1)
            Gm_temp=cp.column_stack((Gm1,Gm2,Gm3))
            Gm=Gm_temp
        elif i==count[-1]:
            C11=c11(phi[i:],theta[i:],r[i:],mesh)
            C12=c12(phi[i:],theta[i:],r[i:],mesh)
            C13=c13(phi[i:],theta[i:],r[i:],mesh)
            X1=cp.column_stack((C11,C12,C13))
            del C11,C12,C13
            Gm1=cp.dot(X1,x0).reshape(-1,1)
            del X1
            C21=c21(phi[i:],theta[i:],r[i:],mesh)
            C22=c22(phi[i:],theta[i:],r[i:],mesh)
            C23=c23(phi[i:],theta[i:],r[i:],mesh)
            X2=cp.column_stack((C21,C22,C23))
            del C21,C22,C23
            Gm2=cp.dot(X2,x0).reshape(-1,1)
            del X2
            C31=c31(phi[i:],theta[i:],r[i:],mesh)
            C32=c32(phi[i:],theta[i:],r[i:],mesh)
            C33=c33(phi[i:],theta[i:],r[i:],mesh)
            X3=cp.column_stack((C31,C32,C33))
            del C31,C32,C33
            Gm3=cp.dot(X3,x0).reshape(-1,1)
            del X3
            Gm_temp=cp.column_stack((Gm1,Gm2,Gm3))
            Gm=cp.vstack((Gm,Gm_temp))
        else:
            C11=c11(phi[i:i+step],theta[i:i+step],r[i:i+step],mesh)
            C12=c12(phi[i:i+step],theta[i:i+step],r[i:i+step],mesh)
            C13=c13(phi[i:i+step],theta[i:i+step],r[i:i+step],mesh)
            X1=cp.column_stack((C11,C12,C13))
            del C11,C12,C13
            Gm1=cp.dot(X1,x0).reshape(-1,1)
            del X1
            C21=c21(phi[i:i+step],theta[i:i+step],r[i:i+step],mesh)
            C22=c22(phi[i:i+step],theta[i:i+step],r[i:i+step],mesh)
            C23=c23(phi[i:i+step],theta[i:i+step],r[i:i+step],mesh)
            X2=cp.column_stack((C21,C22,C23))
            del C21,C22,C23
            Gm2=cp.dot(X2,x0).reshape(-1,1)
            del X2
            C31=c31(phi[i:i+step],theta[i:i+step],r[i:i+step],mesh)
            C32=c32(phi[i:i+step],theta[i:i+step],r[i:i+step],mesh)
            C33=c33(phi[i:i+step],theta[i:i+step],r[i:i+step],mesh)
            X3=cp.column_stack((C31,C32,C33))
            del C31,C32,C33
            Gm3=cp.dot(X3,x0).reshape(-1,1)
            Gm_temp=cp.column_stack((Gm1,Gm2,Gm3))
            Gm=cp.vstack((Gm,Gm_temp))
    Gm=Gm.T.reshape(-1,1)
    rr0=Gm - b;
    p0 = p;
m0=x0