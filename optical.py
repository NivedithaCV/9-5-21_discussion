import pandas as pd
import numpy as np
import math
import random
from pandas import ExcelFile
from pandas import ExcelWriter
import matplotlib.pyplot as plt
from scipy import interpolate
# endmember data
xls= pd.ExcelFile('olivine.xlsx')
df1=pd.read_excel(xls,'olivine.txt', usecols=[0,1,2], names=['colA','colB','colC'],header=None)
lmda = df1['colA'].tolist()
lmda=np.asarray(lmda)
n = df1['colB'].tolist()
n=np.asarray(n)
k=df1['colC'].tolist()
k=np.asarray(k)

xl_= pd.ExcelFile('O33E33A33.xlsx')
df_=pd.read_excel(xl_,'O33E33A33.txt', usecols=[0,1], names=['colA','colB'],header=None)
lmda_ = df_['colA'].tolist()
lmda_=np.asarray(lmda_)
n_ = df_['colB'].tolist()
n_=np.asarray(n_)

df2=pd.read_excel(xls,'anorthite.txt', usecols=[0,1,2], names=['colA','colB','colC'],header=None)
lmda2 = df2['colA'].tolist()
lmda2=np.asarray(lmda2)
n2 = df2['colB'].tolist()
n2=np.asarray(n2)
k2=df2['colC'].tolist()
k2=np.asarray(k2)

df3=pd.read_excel(xls,'enstatite.txt', usecols=[0,1,2], names=['colk','colr','cols'],header=None)
lmda3 = df3['colk'].tolist()
lmda3=np.asarray(lmda3)
n3 = df3['colr'].tolist()
n3=np.asarray(n3)
k3=df3['cols'].tolist()
k3=np.asarray(k3)

#__________________________________________________________________________________________________________________________________________________

rho=[3.32,3.20,2.73]
m=[0.33,0.33,0.33]
D=random.sample(range(45,75),3)
D=np.asarray(D)*10**(-6)
print(D)

#mass abundance , density, grain size of end member
def fractional_relative_cross_section(m,rho,D,i):
    k=0; sum=0
    sig_i=m[i]/(rho[i]*D[i])
    for j in rho:
        sig_k=m[k]/(rho[k]*D[k])
        sum=sum+sig_k
        k=k+1
    f_i=sig_i/sum
    return(f_i)

    sig_i=m[i]/(rho[i]*D[i])

#def Hapke_forwar(fi,Di,n,k,lambda):pass
def sswa_endmember(n,k,lamda,D,i):
    alpha=(4*(math.pi)*k)/lamda

    s=0
    r_i=(1-(np.sqrt(alpha/(alpha+s)) ))/(1+(np.sqrt(alpha/(alpha+s)) ))
    Da=(2/3)*D[i]*(np.power(n,2)-((1/n)*np.power((np.power(n,2)-1),1.5)))
    theta=(r_i+np.exp(-np.sqrt(alpha*(alpha+s)*Da)))/(1+r_i*np.exp(-np.sqrt(alpha*(alpha+s)*Da)))
    Se=((np.power((n-1),2)+np.power(k,2))/(np.power((n+1),2)+np.power(k,2)) ) +0.05
    Si=1.014-( 4/n*np.power((n+1),2) )
    w=Se+(1-Se)*((1-Si)/(1-Si*theta))*theta
    return(w)
w_0=sswa_endmember(n,k,lmda,D,0)
w_1=sswa_endmember(n2,k2,lmda2,D,1)
w_2=sswa_endmember(n3,k3,lmda3,D,2)


f_0=fractional_relative_cross_section(m,rho,D,0)
f_1=fractional_relative_cross_section(m,rho,D,1)
f_2=fractional_relative_cross_section(m,rho,D,2)

# interpolation of SSA
minim=[np.amin(lmda),np.amin(lmda2),np.amin(lmda3)]
minspectra=[np.amax(minim)]
print(minspectra)
maxim=[np.amax(lmda),np.amax(lmda2),np.amax(lmda3)]
maxspectra=[np.amin(maxim)]
print(maxspectra)

dif=np.diff(lmda)
diff=np.mean(dif)/2
print(diff)
lam=np.arange(minspectra[0],maxspectra[0],diff)
f0=interpolate.interp1d(lmda,w_0,kind='cubic')
f1=interpolate.interp1d(lmda2,w_1,kind='cubic')
f2=interpolate.interp1d(lmda3,w_2,kind='linear')
def interpolation(lmda,w_0,lmda2,w_1,lmd3,w_2):pass
w__0=f0(lam)
w__1=f1(lam)
w__2=f2(lam)

w_mix=(f_1*w__1)+(f_0*w__0)+(f_2*w__2)
#SSA to reflectance
def Chandrashekar(w,x):
    gamma =np.sqrt(1-w)
    r_0=(1-gamma)/(1+gamma)
    H= 1/(1-w*x*(r_0+((1-2*r_0*x)/2)) *(math.log((1+x)/x)))
    return(H)
a=math.pi/6
H_mu=Chandrashekar(w_mix,math.cos(a))
H_m0=Chandrashekar(w_mix,math.cos(0))

r=(w_mix*H_mu*H_m0)/(4*(math.cos(a)+math.cos(0)))

plt.plot(lam,w_mix)
plt.xlabel("wavelength")
plt.ylabel("SSA of mix")
plt.title("Olivine33% and anorthite33% and enstatite33%")
plt.show()



plt.plot(lam,r)
plt.xlabel("wavelength")
plt.ylabel("reflectance of mix")
plt.title("Olivine33% and anorthite33% and enstatite33%")
plt.show()
