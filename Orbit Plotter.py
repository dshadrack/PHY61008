# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import numpy as np
import csv
AU = 1.5e11
file = open('F:/Fortran/Projects/3 Body Problem/Earth Motion.csv')
xe = np.array([])
ye = xe
ze = xe
t = xe
se = xe
xs = xe
ys = xe
zs = xe
ss = xe
xj = xe
yj = xe
zj = xe
sj = xe
xm = xe
ym = xe
zm = xe
sm = xe
xst = xe
yst = xe
zst = xe
sst = xe
xu = xe
yu = xe
zu = xe
su = xe
xn = xe
yn = xe
zn = xe
sn = xe

sstar = xe
#y = x**2
look = csv.reader(file,delimiter = ',')
for row in look:
    t = np.append(t,float(row[0]))
    xe = np.append(xe,float(row[1]))
    ye = np.append(ye,float(row[2]))
    ze = np.append(ze,float(row[3]))
    se = np.append(se,float(row[4]))
file.close()

file = open('F:/Fortran/Projects/3 Body Problem/Sun Motion.csv')
look = csv.reader(file,delimiter = ',')
for row in look:
    xs = np.append(xs,float(row[1]))
    ys = np.append(ys,float(row[2]))
    zs = np.append(zs,float(row[3]))
    ss = np.append(ss,float(row[4]))
file.close()

file = open('F:/Fortran/Projects/3 Body Problem/Jupiter Motion.csv')
look = csv.reader(file,delimiter = ',')
for row in look:
    xj = np.append(xj,float(row[1]))
    yj = np.append(yj,float(row[2]))
    zj = np.append(zj,float(row[3]))
    sj = np.append(sj,float(row[4]))
file.close()

file = open('F:/Fortran/Projects/3 Body Problem/Mars Motion.csv')
look = csv.reader(file,delimiter = ',')
for row in look:
    xm = np.append(xm,float(row[1]))
    ym = np.append(ym,float(row[2]))
    zm = np.append(zm,float(row[3]))
    sm = np.append(sm,float(row[4]))
file.close()

file = open('F:/Fortran/Projects/3 Body Problem/Saturn Motion.csv')
look = csv.reader(file,delimiter = ',')
for row in look:
    xst = np.append(xst,float(row[1]))
    yst = np.append(yst,float(row[2]))
    zst = np.append(zst,float(row[3]))
    sst = np.append(sst,float(row[4]))
file.close()

file = open('F:/Fortran/Projects/3 Body Problem/Uranus Motion.csv')
look = csv.reader(file,delimiter = ',')
for row in look:
    xu = np.append(xu,float(row[1]))
    yu = np.append(yu,float(row[2]))
    zu = np.append(zu,float(row[3]))
    su = np.append(su,float(row[4]))
file.close()

file = open('F:/Fortran/Projects/3 Body Problem/Neptune Motion.csv')
look = csv.reader(file,delimiter = ',')
for row in look:
    xn = np.append(xn,float(row[1]))
    yn = np.append(yn,float(row[2]))
    zn = np.append(zn,float(row[3]))
    sn = np.append(sn,float(row[4]))
file.close()

file = open('F:/Fortran/Projects/3 Body Problem/star motion.csv')
look = csv.reader(file,delimiter = ',')
for row in look:
    sstar = np.append(sstar,float(row[1]))
file.close()







fig = plt.figure(figsize =(20,20))
axis = fig.add_subplot(111, projection = '3d',)
#plt.axis("off")
limit= 20
axis.set_xlim(-1*limit,limit)
axis.set_ylim(-1*limit,limit)
axis.set_zlim(-1*limit,limit)

plt.plot(xs/AU,ys/AU,zs/AU,c='yellow')
plt.plot(xj/AU,yj/AU,zj/AU,'-',c='orange')
plt.plot(xm/AU,ym/AU,zm/AU,'-',c='r')
plt.plot(xe/AU,ye/AU,ze/AU,'-',c='b')
plt.plot(xst/AU,yst/AU,zst/AU,'-',c='magenta')
plt.plot(xu/AU,yu/AU,zu/AU,'-',c='cyan')
plt.plot(xn/AU,yn/AU,zn/AU,'-',c='purple')
#plt.plot(2*x/AU,y/AU,z/AU, c = 'r')
#axis.view_init(0,-90)
plt.show()
#plt.plot(t,se/AU)
#print(t)
#print(x)
#print(row)
print(max(se/AU),min(se/AU))
print(max(sj/AU),min(sj/AU))
print(max(sm/AU),min(sm/AU))
print(max(sst/AU),min(sst/AU))
print(max(su/AU),min(su/AU))
print(max(sn/AU),min(sn/AU))
ee = (max(se/AU) - min(se/AU))/(max(se/AU) + min(se/AU))
ej = (max(sj/AU) - min(sj/AU))/(max(sj/AU) + min(sj/AU))
em = (max(sm/AU) - min(sm/AU))/(max(sm/AU) + min(sm/AU))
est = (max(sst/AU) - min(sst/AU))/(max(sst/AU) + min(sst/AU))
eu = (max(su/AU) - min(su/AU))/(max(su/AU) + min(su/AU))
en = (max(sn/AU) - min(sn/AU))/(max(sn/AU) + min(sn/AU))
Eccentricities = [ee,ej,em,est,eu,en]


i = 1
for e in Eccentricities:
    if e > 0.1:
        print("Perturbation Detected")
        print(i,e)
    i += 1
i = 1
print(min(sstar)/AU)
#%%
#fig,axis = plt.subplots(figsize = (20,20))
#plt.plot(t/(365*24),se/AU)
#axis.set_ylim(,0.9974)
#plt.plot(t/365,sj/AU)
#plt.show()