#!/usr/bin/env python
# coding: utf-8

# In[1]:


import sys

import math

import matplotlib.pyplot as plt

import numpy as np
# Version 5
#read data: z, sample, conodont zone, d13Ccon, taxon, d13Ccarb, d18Ocarb


fname=input("Datafile=")

z=[]

sample=[]

zon=[]

dccon=[]

tax=[]

dccarb=[]

docarb=[]

ni=0

udcd=[]

#open file fname (comma delimited)

input1 = open(fname, "r")

text=input1.readlines()

input1.close()

i=0

for line in text:

    items = line.strip("()\n").split(",")

    z.append(float(items[0]))

    sample.append((items[1]))

    zon.append(items[2])

    dccon.append(float(items[3]))

    tax.append(items[4])

    dccarb.append(float(items[5]))

    docarb.append(float(items[6]))

    udcd.append(float(items[5])-float(items[3]))

    ni=ni+1


#number of records is ni

#z - position of the sample in the succession

#sample - sample name

#zon - conodont Zone

#dccon - d13Ccon value (PDB)

#tax - conodont taxon

#dccarb - d13Ccarb (PDB) 

#docarb - d18Ocarb (SMOW)

#udcd - difference between d13Ccarb and d13Ccon - delta 13C(car-con)

print(ni)

umic=[]

umac=[]

umim=[]

umam=[]

umm=[]

umad=[]

umid=[]


#Input of d18O of sea water  PDB for water temperature
do1=input("d18Oseawater PDB (about -1.2)=")
dow=float(do1)
t=[]

print("Max DELTA="+str(max(udcd)))

# Input probable Ef
ef=float(input("Ef(20-29)="))

for e in range(ni):
#SMOW--PDB
    dop=(docarb[e]-30.91)/1.03092
# temperature calculation
    t.append(15.7-4.36*(dop-dow)+0.12*(dop-dow)*(dop-dow))
    
    dcd=udcd[e]

# mmu - growth rate of phytoplankton

    mmu=0.0

# delt - computed value of delta 13C(car-con)

    delt=0.0

# co - CO2 content in water

    co=0.0

# dcph - d13C of phytoplankton

    dcph=0.0

# ddic - d13C of DIC 

    ddic=0.0

    m=[]

    d=[]

    c=[]

    msel=[]

    csel=[]

    dcsel=[]

    for i  in range(200, 900, 10):

        for j in range(1,170,1):

            mmu=j/100.0

            m.append(mmu)

# According to Henry's Law CO2aq=CO2at/29.41

            co=float(i)/29.41

            c.append(co)

# ddic=dccarb+1.2 (Hartke et al., 2021)

# dco=ddic-(-0.095t + 10.72) (Yoshioka, 1997)

# Ef=28

# Ep=Ef-182.0*(mmu/co) (Hartke et al., 2021)

# dcph=dco-Ep

# dccon(calculated)= dcph-0.113*dcph-1.916=0.887*dcph-1.916 (Caut et al., 2009) invertebrates
#
# dccon(calculated)= 0.787*dcph-2.848 (Caut et al., 2009) vertebrates

# delt=dccar-dccon(calculated)=dccar-(0.887*dcph-1.916)
#

            ddic=dccarb[e]+1.2

            dco=ddic+0.095*t[e]-10.72

            Ep=ef-182.0*(mmu/co)

            dcph=dco-Ep

            delt=dccarb[e]-0.787*dcph+2.848

            d.append(delt)

#possible individual variations are assumed as about +-0.5 promille

            if abs(delt-dcd)<0.5:

                csel.append(co)

                msel.append(mmu)

                dcsel.append(dcph)

    mico=min(csel)

    mimu=min(msel)

    maco=max(csel)

    mamu=max(msel)

    midc=min(dcsel)

    madc=max(dcsel)
     
    memu=sum(msel)/len(msel)

#    print("PhPLK d13C from "+str(midc)+" to "+str(madc))
#    print("PHPLK growth from "+str(mimu)+" to "+str(mamu))

    umic.append(mico)

    umac.append(maco)

    umim.append(mimu)

    umam.append(mamu)

#    umm.append((mamu+mimu)/2.0)
    umm.append(memu)

    umad.append(madc)

    umid.append(midc)


adcd=sum(udcd)/ni
tt=sum(t)/ni
tt=int(tt)

fig, ax = plt.subplots()

ax.set_xlabel("z")

ax.set_ylabel("Delta d13C(car-con)")

ax.plot(z,udcd, 'r*')

ax.axhline(adcd, ls='-', color='g')

fig.tight_layout()

plt.show()


figa, axx = plt.subplots()

axx.set_xlabel("z")

axx.set_ylabel("PHPLK growth rate (min - max); Mean t="+str(tt)+" Ef="+str(ef))

axx.plot(z,umim, 'g.')

axx.plot(z,umam,'r.')

#axx.plot(z,umm,'b*')

figa.tight_layout()

plt.show()



figa, a = plt.subplots()

a.set_xlabel("z")

a.set_ylabel("Mean PHPLK growth rate;  Mean t="+str(tt)+" Ef="+str(ef))

a.plot(z,umm,'b*')

figa.tight_layout()

plt.show()



figb, ay = plt.subplots()

ay.set_xlabel("z")

ay.set_ylabel("CO2")

ay.plot(z,umic, 'g.')

ay.plot(z,umac, 'g.')

figb.tight_layout()

plt.show()


figt, axt = plt.subplots()

axt.set_xlabel("z")

axt.set_ylabel("t; Mean t="+str(tt)+" Ef="+str(ef))

axt.plot(z,t, 'g.')

figt.tight_layout()

plt.show()


figc, aaa=plt.subplots()

aaa.set_ylabel("t")

aaa.set_xlabel("PHPLK growth rate (mean); Ef="+str(ef))

aaa.plot(umm,t, 'g.')

p1, p0 = np.polyfit(umm,t, deg=1)  # slope, intercept

aaa.axline(xy1=(0, p0), slope=p1, color='r', lw=2)

figb.tight_layout()

plt.show()

#correlation matrix dccon, dccarb, docarb, umm

arr=np.array([dccon, dccarb, docarb, umm])

r = np.corrcoef(arr).round(decimals=2)

figd, ag = plt.subplots()

im = ag.imshow(r)

im.set_clim(-1, 1)

ag.grid(False)

ag.xaxis.set(ticks=(0, 1, 2, 3), ticklabels=('d13Ccon', 'd13Ccar', 'd18Ocar', "PLKgrowth"))

ag.yaxis.set(ticks=(0, 1, 2, 3), ticklabels=('d13Ccon', 'd13Ccar', 'd18Ocar', "PLKgrowth"))



for i in range(4):

    for j in range(4):

        ag.text(j, i, r[i, j], ha='center', va='center', color='r')


cbar = ag.figure.colorbar(im, ax=ag, format='% .2f')

plt.show()


#Output - save results in csv file

oname=input("Output Datafile=")

output2 = open(oname, "w")

text=" "

text="Position; Zone; Sample; Taxon; d13Ccon; d13Ccarb; d18Ocarb; delta d13C; PhPLK growth min; PhPLK growth max; PhPLK growth mean; d13Cph min; d13Cph max \n"

output2.write(text) 

for j in range(ni):

    text=(str(z[j])+";"+zon[j]+";"+sample[j]+";"+tax[j]+";"  +str(dccon[j])+";"+str(dccarb[j])+";"+str(docarb[j])+";"+str(udcd[j])+";"+str(umim[j])+";"+str(umam[j])+";"+str(umm[j])+";"+str(umid[j])+";"+str(umad[j])+"\n")
    output2.write(text)

output2.close()



# In[3]:


import sys
import math
import random
import matplotlib.pyplot as plt
import numpy as np
#read data filename 3a.txt
fname="3a.txt"

sample=[]
co=[]
dccarb=[]
docarb=[]
dccon=[]
deltac=[]
#open file fname (comma delimited)
input1 = open(fname, "r")
text=input1.readlines()
input1.close()
i=0
for line in text:
    items = line.strip("()\n").split(",")
    sample.append((items[0]))
    co.append(float(items[1])/29.41)
    dccarb.append(float(items[2]))
    docarb.append(float(items[3]))
    dccon.append(float(items[4]))
    i=i+1
print(i)

#initial calculation
for j in range(i):
    
        deltac.append(dccarb[j]-dccon[j])

nd=i
md=0.0
md=(sum(deltac)/i)
print("Mean DELTA")
print(md)


fig, ax = plt.subplots()
ax.set_xlabel("Delta 13C(car-con)")
ax.set_ylabel("d13Ccarb")
ax.plot(deltac,dccarb, 'g*')
fig.tight_layout()
plt.show()


fig, ay=plt.subplots()
ay.set_xlabel("Delta 13C(car-con)")
ay.set_ylabel("d18Ocarb")
ay.plot(deltac,docarb, 'g*')
fig.tight_layout()
plt.show()


fig, az=plt.subplots()
az.set_xlabel("Delta 13C(car-con)")
az.set_ylabel("d13Ccon")
az.plot(deltac,dccon, 'g*')

p1, p0 = np.polyfit(deltac,dccon, deg=1)  # slope, intercept

az.axline(xy1=(0, p0), slope=p1, color='r', lw=2)
fig.tight_layout()
plt.show()

nb=0.0
nb=int(2.0*(nd**0.3333))
fig, axx= plt.subplots()
n, bins, patches = axx.hist(deltac, nb, density=True)
axx.set_xlabel("Delta 13C(car-con)  " + "   N="+str(nd))
fig.tight_layout()
plt.show()



#correlation matrix dccon, dccarb, docarb, deltac

arr=np.array([dccon, dccarb, docarb, deltac])

r = np.corrcoef(arr).round(decimals=2)

figd, ag = plt.subplots()

im = ag.imshow(r)

im.set_clim(-1, 1)

ag.grid(False)

ag.xaxis.set(ticks=(0, 1, 2, 3), ticklabels=('d13Ccon', 'd13Ccar', 'd18Ocar', "Delta C"))

ag.yaxis.set(ticks=(0, 1, 2, 3), ticklabels=('d13Ccon', 'd13Ccar', 'd18Ocar', "Delta C"))



for i in range(4):

    for j in range(4):

        ag.text(j, i, r[i, j], ha='center', va='center', color='r')


cbar = ag.figure.colorbar(im, ax=ag, format='% .2f')

plt.show()


# In[ ]:




