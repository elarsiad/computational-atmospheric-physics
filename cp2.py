#Computational Problem 2
import math
import numpy as np
import matplotlib.pyplot as plt

## Read the atmosphere layers effective ztp from problem 1 of CP1
f = open('data\ztp_eff.txt','r')
lines = f.readlines()

z = []
t = []
p = []

nlayer = 0
for line in lines:
    z.append(float(line[0:4]))
    t.append(float(line[5:11]))
    p.append(float(line[12:23]))
    # print(all_z_eff[count],all_t_eff[count],all_p_eff[count])
    nlayer += 1
# print(type(all_z_eff),type(all_t_eff),type(all_p_eff))
f.close()

SZAs = range(0,81,1)
mr = [[0 for i in range(nlayer)],[0 for i in range(nlayer)]] # shape: 2 x 10
incang = [[[0 for i in SZAs] for i in range(nlayer)],[[0 for i in SZAs] for i in range(nlayer)]] # shape: 2 x 10 x 81

wls = [350e-9,1000e-9] # wavelength

## PROBLEM 1 start computation
PL = [[[0 for i in SZAs] for i in range(nlayer)],[[0 for i in SZAs] for i in range(nlayer)]] # shape: 2 x 10 x 81
PLsum = [[0 for i in SZAs],[0 for i in SZAs]] # shape: 2 x 81

for j in range(2):
    wl = wls[j] * 1e6 # convert wavelength to micron
    for SZA in SZAs:
        for i in reversed(range(nlayer)):
            mr[j][i] = (6432.8 + 2949810/(146-wl**(-2)) + 25540/(41-wl**(-2))) / (1013.25/288.15) * (p[i]/t[i]) * 1e-8 + 1
            if i == 9:
                incang[j][i][SZA] = SZA
            else:
                incang[j][i][SZA] = np.arcsin(mr[j][i+1]*math.sin(incang[j][i+1][SZA]/180*math.pi)/mr[j][i]) *180/math.pi
            PL[j][i][SZA] = 5/math.cos(incang[j][i][SZA]/180*math.pi)
            PLsum[j][SZA] = PLsum[j][SZA] + PL[j][i][SZA]
        PLsum[j][SZA] = 50/math.cos(SZA/180*math.pi) - PLsum[j][SZA]

# plot mr (refractive index)
plt.figure(4)
plt.plot(mr[0][:],z)
plt.plot(mr[1][:],z,'--')
plt.legend(['350 nm','1000 nm'])
plt.title('Refractive Index vs height')
plt.ylabel('Height [km]')
plt.xlabel('mr')
plt.show()

# plot Path Length difference from refraction
plt.figure(1)
plt.plot(PLsum[0])
plt.plot(PLsum[1],'--')
plt.title('Refracted path length difference')
plt.ylabel('Path Length [km]')
plt.xlabel('SZA')
plt.legend(['350 nm','1000 nm'])
plt.show()

print(PLsum[0][60])
print(PLsum[1][60])
print(PLsum[1][80]-PLsum[0][80])

## PROBLEM 2 start computation
UWbot = [[0 for i in SZAs],[0 for i in SZAs]] # shape: 2 x 81
UWsum = [[0 for i in SZAs],[0 for i in SZAs]] # shape: 2 x 81
UW = [[[0 for i in SZAs] for i in range(nlayer)],[[0 for i in SZAs] for i in range(nlayer)]] # shape: 2 x 10 x 81
DW = [[[0 for i in SZAs] for i in range(nlayer)],[[0 for i in SZAs] for i in range(nlayer)]] # shape: 2 x 10 x 81

Av = 6.02295 * 1e23 #(/mol) Avogadro number
M = 28.97 #(g/mol) molecular weight of air
delt = 0.035
fdelt = (6+3*delt)/(6-7*delt)
A = 1.0 #Albedo

for j in range(2):
    wl = wls[j] * 1e2 # convert wavelength to cm
    for SZA in SZAs:
        for i in reversed(range(nlayer)):
            rho = p[i]/t[i]/(287.058*10)
            Ns = rho*Av/M
            sigs = 8*math.pi**3*(mr[j][i]**2-1)**2/(3*wl**4*Ns**2)*fdelt
            if i == 9:
                Iin = 1*math.cos(SZA/180*math.pi)
            else:
                Iin = DW[j][i+1][SZA]*math.cos(incang[j][i+1][SZA]/180*math.pi)/math.cos(incang[j][i][SZA]/180*math.pi)
            #scattered light
            Iscat = Iin*(1-math.exp(-sigs*Ns*PL[j][i][SZA]*1e5))
            #source term SZA of the scattered light as upwelling
            UW[j][i][SZA] = Iscat*3/4*(1+math.cos((180-incang[j][i][SZA])/180*math.pi)**2)/(4*math.pi)
            UWbot[j][SZA] = UWbot[j][SZA] + Iscat*3/4*(1+math.cos((180-incang[j][i][SZA])/180*math.pi)**2)/(4*math.pi)
            #Incoming - scattered + forward scattered as downwelling
            DW[j][i][SZA] = Iin - Iscat + Iscat*3/4*(1+math.cos(0/180*math.pi)**2)/(4*math.pi)
        UWbot[j][SZA] = UWbot[j][SZA] + DW[j][i][SZA]*A/math.pi

for j in range(2):
    wl = wls[j] * 1e2 # convert wavelength to cm
    for SZA in SZAs:
        UWsum[j][SZA] = UWbot[j][SZA]
        for i in range(nlayer):
            rho = p[i]/t[i]/(287.058*10)
            Ns = rho*Av/M
            sigs = 8*math.pi**3*(mr[j][i]**2-1)**2/(3*wl**4*Ns**2)*fdelt
            Iscat = UWsum[j][SZA]*(1-math.exp(-sigs*Ns*PL[j][i][SZA]*1e5))
            UWsum[j][SZA] = UWsum[j][SZA] - Iscat + UW[j][i][SZA]

# plot upwelling radiation to SZA
plt.figure(2)
plt.plot(UWsum[0])
plt.plot(UWsum[1],'--')
plt.title('upwelling radiation')
plt.ylabel('Fraction of I0')
plt.xlabel('SZA')
plt.axis([0,80,0,0.35])
plt.legend(['350 nm','1000 nm'])
plt.show()

# plot downwelling radiation to SZA
plt.figure(3)
plt.title('downwelling radiation per layer')
for i in reversed(range(nlayer)):
    plt.subplot2grid((10,2),((9-i)%5*2,(9-i)//5))
    plt.plot(DW[0][i])
    plt.plot(DW[1][i],'--')
    plt.title('layer '+str(int(z[i]))+' km', fontsize=10)
    #plt.ylabel('Fraction')
    plt.xlabel('SZA',fontsize=8)
    plt.xticks(fontsize=8)
    plt.yticks(fontsize=8)
    plt.axis([0,80,0,1.1])
    if i == 5:
        plt.legend(['350 nm','1000 nm'])
plt.show()

# plot downwelling radiation to SZA in contour
plt.figure(5)
plt.title('Downwelling 350 nm')
plt.contourf(SZAs,z,DW[0],levels=[0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1], ticks=[0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1],
colors=[[0,0,0],[0.2,0,0],[0.4,0,0],[0.6,0,0],[0.8,0,0],[1,0,0],[1,0.2,0.2],[1,0.4,0.4],[1,0.6,0.6],[1,0.8,0.8],[1,1,1]])
plt.colorbar()
plt.ylabel('Altitude [km]')
plt.xlabel('SZA')
plt.show()
plt.figure(6)
plt.title('Downwelling 1000 nm')
plt.contourf(SZAs,z,DW[1],levels=[0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1], ticks=[0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1],
colors=[[0,0,0],[0.2,0,0],[0.4,0,0],[0.6,0,0],[0.8,0,0],[1,0,0],[1,0.2,0.2],[1,0.4,0.4],[1,0.6,0.6],[1,0.8,0.8],[1,1,1]])
plt.colorbar()
plt.ylabel('Altitude [km]')
plt.xlabel('SZA')
plt.show()