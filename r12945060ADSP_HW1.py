# -*- coding: utf-8 -*-
"""
Created on Wed Mar 20 18:21:12 2024

@author: cdpss
"""

import numpy as np
import matplotlib.pyplot as plt

##step one : choose initial extreme frequence
def choose_init_point(length,transitionband,fs):
    k = (length-1)//2
    init_point = np.linspace(0,(transitionband[0]-10)/fs,int((k+2)/2))
    init_point = np.append(init_point,np.linspace((transitionband[1]+10)/fs,0.5,int((k+2)/2)))
    return init_point,k

##step two : calculate the matrix
def cal_matrix(k,extremepoint):
    A = np.ones([k+2,k+1])
    W = np.ones([k+2,1])
    S = np.zeros_like(W)
    H = np.zeros_like(W)
    for i in range(k+2):
        for j in range(k):
            A[i,j+1] = np.cos(2*np.pi*(j+1)*extremepoint[i])
        if extremepoint[i] < transitionband[0]/fs:
            W[i] = 1 / weighted[0]*((-1)**i)
            H[i] = 1
        elif extremepoint[i] > transitionband[1]/fs:
            W[i] = 1 / weighted[1] * ((-1)**i)
            H[i] = 0
    AW = np.append(A,W, axis=1)
    S = np.dot(np.linalg.inv(AW) , H)
    return AW,H,S

#step three : calculate new error function and corresponding extremepoint
def cal_err(S,k):
    R = 0
    err = []
    for i in range(len(S)-1):
        R += S[i]*np.cos(2*np.pi*F*i)
    err = (R-Hdf)*Wf
    extreme_p, value = find_extremep(err, k)
    return err,extreme_p,value,R

def find_extremep(err,k):
    err = np.insert(err,0,0)
    err = np.append(err,0)
    value = []
    extreme_point = []
    for i in range(1,len(err)-1):
        if err[i]>err[i-1] and err[i]>err[i+1]:   #local max
            value.append(err[i])
            extreme_point.append(F[i-1])
        elif err[i]<err[i-1] and err[i]<err[i+1]: #local min
            value.append(err[i])
            extreme_point.append(F[i - 1])
    if len(extreme_point) > k+2:
        extreme_point.pop(-1)
        value.pop(-1)
    # print(len(extreme_point))
    return extreme_point,value

## constant parameter
length = 17
passband = [0,1200]
transitionband = [1200,1500]
fs = 6000
weighted = [1,0.6]
F = np.linspace(0, 0.6, 50001)
Hdf = np.ones_like(F)
Hdf[F > 1350/fs] = 0
Wf = np.zeros_like(F)
Wf[F < 1200/fs ] = 1
Wf[F > 1500/fs ] = 0.6
E0 = []

#step one
init_point, k = choose_init_point(length, transitionband, fs)
#step two
AW,H,S = cal_matrix(k,init_point)
#step three
err_func, extremepoint ,value,_= cal_err(S,k)
# print(len(extremepoint))
E0.append(max(abs(err_func)))

#step four : calculate max error and loop
while True:
    AW, H, S = cal_matrix(k, extremepoint)
    err_func, extremepoint, value,R = cal_err(S, k)
    E0.append(max(abs(err_func)))

    # step five
    if abs(E0[-1]-E0[-2]) < 0.0001:
        print(E0)
        fig, ax = plt.subplots(1, 2)
        ax[0].set_title("Frequency response")
        ax[0].set_xlabel("Normalized frequence(Hz)")
        ax[0].set_ylabel("Amplitude")
        ax[0].plot(F, R)
        ax[0].plot(F, Hdf)
        break

#step six : impulse response
S = S[:-1]
h = np.zeros(len(S)*2-1)
h[k] = S[0]
for i in range(1,k+1):
    h[k+i] = h[k-i] = S[i]/2
#plt.figure()

ax[1].set_title("Impulse response")
ax[1].set_xlabel("Position")
ax[1].set_ylabel("Amplitude")
ax[1].plot(h) #, label="123")
ax[1].legend()
plt.show()



