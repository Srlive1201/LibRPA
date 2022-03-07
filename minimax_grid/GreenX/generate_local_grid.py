import os
from scipy.optimize import minimize,Bounds,least_squares
import numpy as np
import matplotlib.pyplot as plt
import math
import sys

def read_grid(N,grid_type):
    if grid_type == 'f':
        fm=str(N)+"_freq_points.dat"
    elif grid_type == 't':
        fm=str(N)+"_time_points.dat"
    with open(path+fm,'r') as ff:
        while True:
            lines=ff.readline()
            if lines[:6]=='Erange':
                ss,lb,ub=lines.split()
                if R>=eval(lb) and R<=eval(ub):
                    lines=next(ff)
                    i=2*N
                    with open('local_'+fm,'w') as fw:
                        while i>0:
                            lines=next(ff)
                            print(lines.strip(),file=fw)
                            if grid_type == 'f':
                                omega.append(eval(lines.strip()))
                            elif grid_type == 't':
                                tau.append(eval(lines.strip()))
                            i-=1
                    break

def cosine_time2freq_grid(R):
    def tau_fun(p,x,tau):
        sum=0
        i=0
        while i< len(p):
            sum+=p[i]*np.exp(-1*x*abs(tau[i]))
            i=i+1
        return sum

    def omega_fun(x,freq):
        return (2*x)/(x**2+freq**2)
    def error (p,x,y,tau):                    # 拟合残差
        return tau_fun(p,x,tau)-y
    num_points_per_magnitude=200
    num_x_node=(int(math.log10(R))+1)*num_points_per_magnitude
    num_x_node=max(num_x_node,len(tau)/2)

    multipicator=(R)**(1.0/(num_x_node-1.0))
    xvalue=[multipicator**(i) for i in range(0,num_x_node)]
    x=np.array(xvalue)

    p0=[]
    i=0
    while i<len(tau)/2:
        #p0.append(np.exp(0.1*i)-1)
        p0.append(1)
        i+=1
    print(p0)
    with open("INFO.txt","at") as f1:
        print("Time_to_freq_transform   N:  ",N,"   Range_x: ",np.min(x),np.max(x),file=f1)
    with open(str(N)+"_time2freq_grid.txt","wt") as f0:
        omega_g=omega[:N]
        for i_freq in omega_g:
            print("Freq: ",i_freq,file=f0)
            y=omega_fun(x,i_freq)
            para_i=least_squares(error, p0,args=(x,y,tau))
            for gamma in para_i.x:
                print(gamma,file=f0)
            with open("INFO.txt","at") as f1:
                print("Freq_point: ",i_freq,"         transform_err: ",np.max(abs(para_i.fun)),file=f1)
    with open("INFO.txt","at") as f1:
        print("Finish!",file=f1)

if __name__=='__main__':
    N=int(sys.argv[1])
    R=float(sys.argv[2])
    print("N: ",N)
    print("R: ",R)
    omega=[]
    tau=[]
    path=os.path.dirname(__file__) + "/"
    freq_grid=[]
    read_grid(N,'f')
    read_grid(N,'t')
    cosine_time2freq_grid(R)
