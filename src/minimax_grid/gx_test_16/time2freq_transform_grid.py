from scipy.optimize import minimize,Bounds,least_squares
import numpy as np
import matplotlib.pyplot as plt
import math
#from scipy.optimize import least_squares

def Fun_minimax(p,x,tau,freq):                        # 定义拟合函数形式
    def eta(p,xi,freq):
        sum=0.0
        i=0
        while i<len(p):
            sum+=p[i]*np.exp(-1*xi*abs(tau[i]))
            i=i+1
        return Ori_Fun(xi,freq)-sum
    #g=reduce(lambda a,b:(1/a-sum((4*x[:int(len(x)/2)]*a**2)/(a**2+x[int(len(x)/2):]**2)**2))+(1/b-sum((4*x[:int(len(x)/2)]*b**2)/(b**2+x[int(len(x)/2):]**2)**2)),args)
    #v=lambda x:1/a-sum((4*x[:int(len(x)/2)]*a**2)/(a**2+x[int(len(x)/2):]**2)**2) 
    list_a=[]
    for a in x:
        list_a.append(abs(eta(p,a,freq)))
    # while i< len(p)/2:
        # sum+=(4*p[i]*x**2)/(x**2+p[i+int(len(p)/2)]**2)**2
        # i=i+1
    return max(list_a)

def Fun(p,x,tau):                        # 定义拟合函数形式
    sum=0
    i=0
    #a1,a2,a3,a4,a5,a6 = p
    while i< len(p):
        sum+=p[i]*np.exp(-1*x*abs(tau[i]))
        i=i+1
    return sum
    
def error (p,x,y,tau):                    # 拟合残差
    return Fun(p,x,tau)-y 
def Ori_Fun(x,freq):
    return (2*x)/(x**2+freq**2)
def main_trans(R=90):
    omega=[]
    womega=[]
    tau=[]
    wtau=[]
    with open("freq_grid.txt","rt") as f:
        for line in f.readlines():
            a,b=line.split()
            omega.append(float(a))
            womega.append(float(b))
   
    with open("time_grid.txt","rt") as f:
        for line in f.readlines():
            c,d=line.split()
            tau.append(float(c))
            wtau.append(float(d))
            
    num_points_per_magnitude=200
    num_x_node=(int(math.log10(R))+1)*num_points_per_magnitude
    num_x_node=max(num_x_node,len(tau))

    multipicator=(R)**(1.0/(num_x_node-1.0))
    x=[multipicator**(i) for i in range(0,num_x_node)]
    # R_ori=np.sqrt(10*R)        
    # ori_x = np.linspace(3.162,R_ori,10000)
    # #x=np.exp(0.2*ori_x)-1
    # x=0.1*ori_x**2


    
    p0=[]
    i=0
    while i<len(wtau):
        p0.append(np.exp(0.1*i)-1)
        i+=1
    
    with open("INFO.txt","at") as f1:
        print("Time_to_freq_transform      Range_x: ",np.min(x),np.max(x),file=f1)
    with open("time2freq_transform_grid.txt","wt") as f0:
        for i_freq in omega:
            print("Freq: ",i_freq,file=f0)
            y=Ori_Fun(x,i_freq)
            #para_i=least_squares(error, p0, bounds=[0,np.inf],args=(x,y,tau))
            #para_min =minimize(Fun_minimax, para_i.x ,args=(x,tau,i_freq),bounds=Bounds(0,np.inf)) # 进行拟合
            para_i=least_squares(error, p0, args=(x,y,tau))
            #para_min =minimize(Fun_minimax, para_i.x ,args=(x,tau,i_freq)) # 进行拟合
            #print("err: ",np.max(abs(para_i.fun)),"       ",np.max(abs(para_min.fun))) 
            for gamma in para_min.x:
                print(gamma,file=f0)
            with open("INFO.txt","at") as f1:
                print("Freq_point: ",i_freq,"         transform_err: ",np.max(abs(para_min.fun)),file=f1)
    with open("INFO.txt","at") as f1:
        print("Finish!",file=f1)
   
if __name__=='__main__':
   main_trans()
