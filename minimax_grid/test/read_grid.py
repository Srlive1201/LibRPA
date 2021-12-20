from scipy.optimize import minimize,Bounds,least_squares
import numpy as np
import matplotlib.pyplot as plt
#from scipy.optimize import least_squares
def Fun_minimax(x,args):                        # 定义拟合函数形式
    def eta(x,a):
        sum_x=0.0
        weight=x[:int(len(x)/2)]
        omega=x[int(len(x)/2):]
        i=0
        while i<int(len(x)/2):
            sum_x+=(4*weight[i]*a**2)/(a**2+omega[i]**2)**2
            i=i+1
        return 1/a-sum_x
    #g=reduce(lambda a,b:(1/a-sum((4*x[:int(len(x)/2)]*a**2)/(a**2+x[int(len(x)/2):]**2)**2))+(1/b-sum((4*x[:int(len(x)/2)]*b**2)/(b**2+x[int(len(x)/2):]**2)**2)),args)
    #v=lambda x:1/a-sum((4*x[:int(len(x)/2)]*a**2)/(a**2+x[int(len(x)/2):]**2)**2) 
    list_a=[]
    for a in args:
        list_a.append(abs(eta(x,a)))
    # while i< len(p)/2:
        # sum+=(4*p[i]*x**2)/(x**2+p[i+int(len(p)/2)]**2)**2
        # i=i+1
    return max(list_a)

def Fun(p,x):                        # 定义拟合函数形式
    sum=0
    i=0
    #a1,a2,a3,a4,a5,a6 = p
    while i< len(p)/2:
        sum+=(4*p[i]*x**2)/(x**2+p[i+int(len(p)/2)]**2)**2
        i=i+1
    return sum
    #return a1*x**2+a2*x+a3
    #return (4*a1*x**2)/(x**2+a2**2)**2+(4*a3*x**2)/(x**2+a4**2)**2+(4*a5*x**2)/(x**2+a6**2)**2
def error (p,x,y):                    # 拟合残差
    return Fun(p,x)-y 
def Ori_Fun(x):
    return 1/x
def main():
    #ori_x = np.linspace(3.48,50,10000)  # 创建时间序列
    ori_x=np.linspace(3.1623,1000,10000)
    #x=np.exp(0.2*ori_x)-1
    x=0.1*ori_x**2
    #p_value = [0.1,5,10] # 原始数据的参数
    #noise = np.random.randn(len(x))  # 创建随机噪声
    para_yy=[]
    y = Ori_Fun(x) # 加上噪声的序列
    with open("yy_grid.txt","rt") as f:
        for line in f.readlines():
            para_yy.append(float(line.strip()))
            
    plt.figure
    plt.plot(x,abs(error(para_yy,x,y)),'r', label = 'yy minimax error')
    #plt.plot(x,abs(error(para_mim.x,x,y)),'-b', label ='minimax error')
    plt.legend()
    plt.show()
  
 
if __name__=='__main__':
   main()