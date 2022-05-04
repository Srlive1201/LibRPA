from scipy.optimize import minimize,Bounds,least_squares
import numpy as np
import matplotlib.pyplot as plt
#from scipy.optimize import least_squares
def Fun_minimax(x,args):                        # 定义拟合函数形式
    def eta(x,a):
        sum_x=0.0
        weight=x[:int(len(x)/2)]
        tau=x[int(len(x)/2):]
        i=0
        while i<int(len(x)/2):
            sum_x+=weight[i]*np.exp(-2*a*abs(tau[i]))
            i=i+1
        return 1/(2*a)-sum_x
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
        sum+=p[i]*np.exp(-2*x*abs(p[i+int(len(p)/2)]))
        i=i+1
    return sum
    #return a1*x**2+a2*x+a3
    #return (4*a1*x**2)/(x**2+a2**2)**2+(4*a3*x**2)/(x**2+a4**2)**2+(4*a5*x**2)/(x**2+a6**2)**2
def error (p,x,y):                    # 拟合残差
    return Fun(p,x)-y 
def Ori_Fun(x):
    return 1/(2*x)
def main_time(Npoint=8,R=90):
    #ori_x = np.linspace(2,50,20000)  # 创建时间序列
    R_ori=np.sqrt(10*R)
    ori_x = np.linspace(3.162,R_ori,10000)
    #x=np.exp(0.2*ori_x)-1
    x=0.1*ori_x**2
    #p_value = [0.1,5,10] # 原始数据的参数
    #noise = np.random.randn(len(x))  # 创建随机噪声
    y = Ori_Fun(x) # 加上噪声的序列
    p0=[]
    #Npoint=8
    i=0
    j=0
    while i < Npoint:
        p0.append(np.exp(0.1*i)-1)
        i=i+1
    while j < Npoint:
        p0.append(np.exp(0.1*j)-1)
        j=j+1
    #p0 = [0.1,0.1,0.1,0.1,0.1,0.1] # 拟合的初始参数设置
    para =least_squares(error, p0, bounds=[0,np.inf],args=(x,y)) # 进行拟合
    #y_fitted = Fun (para.x,x) # 画出拟合后的曲线
 
    k=0
    grid=[]
    while k<len(para.x)/2:
        grid.append([para.x[k+int(len(para.x)/2)],para.x[k]])
        k=k+1
    grid.sort(key=lambda x:x[0])
 #   print("grid:",grid)
    print("Leastsq-grid:")
    for a in grid:
        print(a[0],",     ",a[1],",")
        
        
    para_mim =minimize(Fun_minimax, para.x ,args=x,bounds=Bounds(0,np.inf)) # 进行拟合
   
    q=0
    grid_mim=[]
    while q<len(para.x)/2:
        grid_mim.append([para_mim.x[q+int(len(para_mim.x)/2)],para_mim.x[q]])
        q=q+1
    grid_mim.sort(key=lambda x:x[0])
    print("Minimax-grid:")
    for b in grid_mim:
        print(b[0],",     ",b[1],",")
    
    with open("time_grid.txt", "wt") as f:  
        for b in grid_mim:
            print(b[0],"         ",b[1],file=f)
            
    with open("INFO.txt","at") as f1:
        print("time_grid_num: ",Npoint,"      Range_x: ",np.min(x),np.max(x),"          err: ",np.max(abs(para_mim.fun)),file=f1)
'''    
    y_fitted_mini=Fun(para_mim.x,x)
    plt.figure
    plt.plot(x,abs(error(para.x,x,y)),'r', label = 'leastsq error')
    plt.plot(x,abs(error(para_mim.x,x,y)),'-b', label ='minimax error')
    plt.legend()
    plt.show()
'''  
        
 
if __name__=='__main__':
   main_time()
