import matplotlib.pyplot as plt
import sys
import numpy as np

t=[]
w=[]

def reverse_x_fun(x_lsit):
    return [1.0/(i) for i in x_lsit]

def reverse_2x_fun(x_lsit):
    return [1.0/(2*i) for i in x_lsit]

def exp_fit_fun(t_list,w_list,x_list):
    def exp_fit(t_list,w_list,x_value):
        sum_x=0.0
        i=0
        while i<len(t_list):
            sum_x+=w_list[i]*np.exp(-2*x_value*abs(t_list[i]))
            i+=1
        return sum_x
    
    return [exp_fit(t_list,w_list,x) for x in x_list]

def omega_fit_fun(f_list,w_list,x_list):
    def omega_fit(f_list,w_list,x_value):
        sum_x=0.0
        i=0
        while i<len(f_list):
            sum_x+=(4*w_list[i]*x_value**2)/(x_value**2+f_list[i]**2)**2
            i+=1
        return sum_x
    
    return [omega_fit(f_list,w_list,x) for x in x_list]





if __name__=="__main__":
    file =sys.argv[1]
    grid_type=file.split('/')[-1][0]
    print(grid_type)
    with open(file) as fs:
        lines=fs.readlines()
        for line in lines:
            a,b=line.split()
            t.append(eval(a))
            w.append(eval(b))
    print(t)
    print(w)
    x_list=np.linspace(1,100000,300000)
    if grid_type=='t':
        fit_curve=exp_fit_fun(t,w,x_list)
        goal=reverse_2x_fun(x_list)
    elif grid_type=='f':
        fit_curve=omega_fit_fun(t,w,x_list)
        goal=reverse_x_fun(x_list)
    
    err=[abs(fit_curve[i]-goal[i]) for i in range(0,len(x_list))]
    print('MAX Err: ',max(err))
    #goal=reverse_x_fun(x_list)
    #goal2=reverse_2x_fun(x_list)
    
    plt.plot(x_list,err,label='err')
    #plt.plot(x_list,fit_curve,label='fit_curve')
    #plt.plot(x_list,goal,label='1/x')
    #plt.plot(x_list,goal2,label='1/2x')
    plt.legend()
    plt.show()
