from __future__ import division
import numpy as np
from numpy.linalg import inv
from numpy import array,dot

def MultiDimNewtRaph(f,x0,dx=1e-6,args=(),ytol=1e-5,w=1.0,JustOneStep=False):
    """
    A Newton-Raphson solver where the Jacobian is always re-evaluated rather than
    re-using the information as in the fsolve method of scipy.optimize
    """
    #Convert to numpy array, force type conversion to float
    x=np.array(x0,dtype=np.float)
    error=999
    J=np.zeros((len(x),len(x)))
    
    #If a float is passed in for dx, convert to a numpy-like list the same shape
    #as x
    if isinstance(dx,int) or isinstance(dx,float):
        dx=dx*np.ones_like(x0)
    
    r0=array(f(x,*args))
    while abs(error)>ytol:
        #Build the Jacobian matrix by columns
        for i in range(len(x)):
            epsilon=np.zeros_like(x)
            epsilon[i]=dx[i]
            J[:,i]=(array(f(x+epsilon,*args))-r0)/epsilon[i]
        v=np.dot(-inv(J),r0)
        x=x+w*v
        #Calculate the residual vector at the new step
        r0=f(x,*args)
        error = np.max(np.abs(r0))
        #Just do one step and stop
        if JustOneStep==True:
            return x
    return x
        
def Broyden(f,x0,dx=1e-5,args=(),ytol=1e-5,w=1.0,itermax=10,JustOneStep=False):
    """
    Broyden's method
    
    If f returns ``None``, then the computation is stopped, and a list the size of x0 is returned with all ``None`` values
    """
    x0=np.array(x0,dtype=np.float)
    x1=x0.copy()
    error=999
    A1=np.zeros((len(x0),len(x0)))
    A0=np.zeros((len(x0),len(x0)))
    
    #If a float is passed in for dx, convert to a numpy-like list the same shape
    #as x
    if isinstance(dx,float) or isinstance(dx,int):
        dx=dx*np.ones_like(x0)
    f0=f(x0,*args)
    F0=array(f0)
    if f0 is None:
        return [None,None]
    iter=1
    x1=x0.copy()
    F1=F0.copy()
    while abs(error)>ytol:
        if iter==1:
            #Build the Jacobian matrix by columns
            for i in range(len(x0)):
                epsilon=np.zeros_like(x0)
                epsilon[i]=dx[i]
                fplus = f(x0+epsilon,*args)
                if fplus is None:
                    return [None,None]
                A0[:,i]=(array(fplus)-F0)/epsilon[i]
            #Get the difference vector
            x1=x0-np.dot(inv(A0),F0)
            #Just do one step and stop
            if JustOneStep==True:
                return x1
            iter+=1
        elif iter>1 and iter<itermax:
            #Jacobian updating parameters
            S=x1-x0
            d=np.dot(S.T,S)
            f1=f(x1,*args)
            F1=array(f1)
            if f1 is None:
                return [None,None]
            Y=F1-F0
            A1=A0+1.0/d*dot((Y-dot(A0,S)),S.T)
            x2=x1-np.dot(inv(A1),F1)
            #Update values
            x0=x1
            x1=x2
            F0=F1
            A0=A1
            error=np.sqrt(np.sum(np.power(F1,2)))
            iter+=1
        else:
            raise ValueError('Reached maximum number of iterations without getting below ytol RMS error')
            return np.nan*np.ones_like(x0)
    return x1
                
if __name__=='__main__':
    pass
##     def OBJECTIVE(x):
##         return [x[0]**2-2*x[1]-3.7, x[0]+x[1]**2-1]
##     from scipy.optimize import fsolve
##     _x=fsolve(OBJECTIVE,[-1.0,1.0]); print _x
##     print OBJECTIVE(_x)
##     _x=MultiDimNewtRaph(OBJECTIVE,[-1,1],ytol=1e-10); print _x
##     print OBJECTIVE(_x)
#    def f(x):
#        print '.',
#        return [1-4*x[0]+2*x[0]**2-2*x[1]**3,-4+x[0]**4+4*x[1]+4*x[1]**4]
#        
#    _x=Broyden(f,[0.1,0.7],ytol=1e-8); print _x
#    print f(_x)
#    from scipy.optimize import fsolve
#    _x=fsolve(f,[0.1,1.0]); print _x
#    print f(_x)