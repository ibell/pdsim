#cython: embedsignature=True

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
        
def Broyden(f, x0, dx=1e-5, args=(), ytol=1e-5, Nfd = 1, w=1.0, wJ=1.0, itermax=10, JustOneStep=False):
    """
    Broyden's method
    
    If f returns ``None``, then the computation is stopped, and a list the size of x0 is returned with all ``None`` values
    
    Parameters
    ----------
    f : function
        Must have a signature of the form::
         
            (x0,*args) --> array
            
        that returns an array-like object the same shape as ``x0``
        
    x0 : array-like
        Initial values for the solver
    args : tuple
        Positional arguments to pass to the objective function 
    ytol : float
        The allowed root-sum-square-error
    itermax : int
        maximum number of iterations allowed
    Nfd : int
        The number of finite difference steps to be used at the beginning
    w : float
        relaxation factor
    wJ : float
        relaxation factor for updating of Jacobian matrix
    JustOneStep : boolean
        If ``True``, stop after one step
    """
    x0=np.array(x0,dtype=np.float)
    x1=x0.copy()
    error=999
    A1=np.zeros((len(x0),len(x0)))
    A0=np.zeros((len(x0),len(x0)))
    
    if isinstance(dx,float) or isinstance(dx,int):
        dx=dx*np.ones_like(x0)
    error_vec=[]
    
    def not_improving():
        # If the error increased in the last step, re-evaluate the Jacobian
        if len(error_vec)>2 and error_vec[-1]>error_vec[-2]:
            return True
        else:
            return False
    
    iter = 1
    while abs(error)>ytol:
        if iter <= Nfd or not_improving():
            # If past the first numerical derivative step, use the value of x
            # calculated below 
            if iter > 1:
                x0 = x1.copy()
            f0=f(x0,*args)
            F0=array(f0)
            if f0 is None:
                return [None,None]
            
            #Build the Jacobian matrix by columns
            for i in range(len(x0)):
                epsilon=np.zeros_like(x0)
                epsilon[i]=dx[i]
                fplus = f(x0+epsilon,*args)
                if fplus is None:
                    return [None,None]
                A0[:,i]=(array(fplus)-F0)/epsilon[i]
            #Get the difference vector
            x1=x0-w*np.dot(inv(A0),F0)
            #Just do one step and stop
            if JustOneStep==True:
                return x1
            iter+=1
            # Flush the error vector to ensure that you don't continually 
            # rebuild the Jacobian
            error_vec=[]
            
        elif iter>1 and iter<itermax:
            #Jacobian updating parameters
            S=x1-x0
            d=np.dot(S.T,S)
            f1=f(x1,*args)
            F1=array(f1)
            if f1 is None:
                return [None,None]
            Y=F1-F0
            A1=A0+w/d*dot((Y-dot(A0,S)),S.T)
            x2=x1-w*np.dot(inv(A1),F1)
            #Update values
            x0=x1
            x1=x2
            F0=F1
            A0=A1
            error=np.sqrt(np.sum(np.power(F1,2)))
            error_vec.append(error)
            print 'error_vec',error_vec
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