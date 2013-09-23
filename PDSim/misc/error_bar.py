import numpy as np
from math import log10

def error_ascii_bar(val, tol, N = 50, min_p = -5, max_p = 2):

    assert log10(tol) > min_p
    
    base_string = list('||'+'.'*N+'||')
    
    if np.isnan(val):
        return ' XXX VAL is NAN XXX '
    #Add the bar for the tolerance
    i = int(N*(log10(tol)-min_p)/(max_p-min_p)+1)
    base_string[i] = '|'
    
    #Add the bar for the value
    i = int(N*(log10(val)-min_p)/(max_p-min_p)+1)
    base_string[i] = '@'
    
    base_string = ''.join(base_string)
    
    return base_string

if __name__=='__main__':
    print error_ascii_bar(1,1e-3)
    print error_ascii_bar(0.5,1e-3)
    print error_ascii_bar(0.1,1e-3)
    print error_ascii_bar(0.0011,1e-3)