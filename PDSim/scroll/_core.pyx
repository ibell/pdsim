

cdef class 
def V_sa(self,theta,full_output=False):
        """
        Wrapper around the Cython code for sa calcs
        
        theta: angle in range [0,2*pi]
        
        Returns
        -------
        
        """
        if full_output==True:
            HTangles = {'1_i':None,'2_i':None,'1_o':None,'2_o':None}
            return ScrollGeo.SA(theta,self.geo)[0:2],HTangles
        else:
            return ScrollGeo.SA(theta,self.geo)[0:2]