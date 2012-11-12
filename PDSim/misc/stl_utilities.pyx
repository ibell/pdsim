from cython.operator cimport dereference as deref, preincrement as inc

cpdef bool is_in_vector(string a, vector[string] v):
    cdef vector[string].iterator it = v.begin()
    while it != v.end():
        if deref(it) == a:
            return True
        #Increment iterator
        inc(it)
    return False
    
cpdef double get_map_sd(map[string, double] m, string s):
    cdef map[string,double].iterator it = m.find(s)
    if it != m.end():
        return deref(it).second
    else:
        return 1e99