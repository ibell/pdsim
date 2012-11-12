from libcpp.string cimport string
from libcpp.vector cimport vector
from libcpp.map cimport map
from libcpp cimport bool

cpdef bool is_in_vector(string a, vector[string] v)
cpdef double get_map_sd(map[string, double] m, string s)