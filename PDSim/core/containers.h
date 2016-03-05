#ifndef __PYX_HAVE__PDSim__core__containers
#define __PYX_HAVE__PDSim__core__containers


/* "PDSim\core\containers.pyx":5
 * cimport cython
 * 
 * cdef public enum STATE_VARS:             # <<<<<<<<<<<<<<
 *     STATE_VARS_TD
 *     STATE_VARS_TM
 */
enum STATE_VARS {
  STATE_VARS_TD,
  STATE_VARS_TM
};

#ifndef __PYX_HAVE_API__PDSim__core__containers

#ifndef __PYX_EXTERN_C
  #ifdef __cplusplus
    #define __PYX_EXTERN_C extern "C"
  #else
    #define __PYX_EXTERN_C extern
  #endif
#endif

#ifndef DL_IMPORT
  #define DL_IMPORT(_T) _T
#endif

#endif /* !__PYX_HAVE_API__PDSim__core__containers */

#if PY_MAJOR_VERSION < 3
PyMODINIT_FUNC initcontainers(void);
#else
PyMODINIT_FUNC PyInit_containers(void);
#endif

#endif /* !__PYX_HAVE__PDSim__core__containers */
