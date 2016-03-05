#ifndef __PYX_HAVE__PDSim__flow__flow_models
#define __PYX_HAVE__PDSim__flow__flow_models


/* "PDSim\flow\flow_models.pyx":13
 * """
 * 
 * cdef public enum:             # <<<<<<<<<<<<<<
 *     TYPE_RADIAL
 *     TYPE_FLANK
 */
enum  {
  TYPE_RADIAL,
  TYPE_FLANK,
  TYPE_DISABLED
};

/* "PDSim\flow\flow_models.pyx":18
 *     TYPE_DISABLED
 * 
 * cdef public enum:             # <<<<<<<<<<<<<<
 *     OUTPUT_VELOCITY
 *     OUTPUT_MA
 */
enum  {
  OUTPUT_VELOCITY,
  OUTPUT_MA
};

#ifndef __PYX_HAVE_API__PDSim__flow__flow_models

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

#endif /* !__PYX_HAVE_API__PDSim__flow__flow_models */

#if PY_MAJOR_VERSION < 3
PyMODINIT_FUNC initflow_models(void);
#else
PyMODINIT_FUNC PyInit_flow_models(void);
#endif

#endif /* !__PYX_HAVE__PDSim__flow__flow_models */
