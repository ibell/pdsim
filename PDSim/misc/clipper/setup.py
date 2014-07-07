
if __name__=='__main__':
    from distutils.core import setup
    from distutils.extension import Extension
    from Cython.Distutils import build_ext

    import sys
    sys.argv += ['build_ext','--inplace']

    ext =	Extension("pyclipper", 
                    sources=["pyclipper.pyx", "clipper.cpp"],
                    language="c++",              # this causes Pyrex/Cython to create C++ source
                    include_dirs=["./../include"],
                    )

    setup(
            ext_modules=[ext],
            cmdclass = {'build_ext': build_ext},
    )


