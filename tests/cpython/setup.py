from distutils.core import setup, Extension

api_dir = "/home/gpenazzi/codes/dev_libnegf/src/api/"
libraries = ["python2.7", "negf_x86_64", "lapack", "f77blas", "cblas", "atlas",
"gomp", "pthread","stdc++"]

# define the extension module
pynegf = Extension('cpynegf', sources=['cpynegf.c'],
            include_dirs = [api_dir],
            libraries = libraries,
            extra_link_args = ["-L"+api_dir, "-Wl,--no-undefined"])

setup(ext_modules=[pynegf])
