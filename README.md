pw2py python module
====================

Description
------------------------------------
This code creates an interface between plane-wave based codes (particularly QE) and python. 
It define several classes and methods which are commonly useful.

Prerequisites:
------------------------------------
* [Python3](https://www.python.org/downloads)
* Python Packages: [f90nml](https://pypi.org/project/f90nml/), [numpy](https://pypi.org/project/numpy/), [pandas](https://pypi.org/project/numpy/), [mendeleev](https://pypi.org/project/mendeleev/)

Install this module:
------------------

```bash
# create directory for local python modules
mkdir -p $(python -m site --user-site) 2> /dev/null

# EITHER create a link (suggested for development)
for f in *.py ; do
    ln -fs $(realpath $f) $(python -m site --user-site)/$f
done

# OR copy files (suggested for permanent copy)
cp *.py $(python -m site --user-site)
```

Examples:
------------------------------------
Try out the examples under the directory [Examples](Examples/)


Author(s)
------------------------------------
Tyler J. Smart


