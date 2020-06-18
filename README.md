pw2py python module
====================

Description
------------------------------------
This code creates an interface between plane-wave based codes (particularly QE) and python. 
It defines several classes and modules which are commonly useful.

Prerequisites:
------------------------------------
* [Python3](https://www.python.org/downloads)
* Python Packages: [f90nml](https://pypi.org/project/f90nml/), [numpy](https://pypi.org/project/numpy/), [pandas](https://pypi.org/project/numpy/), [mendeleev](https://pypi.org/project/mendeleev/), [lxml](https://pypi.org/project/lxml/), [scipy](https://pypi.org/project/scipy/)

Install this module:
------------------

```bash
# create directory for local python modules
mkdir -p $(python -m site --user-site) 2> /dev/null

# EITHER create a link to the package (suggested for development)
ln -fs $(realpath ./pw2py) $(python -m site --user-site)/$f

# OR copy the package (suggested for permanent copy)
cp -r ./pw2py $(python -m site --user-site)
```

Examples:
------------------------------------
Try out the examples under the directory [Examples](Examples/)


Author(s)
------------------------------------
Tyler J. Smart

