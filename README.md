pw2py python package
====================

Description
------------------------------------
This code creates an interface between plane-wave based codes (particularly QE) and python. 
It defines several classes and modules which are commonly useful.

Prerequisites:
------------------------------------
* [Python3](https://www.python.org/downloads)
* See doc/requirements.txt for list of required python packages

Installation:
------------------

To install pw2py, navigate to the root of the repository and execute:

```bash
pip install .
```

To uninstall pw2py execute:

```bash
pip uninstall pw2py
```

To update pw2py execute:

```bash
git pull
pip install .
```

(Development) To symlink the package and scripts to local user site execute:

```bash
mkdir -p $(python3 -m site --user-site)
ln -fs $PWD/pw2py $(python3 -m site --user-site)
mkdir -p ~/.local/bin 2> /dev/null
for f in scripts/* ; do ln -fs $PWD/$f ~/.local/bin ; done
```

Examples:
------------------------------------
Try out the examples under the directory [Examples](Examples/)


Author(s)
------------------------------------
Tyler J. Smart

