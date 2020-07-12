qelint
======

qelint uses pw2py formatting to lint quantum espresso input files


Basic usage
====

To lint a qe input file:
```bash
qelint scf.in
```

Vim usage (with ale)
===

0. Install [vim-qe](https://github.com/tjsmart/vim-qe), a vim plugin for quantum espresso input files.
1. Install [ale](https://github.com/dense-analysis/ale) from their git repo.
2. Find the location of the ale repository. For example:
```
~/.vim/plugged/ale/autoload/ale/fixers
```
3. Use the `install.sh` to install add the qelint script to ale fixers.
```
./install.sh -l ~/.vim/plugged/ale/autoload/ale/fixers  # using symlinks
# OR
./install.sh ~/.vim/plugged/ale/autoload/ale/fixers     # copies source files
```
