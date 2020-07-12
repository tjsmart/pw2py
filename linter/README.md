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
3. Copy or symlink the vimscipt file `qelint.vim` to the ale fixers folder you found above.
```
cp qelint.vim ~/.vim/plugged/ale/autoload/ale/fixers            # copy
ln -fs $PWD/qelint.vim ~/.vim/plugged/ale/autoload/ale/fixers   # symlink
```
4. Locate the ale `registry.vim` file:
```
~/.vim/plugged/ale/autoload/ale/fix/registry.vim
```
5. Open the file and 
