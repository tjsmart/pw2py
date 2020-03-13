pw2py python module
====================


Install this module
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

