#!/usr/bin/env bash

help_and_quit(){
    echo "Description: installs qelint to /usr/local/bin, link to ale, and update vimrc"
    echo "Usage: ./install.sh [-h] [-l] /path/to/ale/"
    exit 0
}

arg_parse(){
    if [ -z "$1" ] || [ "$1" == "-h" ] || [ "$1" == "--help" ]; then
        help_and_quit
    fi
    # TODO parse other options (i.e. linking [l])
}

check_ale(){
    # check provided folder exists
    if [ ! -d "$1" ]; then
        echo "Folder provided does not exist: $1" >&2
        exit 1
    fi
    # check for registry file
    registry_file="$1/autoload/ale/fix/registry.vim"
    if [ ! -f "$registry_file" ]; then
        echo "Registry file not found: '$registry_file'"
    fi
    # check for fixers folder
    fixers_folder="$1/autoload/ale/fixers"
    if [ ! -d "$fixers_folder" ]; then
        echo "Fixers folder not found: '$fixers_folder'"
    fi
}

check_perl(){
    # check for perl
    if ! perl --version > /dev/null ; then
        echo "Command not found: perl" >&2
        echo "perl is required to run this install script" >&2
        exit 1
    fi
}

insert_after_pattern(){
    # check that ptrn is in file
    if ! grep -q "$ptrn" "$1"; then
        echo "Desired pattern was not found in '$1', are you sure this is the right file?"
        exit 1
    fi
    # check if inst is already in file
    if grep -q "qelint" "$1"; then
        echo "qelint already insterted into file '$1' ... skipping"
    else
        # otherwise insert 'inst' after 'ptrn'
        perl -i -pe "s/$ptrn/$ptrn$inst/;" $1
    fi
}

cp_or_ln_qelint(){
    # check for file
    echo "TODO"
    # TODO cp qelint.vim $fixers_folder
}


ptrn="let s:default_registry = {"
inst="\n\
\\\   'qelint': {\n\
\\\       'function': 'ale#fixers#qelint#Fix',\n\
\\\       'suggested_filetypes': ['espresso'],\n\
\\\       'description': 'Format quantum espresso input files',\n\
\\\   },"

arg_parse $@
check_ale $@
check_perl
insert_after_pattern $registry_file
cp_or_ln_qelint $fixers_folder
echo "Almost done! Don't forget to add these lines to your ~/.vimrc!"
echo "let g:ale_fixers = {'espresso': ['qelint']}"
echo "let g:ale_fix_on_save = 1"
