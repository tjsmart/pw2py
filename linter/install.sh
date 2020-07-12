#!/usr/bin/env bash

help_and_quit(){
    echo "Description: installs qelint to /usr/local/bin, link to ale, and update vimrc"
    echo "Usage: ./install.sh [-h] [-l] /path/to/ale/"
    exit 0
}

arg_parse(){
    if [ -z "$1" ]; then
        help_and_quit
    fi
    uselink=false
    while [ -n "$1" ]; do
        if [ "$1" == "-h" ] || [ "$1" == "--help" ]; then
            help_and_quit
        elif [ "$1" == "-l" ] || [ "$1" == "--link" ]; then
            uselink=true
            shift
        else
            path_to_ale="$1"
            break
        fi
    done
    if [ -z $path_to_ale ]; then
        echo "Please provide the path to ale" >&2
        exit 1
    fi
}

check_ale(){
    # check provided folder exists
    if [ ! -d "$path_to_ale" ]; then
        echo "Folder provided does not exist: $path_to_ale" >&2
        exit 1
    fi
    # check for registry file
    registry_file="$path_to_ale/autoload/ale/fix/registry.vim"
    if [ ! -f "$registry_file" ]; then
        echo "Registry file not found: '$registry_file'" >&2
        exit 1
    fi
    # check for fixers folder
    fixers_folder="$path_to_ale/autoload/ale/fixers"
    if [ ! -d "$fixers_folder" ]; then
        echo "Fixers folder not found: '$fixers_folder'" >&2
        exit 1
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
        echo "Desired pattern was not found in '$1', are you sure this is the right file?" >&2
        exit 1
    fi
    # check if inst is already in file
    if grep -q "qelint" "$1"; then
        echo "qelint already insterted into file '$1' ... skipping" >&2
    else
        # otherwise insert 'inst' after 'ptrn'
        perl -i -pe "s/$ptrn/$ptrn$inst/;" $1
    fi
}

cp_or_ln_qelint(){
    # check for file to copy
    local vimfile="qelint.vim"
    if [ ! -f $vimfile ]; then
        echo "File 'qelint.vim' not found, exiting" >&2
        exit 1
    fi
    # check for file destination
    local file_dest="$fixers_folder/$vimfile"
    if $uselink ; then
        ln -fs $PWD/$vimfile $file_dest
    else
        cp $vimfile $file_dest
    fi
}

configure_ale(){
    check_ale
    check_perl
    insert_after_pattern $registry_file
    cp_or_ln_qelint $fixers_folder
}

install_qelint(){
    if $uselink ; then
        ln -fs $PWD/qelint /usr/local/bin
    else
        cp qelint /usr/local/bin
    fi
}


ptrn="let s:default_registry = {"
inst="\n\
\\\   'qelint': {\n\
\\\       'function': 'ale#fixers#qelint#Fix',\n\
\\\       'suggested_filetypes': ['espresso'],\n\
\\\       'description': 'Format quantum espresso input files',\n\
\\\   },"

arg_parse $@
echo "Installing qelint to /usr/local/bin ..."
install_qelint
echo "Configuring users ale ..."
configure_ale
echo "Script ran successfully! :)"
echo "--------------------------"

echo "Almost done! Don't forget to add these lines to your ~/.vimrc!"
echo "    let g:ale_fixers = {'espresso': ['qelint']}"
echo "    let g:ale_fix_on_save = 1"
