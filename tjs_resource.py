#!/usr/bin/env python

'''
basic common functions used in many tjs-toolkit programs
'''

import sys
import os
import signal

def handler(signal_received, frame):
    sys.stderr.write('\nSIGINT or CTRL-C detected. Exiting gracefully\n')
    sys.exit(1)


def die(message):
    sys.stderr.write("Error: {}\n".format(message))
    sys.exit(1)


def warn(message):
    sys.stderr.write('Warning: {}\n'.format(message))
    return None


def checkFile(filename):
    if not os.path.isfile(filename):
        die('File \'{}\' does not exist'.format(filename))


def checkFolder(foldername):
    if not os.path.isdir(foldername):
        die('Folder \'{}\' does not exist'.format(foldername))


def grep(var, filename):
    # Only intended for single instance! Returns error if var not found
    with open(filename) as f:
        for line in f.readlines():
            if var in line:
                return line
    # line not found:
    die("in grep: '{}' not found in '{}'\n".format(var, filename))

