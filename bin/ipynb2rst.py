#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Uses nbconvert to recursively convert all ipynbs in a directory to .rst."""
# TODO: catch conversion errors (right now they pass silently)
# Doesn't even use IPython API (TODO!)
import os
import subprocess
import sys


# Build docs
def ipynb_to_rst(directory, filename):
    """Converts a given file in a directory to an rst in the same directory."""
    os.chdir(directory)
    subprocess.Popen(["ipython", "nbconvert", "--to", "rst",
                      filename],
                     stdout=subprocess.PIPE,
                     stderr=subprocess.PIPE,
                     cwd=directory)


def convert_ipynbs(directory):
    """Recursively converts all ipynb files in a directory into rst files in
    the same directory."""
    # The ipython_examples dir has to be in the same dir as this script
    for root, subfolders, files in os.walk(os.path.abspath(directory)):
        for f in files:
            if ".ipynb_checkpoints" not in root:
                if f.endswith("ipynb"):
                    ipynb_to_rst(root, f)


if __name__ == "__main__":
    # Convert notebooks from ipynb to rst
    if len(sys.argv) != 2:
        print "\nipy2nb:\n=======\nUsage: ipynb2rst <directory>\n"
    else:
        directory = sys.argv[1]
        convert_ipynbs(directory)
