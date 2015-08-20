#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Builds sphinx docs for pymbt."""
import os
import subprocess


DOCSDIR = os.path.abspath(os.path.join(os.path.dirname(__file__), "../docs"))


def build_docs(directory):
    """Builds sphinx docs from a given directory."""
    os.chdir(directory)
    process = subprocess.Popen(["make", "html"], cwd=directory)
    process.communicate()


if __name__ == "__main__":
    # Build sphinx docs (produces html)
    build_docs(DOCSDIR)
