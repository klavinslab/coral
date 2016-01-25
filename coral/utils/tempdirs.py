# -*- coding: utf-8
'''Temporary directory helpers for scripts that call command line
applications. '''
import os
import shutil
import tempfile


def tempdir(fun):
    '''For use as a decorator of instance methods - creates a temporary dir
    named self._tempdir and then deletes it after the method runs.

    :param fun: function to decorate
    :type fun: instance method

    '''
    def wrapper(*args, **kwargs):
        self = args[0]
        if os.path.isdir(self._tempdir):
            shutil.rmtree(self._tempdir)
        self._tempdir = tempfile.mkdtemp()
        # If the method raises an exception, delete the temporary dir
        try:
            retval = fun(*args, **kwargs)
        finally:
            shutil.rmtree(self._tempdir)
        if os.path.isdir(self._tempdir):
            shutil.rmtree(self._tempdir)
        return retval
    return wrapper
