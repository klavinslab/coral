#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Uses nbconvert to convert all ipynbs in the ipython_examples dir until I
can figure out a way to use a proper sphinx extension. This is really hacky and
should not be used as a secure production server."""
# TODO: catch conversion errors (right now they pass silently)
# Doesn't even use IPython API (TODO!)
import os
from tornado import web, ioloop, httpserver
from ipynb2rst import convert_ipynbs
from build_sphinx_docs import build_docs


DOCSDIR = os.path.abspath(os.path.join(os.path.dirname(__file__), "../docs"))
ROOT = os.path.abspath(os.path.join(DOCSDIR, "_build/html"))
STATIC = os.path.join(ROOT, "_static")


# Web server (Tornado) classes
class MainHandler(web.RequestHandler):
    def get(self, path, *args):
        uri = self.request.uri
        if uri == "/" or uri == "index.html":
            self.render("index.html")
        else:
            self.render(self.request.uri[1:])


class Application(web.Application):
    def __init__(self):
        handlers = [(r"^(.(?!_static))*$", MainHandler),
                    (r"/_static/(.*)", web.StaticFileHandler,
                     {"path": STATIC})]
        settings = {"template_path": ROOT,
                    "static_url_prefix": "_static"}
        web.Application.__init__(self, handlers, **settings)


if __name__ == "__main__":
    # Convert notebooks from ipynb to rst
    convert_ipynbs(DOCSDIR)
    # Build sphinx docs (produces html)
    build_docs(DOCSDIR)
    # Launch server
    applicaton = Application()
    http_server = httpserver.HTTPServer(applicaton)
    http_server.listen(3089)

    ioloop.IOLoop.instance().start()
