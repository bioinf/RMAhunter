#!/usr/bin/env python
# -*- coding: utf-8 -*-
# ---------------------------------------------------------------------------- #

import os, base64

dir = os.path.dirname(os.path.realpath(__file__)) + '/../web/client/'

css = '<style>' + open(dir + 'style.css', "r").read() + '</style>'
icn = 'data:image/x-icon;base64,' + base64.b64encode(open(dir + 'favicon.ico', "r").read())
jst = '<script>' + open(dir + 'app.js', "r").read() + '</script>'

template = open(dir + 'template.html', "r").read()
template = template.replace('<link rel="stylesheet" href="client/style.css">', css)
template = template.replace('<script src="client/app.js"></script>', jst)
template = template.replace('client/favicon.ico', icn)
template = ' '.join(template.replace('\n', '').split())

r = open(dir + '../index.html', 'w+')
r.write(template)
r.close()

print "Ok"
