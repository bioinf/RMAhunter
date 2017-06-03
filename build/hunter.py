#!/usr/bin/env python
# -*- coding: utf-8 -*-
# ---------------------------------------------------------------------------- #

import sys, os, time
import json, base64
from datetime import datetime
from subprocess import Popen

def fail(error = False) :
    if error : 
        sys.stdout.write('')
        print '\033[0;31mError:'
        print '  ' + error + '\033[0m'

    print '\033[0;37m'
    print 'Usage:'
    print '  \033[1;37m' + sys.argv[0] + '\033[0;37m [input vcf-file] [Optional arguments]'
    print ''
    print 'Optional arguments:'
    print '  \033[1;37m-f\033[0;37m  Path to `input.vcf` file'
    print '  \033[1;37m-c\033[0;37m  Report coding only [0 or 1]. Default: 1'
    print '  \033[1;37m-m\033[0;37m  Allelic frequency cutoff. Default: 0.01'
    print '  \033[1;37m-o\033[0;37m  Output dir name'
    print ''
    print 'Examples:'
    print '  \033[1;37m' + sys.argv[0] + '\033[0;37m input.vcf'
    print '  \033[1;37m' + sys.argv[0] + '\033[0;37m -f input.vcf -c 0 -m 0.05 -o results'
    print '\033[0m'
    sys.exit(1)


# ---------------------------------------------------------------------------- #
# Args parse

try:
    argx = {
        '-f' : sys.argv[1],
        '-b' : '/dev/null',
        '-c' : '1',
        '-m' : '0.01',
        '-o' : 'output_' + datetime.now().strftime("%Y-%m-%d_%H-%M-%S") 
    }
    for i in range(0, len(sys.argv), 2) :
        k, v = [sys.argv[i - 1], sys.argv[i]]
        if k in argx : argx[k] = v
except:
    fail('Аргументы где?')

if not os.path.exists(argx['-f']) :
    fail('vcf file not found (' + argx['-f'] + ')')

if not os.path.exists(argx['-o']):
    os.makedirs(argx['-o'])

if not os.path.exists(argx['-o'] + '/data/'):
    os.makedirs(argx['-o'] + '/data/')

# ---------------------------------------------------------------------------- #
# Get header
# CHROM POS ID REF ALT QUAL FILTER INFO FORMAT sample1 sample2 ...

header = False
params = 9

f = open(argx['-f'])
for line in f:
    if line[0] == '#' :
        if line[1] == '#' or header : continue
        header = line.replace('\n', '').split('\t')[params:]
        break

if not header :
    fail('vcf file format is not valid')


# ---------------------------------------------------------------------------- #
# Separate cols

zygnum = {'1/1' : 1, '0/1' : 2, '0/0' : 3, './1' : 4, './.' : 5}
execrun = []
execsplit = []
dir = os.path.dirname(os.path.realpath(__file__)) + '/'

# header = header[0:2]

for n, sample in enumerate(header) :
    key = argx['-o'] + '/data/_e' + str(n)
    tmp = key + '.xvcf'
    app = dir + 'exec/hunter'
    sdf = dir + 'data/sdf.csv ' + dir + 'data/sdf_plus.csv'
    execrun.append( (' ').join([app, tmp, argx['-b'], argx['-c'], argx['-m'], key, sdf]) )
    for i in ['1','2','3','4'] :
        execsplit.append('cat ' + (key + '.t' + i) + ' | sort | awk {\'print $2\'} > ' + key + '.ts' + i)

    w = open(tmp, 'w+')
    f = open(argx['-f'])
    for line in f:
        if line[0] == '#' : continue
        line = line.replace('\n', '').split('\t')
        chr = line[0].replace('X', '23').replace('Y', '24').replace('M', '25')
        alt = line[4].split(',')
        for alt_x in alt :
            zyg = line[n + params].split(':')[0]
            zyg = zygnum[zyg] if zyg in zygnum else 0
            w.write(('\t').join([chr, line[1], line[3], alt_x, str(zyg)]) + '\n')
    w.close()


# ---------------------------------------------------------------------------- #
# Run HUNTER

def Parallel(exec_list, max_proc, procx = []) :
    def active(procx, cnt = 0):
        for proc in procx : 
            if proc.poll() is None : cnt += 1
        return cnt

    while True :
        time.sleep(0.05)
        if active(procx) + len(exec_list) == 0 : break
        if active(procx) >= max_proc : continue
        if len(exec_list) == 0 : continue
        cmd = exec_list.pop()
        print ' > ' + cmd
        procx.append(Popen(cmd, shell=True))

    return [p.wait() for p in procx]

Parallel(execrun, 4)
Parallel(execsplit, 4)


# ---------------------------------------------------------------------------- #
# Make HTML report

css = '<style>' + open(dir + 'web/client/style.css', "r").read() + '</style>'
ico = 'data:image/x-icon;base64,' + base64.b64encode(open(dir + 'web/client/favicon.ico', "r").read())
jsx = '<script>' + open(dir + 'web/client/app.js', "r").read() + '</script>'

template = open(dir + 'web/index.html', "r").read()
template = template.replace('<link rel="stylesheet" href="client/style.css">', css)
template = template.replace('<script src="client/app.js"></script>', jsx)
template = template.replace('client/favicon.ico', ico)

for n, sample in enumerate(header) :
    data = { 'samples' : header, 'table' : {}, 'current' : sample, 'counts' : [0,0,0,0] }
    for i in ['1','2','3','4'] :
        data['table'][i] = open(argx['-o'] + '/data/_e' + str(n) + '.ts' + i, "r").read()
        data['counts'][int(i) - 1] = len(data['table'][i].split("\n")) - 1

    r = open(argx['-o'] + '/sample_' + str(n) + '.html', 'w+')
    r.write(template.replace('/* data-local */', 'var local_data = ' + json.dumps(data)))
    r.close()

print "Done"
