#!/usr/bin/env python
# -*- coding: utf-8 -*-
# ---------------------------------------------------------------------------- #

import sys, os, time
import multiprocessing
import json, base64

from array import array
from datetime import datetime

from subprocess import Popen

def fail(error = False) :
    if error : 
        print '\033[0;31m'
        print 'Error:'
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
    print '  \033[1;37m-t\033[0;37m  Number of threads. Default: ' + str(multiprocessing.cpu_count())
    print ''
    print 'Examples:'
    print '  \033[1;37m' + sys.argv[0] + '\033[0;37m input.vcf'
    print '  \033[1;37m' + sys.argv[0] + '\033[0;37m -f input.vcf -c 0 -m 0.05 -o results'
    print '\033[0m'
    sys.exit(1)


class Write :
    def __init__(self, prefix_dir = "") :
        self.prefix_dir = prefix_dir
        self.files = {}

    def write(self, name) :
        f = open(self.prefix_dir + name, 'a+')
        f.write(self.files[name]['txt'])
        f.close()
        self.files[name]['cnt'] = 0
        self.files[name]['txt'] = ""

    def add(self, name, data) :
        if name not in self.files : 
            self.files[name] = { 'cnt' : 0, 'txt' : '' };
        self.files[name]['txt'] += data + '\n'
        self.files[name]['cnt'] += 1
        if self.files[name]['cnt'] > 100 : 
            self.write(name)

    def write_all(self) :
        for name in self.files :
            self.write(name)
        return [name for name in self.files]

# ---------------------------------------------------------------------------- #
# Args parse

try:
    argx = {
        '-f' : sys.argv[1],
        '-b' : '/dev/null',
        '-t' : multiprocessing.cpu_count(),
        '-c' : '1',
        '-m' : '0.01',
        '-o' : 'output_' + datetime.now().strftime("%Y-%m-%d_%H-%M-%S") 
    }
    first = -1 * len(sys.argv)%2
    for i in range(0, len(sys.argv) - 1, 2) :
        k, v = [sys.argv[i + first], sys.argv[i + first + 1]]
        if k in argx : argx[k] = v
    argx['-m'] = float(argx['-m'])
except:
    fail('Аргументы где?')

if not os.path.exists(argx['-f']) :
    fail('vcf file not found (' + argx['-f'] + ')')

if not os.path.exists(argx['-o']):
    os.makedirs(argx['-o'])

if not os.path.exists(argx['-o'] + '/data/'):
    os.makedirs(argx['-o'] + '/data/')

fx = Write(argx['-o'] + '/data/')

# ---------------------------------------------------------------------------- #
# Get VCF data
# CHROM POS ID REF ALT QUAL FILTER INFO FORMAT sample1 sample2 sample3 ...
def zygnum(v):
    v = v.split(':')[0]
    z = {'1/1' : 1, '0/1' : 2, '0/0' : 3, './1' : 4, './.' : 5}
    return z[v] if v in z else 0

params = 9
vcf_data = {}
vcf_header = False

f = open(argx['-f'])
for line in f:
    line = line.replace('\n', '').split('\t')

    # Get header
    if line[0][0] == '#' :
        if line[0][1] != '#' and not vcf_header : vcf_header = line[params:]
        continue

    # Get lines
    if not vcf_header : continue

    # CHR POS REF ALT sample1_zyg sample2_zyg sample3_zyg ...
    zyg = array('B', [zygnum(line[n + params]) for n,s in enumerate(vcf_header)])
    for alt in line[4].split(',') :
        key = (':').join([line[0], line[1], line[3], alt])
        vcf_data[key] = zyg

f.close()

#print [len(vcf_data), i]
#sys.exit(0)

# ---------------------------------------------------------------------------- #
dir = os.path.dirname(os.path.realpath(__file__)) + '/'
# zyg = ['X', '1/1', '0/1', '0/0', './1', './.'];

sdf_plus = {}
with open(dir + 'data/sdf_plus.csv') as sdfp:
    next(sdfp)
    for line in sdfp:
        # Chromosome, Position, Ref, Alt, rsID,
        # Gene, Annotated_Type, Real_Type, 1000G_AF, ExAC_AF,
        # ESP_AF, PROVEAN_score, PROVEAN_prediction, Polyphen_score, Polyphen_prediction,
        # SIFT_score, SIFT_prediction, Relative,Gained
        e = line.replace('\n', '').split(',')
        skey = (':').join([e[0], e[17], e[18], e[19], e[20]])
        sdf_plus[skey] = e[0:7] + [''] + e[7:11] + ['',''] + e[11:]

with open(dir + 'data/sdf.csv') as sdf:
    next(sdf)
    for line in sdf:
        # Chromosome, Position, Ref, Alt, ID,
        # Gene, Type, Substitution, UniProt_ID, 1000G_AF,
        # ExAC_AF, ESP_AF, COSMIC_count, ClinVar, PROVEAN_score,
        # PROVEAN_prediction, Polyphen_score, Polyphen_prediction, SIFT_score, SIFT_prediction,
        # Coding
        e = line.replace('\n', '').split(',')
        if len(e) < 20 : continue
        e.append('?/?')

        # 0. Только выбранный класс интересует нас (only coding? = 1)
        if argx['-c'] == '1' and e[20] == 'NO' : continue
        
        maxafs = max(float(e[9]), float(e[10]), float(e[11]))

        # -. Поиск в файле юзера
        key = (':').join(e[0:4])
        zyg = vcf_data[key] if key in vcf_data else -1
        afs = True if maxafs > argx['-m'] else False

        for n, h in enumerate(vcf_header) :
            z = (-1) if zyg == -1 else zyg[n]
            n = str(n)

            # => 3,4 таблица
            if z == 1 :
                # => 3,4 таблица
                s0 = (':').join(e[0:4]) + ':0'
                if s0 in sdf_plus : 
                    t_key = (':').join(sdf_plus[s0][0:4])
                    if t_key in vcf_data : fx.add('tbl.'+n+'.t3', '0\t' + (',').join(sdf_plus[s0]))
                    
                s1 = (':').join(e[0:4]) + ':1'
                if s1 in sdf_plus : 
                    t_key = (':').join(sdf_plus[s1][0:4])
                    if t_key in vcf_data : fx.add('tbl.'+n+'.t4', '0\t' + (',').join(sdf_plus[s1]))

            if afs : continue

            # => 2 таблица

            # - Что находится в файле юзера в зиготности 1/1
            # сохраняется в список 2 (оставляем зиготность 1/1) 
            if z == 1 :
                e[21] = '1/1'
                fx.add('tbl.'+n+'.t2', str(maxafs) + '\t' + (',').join(e))

            # => 1 таблицa

            # - Что не находится в файле юзера сохраняется в список 1 c генотипом ./.
            if z == -1 :
                e[21] = './.'
                fx.add('tbl.'+n+'.t1', str(maxafs) + '\t' + (',').join(e))

            # - Что находится в файле юзера в зиготности 0/0
            # сохраняется в список 1 (изменяем зиготность 1/1, помечаем)
            if z == 3 :
                e[21] = '1/1'
                fx.add('tbl.'+n+'.t1', str(maxafs) + '\t' + (',').join(e))

            # - Что находится в файле юзера в зиготности 0/1
            # сохраняется в список 1 (оставляем зиготность 0/1)
            if z == 2 :
                e[21] = '0/1'
                fx.add('tbl.'+n+'.t1', str(maxafs) + '\t' + (',').join(e))

names = fx.write_all()
for name in names :
    f = argx['-o'] + '/data/' + name
    p = Popen('cat ' + f + ' | sort | awk {\'print $2\'} > /tmp/' + name + ' && cat /tmp/' + name + ' > ' + f, shell=True)
    p.wait()


# ---------------------------------------------------------------------------- #
# Make HTML report

css = '<style>' + open(dir + 'web/client/style.css', "r").read() + '</style>'
ico = 'data:image/x-icon;base64,' + base64.b64encode(open(dir + 'web/client/favicon.ico', "r").read())
jsx = '<script>' + open(dir + 'web/client/app.js', "r").read() + '</script>'

template = open(dir + 'web/index.html', "r").read()
template = template.replace('<link rel="stylesheet" href="client/style.css">', css)
template = template.replace('<script src="client/app.js"></script>', jsx)
template = template.replace('client/favicon.ico', ico)

for n, sample in enumerate(vcf_header) :
    data = { 'samples' : vcf_header, 'table' : {}, 'current' : sample, 'counts' : [0,0,0,0] }
    for i in ['1','2','3','4'] :
        try :
            data['table'][i] = open(argx['-o'] + '/data/tbl.' + str(n) + '.t' + i, "r").read()
        except :
            data['table'][i] = ""
        data['counts'][int(i) - 1] = len(data['table'][i].split("\n")) - 1

    r = open(argx['-o'] + '/sample_' + str(n) + '.html', 'w+')
    r.write(template.replace('/* data-local */', 'var local_data = ' + json.dumps(data)))
    r.close()

print "Done"



