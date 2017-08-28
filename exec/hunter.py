#!/usr/bin/env python
# -*- coding: utf-8 -*-
# ---------------------------------------------------------------------------- #

SDF  = 'RMA_Annotations_NoESP.csv'
SDFp = 'RMA_Neighbor_Variants_WithEffs.csv'

# ---------------------------------------------------------------------------- #

import sys, os, time, json

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
    print '  \033[1;37m-v\033[0;37m  Show log [N or Y]. Default: Y'
    print '  \033[1;37m-c\033[0;37m  Report coding only [N or Y]. Default: Y'
    print '  \033[1;37m-g\033[0;37m  Analyze specific gene set [Comma separated list of genes]'
    print '  \033[1;37m-m\033[0;37m  Allelic frequency cutoff. Default: 0.01'
    print '  \033[1;37m-o\033[0;37m  Output dir name'
    print '  \033[1;37m-z\033[0;37m  Show non-calls [N or Y]. Default: Y'
    print ''
    print 'Examples:'
    print '  \033[1;37m' + sys.argv[0] + '\033[0;37m input.vcf'
    print '  \033[1;37m' + sys.argv[0] + '\033[0;37m -f input.vcf -c 0 -m 0.05 -o results'
    print '  \033[1;37m' + sys.argv[0] + '\033[0;37m -f input.vcf -g TLX1NB,USP17L25,TCP11X2,SFRP1,CAP1'
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
        '-c' : 'Y',
        '-g' : '',
        '-z' : 'Y',
        '-v' : 'Y',
        '-r' : 'Y',
        '-m' : '0.01',
        '-o' : 'output_' + datetime.now().strftime("%Y-%m-%d_%H-%M-%S") 
    }
    first = -1 * len(sys.argv)%2
    for i in range(0, len(sys.argv) - 1, 2) :
        k, v = [sys.argv[i + first], sys.argv[i + first + 1]]
        if k in argx : argx[k] = v
    argx['-m'] = float(argx['-m'])
    genes = False
    if argx['-g'] != "All" and argx['-g'] != "" :
        genes = {}
        g = argx['-g'].split(',')
        for g_name in g :
            genes[g_name] = True
except:
    fail('In qua rationes? (Specify arguments. Please)')

if not os.path.exists(argx['-f']) :
    fail('VCF file not found (' + argx['-f'] + ')')

if os.path.exists(argx['-o']):
    fail('Output directory already exists (' + argx['-o'] + ')')

os.makedirs(argx['-o'] + '/data/')

fx = Write(argx['-o'] + '/data/')

# ---------------------------------------------------------------------------- #
# Log
class Log :
    def __init__(self) :
        self.data = ''
        self.space = 27

    def write(self, param, value = "") :
        if argx['-v'] != 'Y' : return
        space = (self.space - len(param)) * ' '
        param += (':' if value != "" else "")

        print '\033[1;37m' + param + space + '\033[0;37m' + str(value) + '\033[0m'
        self.data += param + space + str(value) + '\n'
    
    def export(self) :
        return json.dumps(self.data)

log = Log()    

log.write('Init time', datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
log.write('Input file', argx['-f'])
log.write('Show log info', ('Yes' if argx['-v'] == 'Y' else 'No'))
log.write('Report coding only', ('Yes' if argx['-c'] == 'Y' else 'No'))
log.write('Allelic frequency cutoff', str(argx['-m']))
log.write('Output dir name', argx['-o'])
log.write('Show non-calls', ('Yes' if argx['-z'] == 'Y' else 'No'))


# ---------------------------------------------------------------------------- #
# Get VCF data
# CHROM POS ID REF ALT QUAL FILTER INFO FORMAT sample1 sample2 sample3 ...
def zygnum(v):
    v = v.split(':')[0]
    z = {'1/1' : 1, '0/1' : 2, '0/0' : 3, './1' : 4, './.' : 5}
    return z[v] if v in z else 0

def coverage(format, v):
    F = format.split(':')
    V = v.split(':')
    p = [i for i in range(len(F)) if F[i] == 'AD']
    if len(p) != 0 and p[0] < len(V):
        return V[p[0]]
    p1 = [i for i in range(len(F)) if F[i] == 'RO']
    p2 = [i for i in range(len(F)) if F[i] == 'AO']
    if len(p1) * len(p2) != 0 and p1[0] < len(V) and p2[0] < len(V):
        return V[p1[0]] + ',' + V[p2[0]]
    return '0,0'
    

params = 9
vcf_data = {}
vcf_data_cov = {}
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
    if line[0][0:3] == 'chr' : line[0] = line[0][3:]

    # CHR POS REF ALT sample1_zyg sample2_zyg sample3_zyg ...
    zyg = array('B', [zygnum(line[n + params]) for n,s in enumerate(vcf_header)])
    cov = [coverage(line[8], line[n + params]) for n,s in enumerate(vcf_header)]

    for alt in line[4].split(',') :
        key = (':').join([line[0], line[1], line[3], alt])
        vcf_data[key] = zyg
        vcf_data_cov[key] = cov

f.close()

log.write('Total lines in VCF file', len(vcf_data))


# ---------------------------------------------------------------------------- #
BED = False
if os.path.exists(argx['-b']) :
    BED,c = {},0
    b = open(argx['-b'])
    for line in b:
        c += 1
        e = line.replace('\n', '').split('\t')
        if e[0] not in BED : BED[e[0]] = {}
        F,T = [int(e[1]), int(e[2])]
        if F not in BED[e[0]] : BED[e[0]][F] = T
    b.close()
    if c == 0 : BED = False
    

def inBED(chr, point):
    if chr not in BED : 
        return False
    for F in BED[chr] :
        if point < F : continue
        return point <= T
    return False


# ---------------------------------------------------------------------------- #
dir = os.path.dirname(os.path.realpath(__file__)) + '/'
cnt = [0,0,0,0]

sdf_plus = {}
with open(dir + '../data/' + SDFp) as sdfp:
    next(sdfp)
    for line in sdfp:
        # Chromosome, Position, Ref, Alt, rsID,
        # Gene, Annotated_Type, Real_Type, 1000G_AF, 
        # ExAC_AF, ESP_AF, PROVEAN_score, PROVEAN_prediction, Polyphen_score, Polyphen_prediction,
        # SIFT_score, SIFT_prediction, Relative_Position, Relative_Ref, Relative_Alt, Gained, Annotated_Eff, Real_Eff
        e = line.replace('\n', '').split(',')
        skey = (':').join([e[0], e[17], e[18], e[19]])
        sdf_plus[skey] = e

with open(dir + '../data/' + SDF) as sdf:
    next(sdf)
    for line in sdf:
        # Chromosome, Position, Ref, Alt, ID,
        # Gene, Type, Substitution, UniProt_ID, 1000G_AF,
        # ExAC_AF, ESP_AF, COSMIC_count, ClinVar, PROVEAN_score,
        # PROVEAN_prediction, Polyphen_score, Polyphen_prediction, SIFT_score, SIFT_prediction,
        # Coding
        e = line.replace('\n', '').split(',')
        if len(e) < 20 : continue
        if genes and e[5] not in genes : continue
        
        e.append('?/?')
        e.append('0,0')

        # 0. Только выбранный класс интересует нас (only coding? = Y)
        if argx['-c'] == 'Y' and e[20] == 'NO' : continue
        
        # 0. Только если есть в BED
        if BED and not inBED(e[0], int(e[1])) : continue
        
        maxafs = max(float(e[9]), float(e[10]), float(e[11]))

        # -. Поиск в файле юзера
        key = (':').join(e[0:4])
        zyg = vcf_data[key] if key in vcf_data else (-1)
        afs = True if maxafs > argx['-m'] else False
        
        index_file = {
            't1' : [],
            't2' : [],
            't3' : [],
            't1_line' : '',
            't2_line' : '',
            't3_line' : ''
        }

        for n, h in enumerate(vcf_header) :
            z = -1
            if zyg != -1 :
                e[22] = vcf_data_cov[key][n]
                z = zyg[n]

            # => 3,4 таблица - теперь просто таблица 3
            
            # → Что находится в файле юзера в зиготности 1/1 проверяем, 
            # есть ли в файле sdf_plus.csv по колонкам Relative_Position, Relative_Ref, Relative_Alt
            # Если есть, это потенциальное место ошибки, смотрим на колонки Position, Ref, Alt 
            # и ищем их в пользовательском файле, если находим — добавляем в таблицу 3
            ztxt = ['','1/1','0/1']
            if z == 1 or z == 2 :
                skey = (':').join(e[0:4])
                if skey in sdf_plus : 
                    t_key = (':').join(sdf_plus[skey][0:4])
                    if t_key in vcf_data : 
                        if vcf_data[t_key][n] == 1 or vcf_data[t_key][n] == 2:
                            index_file['t3'].append(n)
                            index_file['t3_line'] = '0\t' + (',').join(sdf_plus[skey] + [vcf_data_cov[t_key][n], vcf_data_cov[key][n], ztxt[vcf_data[t_key][n]], ztxt[vcf_data[key][n]]])
                            fx.add('tbl.'+str(n)+'.t3', index_file['t3_line'])
                            cnt[3] += 1
            n = str(n)

            if afs : continue

            # => 2 таблица

            # → Что находится в файле юзера в зиготности 1/1
            # сохраняется в список 2 (меняем зиготность на 0/0) 
            if z == 1 :
                e[21] = '0/0'
                index_file['t2'].append(n)
                index_file['t2_line'] = str(maxafs) + '\t' + (',').join(e)
                fx.add('tbl.'+n+'.t2', index_file['t2_line'])
                cnt[2] += 1

            # => 1 таблицa

            # → Что не находится в файле юзера сохраняется в список 1 c генотипом ./.
            if z == -1 :
                e[21] = './.'
                if argx['-z'] == 'Y' :
                    index_file['t1'].append(n)
                    index_file['t1_line'] = str(maxafs) + '\t' + (',').join(e)
                    fx.add('tbl.'+n+'.t1', index_file['t1_line'])
                    cnt[1] += 1

            # → Что находится в файле юзера в зиготности 0/0
            # сохраняется в список 1 (изменяем зиготность 1/1, помечаем (?))
            if z == 3 :
                e[21] = '1/1'
                index_file['t1'].append(n)
                index_file['t1_line'] = str(maxafs) + '\t' + (',').join(e)
                fx.add('tbl.'+n+'.t1', index_file['t1_line'])
                cnt[1] += 1

            # → Что находится в файле юзера в зиготности 0/1
            # сохраняется в список 1 (оставляем зиготность 0/1)
            if z == 2 :
                e[21] = '0/1'
                index_file['t1'].append(n)
                index_file['t1_line'] = str(maxafs) + '\t' + (',').join(e)
                fx.add('tbl.'+n+'.t1', index_file['t1_line'])
                cnt[1] += 1
        
        if len(index_file['t1']) > 0 :
            fx.add('tbl._.t1', index_file['t1_line'] + ',' + ('|').join(map(str, index_file['t1'])))
        if len(index_file['t2']) > 0 :
            fx.add('tbl._.t2', index_file['t2_line'] + ',' + ('|').join(map(str, index_file['t2'])))
        if len(index_file['t3']) > 0 :
            fx.add('tbl._.t3', index_file['t3_line'] + ',' + ('|').join(map(str, index_file['t3'])))

log.write('Total samples', len(vcf_header))
log.write('False negative RMAs',   cnt[1])
log.write('False positive RMAs',   cnt[2])
log.write('Misannotated variants', cnt[3])


names = fx.write_all()
for name in names :
    f = argx['-o'] + '/data/' + name
    p = Popen('cat ' + f + ' | sort | awk {\'print $2\'} > /tmp/' + name + ' && cat /tmp/' + name + ' > ' + f, shell=True)
    p.wait()


# ---------------------------------------------------------------------------- #
# Make HTML report
if argx['-r'] == 'Y' :
    os.makedirs(argx['-o'] + '/samples/')
    
    def det_data(n, sample) :
        data = { 'samples' : vcf_header, 'table' : {}, 'current' : sample, 'counts' : [0,0,0] }
        for i in ['1','2','3'] :
            try :
                data['table'][i] = open(argx['-o'] + '/data/tbl.' + str(n) + '.t' + i, "r").read()
            except :
                data['table'][i] = ""
            data['counts'][int(i) - 1] = len(data['table'][i].split("\n")) - 1
        return data

    template = open(dir + '../web/index.html', "r").read()

    for n, sample in enumerate(vcf_header) :
        r = open(argx['-o'] + '/samples/' + str(n) + '.html', 'w+')
        r.write(template.replace('/* data-local */', 'var local_data = ' + json.dumps(det_data(n,sample))))
        r.close()

    r = open(argx['-o'] + '/index.html', 'w+')
    r.write(template.replace('/* data-local */', 'var local_data = ' + json.dumps(det_data('_', '_'))))
    r.close()

log.write('Done')
