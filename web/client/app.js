"use strict";
var local_data = local_data || false;
var header = header || false;

function Request(url, data, after, method){
    var method = method || 'POST';
    var req = new XMLHttpRequest();
    req.open(method, url, true);
    req.setRequestHeader('Content-Type', 'application/x-www-form-urlencoded; charset=UTF-8');
    req.onreadystatechange = function(e){
        if (req.readyState == 4 && req.status == 200) if (after) (after)(e.target.response || e);
    };
    req.send(data);
}

function $(id){
    return document.getElementById(id);
}

function Tpl(id, data){
    var html = $(id).innerHTML;
    for (var e in data){
        var find = new RegExp("`" + e + "`", "g");
        html = html.replace(find, data[e] == undefined ? '' : data[e]);
    }
    return html;
}

function Log(str){
    $('log').innerHTML += str + "<br/>";
}

function Pred(str){
    if (str != 'Deleterious') return str;
    return '<span style="color:#e00">Deleterious</span>';
}

/* -------------------------------------------------------------------------- */
var BED = {};
var VCF = {};

var key;
var global_key;
var tables = [];

var genes = {};
var process_init = false;
var onpage = 30;

/* -------------------------------------------------------------------------- */
$('result-tables').innerHTML += Tpl('table12', {
    'id' : 1,
    'header' : 'False negative RMAs',
    'desc' : 'False negative RMAs are sites at which your sample has a reference allele in contrast to the majority of global population. Genotype string 1/1 indicates that a minor (reference) allele is present in the homozygous state.'
});
$('result-tables').innerHTML += Tpl('table12', {
    'id' : 2,
    'header' : 'False positive RMAs',
    'desc' : 'False positive RMAs are sites at which your sample has a non-reference allele, similar to the majority of the global population. Genotype string 0/0 indicates that a major (non-reference) allele is present in the homozygous state.'
});
$('result-tables').innerHTML += Tpl('table34', {
    'id' : 3,
    'header' : 'Misannotated variants in the RMA codon',
    'desc' : 'Misannotated variants in the RMA codons. These are variants that fall into the same codon with the RMA site at which an alternative allele is discovered, leading to variant type misannotation.'
});

var init = function(e){
    Log('VCF info uploaded!');

    var x = JSON.parse(e);
    var counts = x[1];
    key = x[0], global_key = x[0];

    $('input').style.display  = 'none';
    $('result').style.display = 'block';
    
    var make_table_rows = function(e, i){
        if (typeof(e) != "string") return [];
        return e.split('\n').map(function(line){
            if (line == '') return '';
            var c = line.split(',');
            if (c[2].length > 3) c[2] = c[2].length + 'bp';
            if (c[3].length > 3) c[3] = c[3].length + 'bp';
            if (c[4].substr(0,3) == '~rs' || c[4].substr(0,2) == 'rs') {
                c[4] = Tpl('ncbi', { rs : c[4].replace('~', '').substr(2), h : c[4] });
            }
            if (i > 2){
                /* RefAFs */
                c[8]  = parseFloat(c[8]).toFixed(4);
                c[9]  = parseFloat(c[9]).toFixed(4);
                c[10] = parseFloat(c[10]).toFixed(4);
                /* Predictions */
                c[11] = Pred(c[11]);
                c[13] = Pred(c[13]);
                c[15] = Pred(c[15]);
                return Tpl('row3', c);
            } else {
                /* RefAFs */
                c[9]  = parseFloat(c[9] ).toFixed(4);
                c[10] = parseFloat(c[10]).toFixed(4);
                c[11] = parseFloat(c[11]).toFixed(4);
                /* Predictions */
                c[15] = Pred(c[15]);
                c[17] = Pred(c[17]);
                c[19] = Pred(c[19]);
                return Tpl('row', c);
            }
        }).filter(function(e){
            return e != '';
        });
    };
    
    tables = counts.map(function(count, i){
        var t = $('t' + (++i));
        t.page = 1;
        t.count = count;
        t.initial_count = count;
        t.tbody = t.childNodes[3];
        t.data = t.tbody.childNodes[1].childNodes[3];
        t.load = function(set_page){
            $('ct' + i).innerHTML = t.count;
            $('dl' + i).style.display  = t.count > 0 ? 'block' : 'none';
            $('dl' + i).href = '/results/' + key + '/data/tbl.0.t' + i;
            
            if (set_page) t.page = set_page;

            if (t.count == 0) {
                t.data.innerHTML = Tpl('row-empty', {});
                t.content = 70 + 45 + 36 + 'px';
            }

            t.maxpages = Math.floor(t.count/onpage) + 1;
            
            var load_part = function(table_part){
                var x = make_table_rows(table_part, i);
                if (x.length != 0) {
                    t.data.innerHTML = x.join('');
                    t.content = x.length * 22 + 45 + 36 + 'px';
                } else {
                    t.data.innerHTML = Tpl('row-empty', {});
                    t.content = 70 + 45 + 36 + 'px';
                }
                $('pgI' + i).innerHTML = t.page;
                $('pgC' + i).innerHTML = t.maxpages;
                if (t.opened) t.tbody.style.height = t.content;
            };

            if (local_data) {
                var cn = (t.page - 1) * onpage;
                var e = local_data.table[i].split('\n').slice(cn, cn + onpage).join('\n');
                load_part(e, i);
                $('dl' + i).href = key + '.t' + i;
                return ;
            } else {
                var p = (t.page - 1) + '';
                if (p.length == 3) p = '0'   + p;
                if (p.length == 2) p = '00'  + p;
                if (p.length == 1) p = '000' + p;
                Request('/results/' + key + '/data/tbl.0.t' + i + '.p' + p, '', function(e){
                    load_part(e, i)
                }, 'GET');
            }
        };
        
        t.load();

        t.childNodes[1].addEventListener('click', function(e){
            if (!t.opened) {
                t.opened = true;
                t.classList.add('opened');
                t.load();
                t.tbody.style.height = t.content;
            } else {
                t.opened = false;
                t.classList.remove('opened');
                t.tbody.style.height = '0px';
            }
        });

        $('pgL' + i).addEventListener('click', function(e){
            t.page = 1;
            t.load();
        });
        $('pgl' + i).addEventListener('click', function(e){
            t.page--;
            if (t.page <= 0) t.page = 1;
            t.load();
        });
        $('pgr' + i).addEventListener('click', function(e){
            t.page++;
            if (t.page > t.maxpages) t.page = t.maxpages;
            t.load();
        });
        $('pgR' + i).addEventListener('click', function(e){
            t.page = t.maxpages;
            t.load();
        });
    
        return t;
    });
};

if (local_data) {
    if (header) {
        var opt = '<option>Select sample ...</option>';
        for (var i in header) opt += '<option value="samples/'+i+'.html">' + header[i] + '</option>';
        $('init-local').innerHTML += '<pre>' + log + '</pre>';
        $('init-local').innerHTML += '<select id="sampleID">' + opt + '</select>';
        $('input').style.display = 'none';
        
    } else {
        var opt = local_data.samples.map(function(e, i){
            var current = local_data.current == e ? ' selected ' : '';
            return '<option value="../samples/'+i+'.html" '+current+'>' + e + '</option>'
        }).join();
        $('top-pane').innerHTML = '<select id="sampleID">' + opt + '</select>';
        init(JSON.stringify(['../data/tbl.' + (local_data.samples.indexOf(local_data.current)), local_data.counts]));
    }

    $('sampleID').onchange = function(){ location.href = this.value; };
}

$('run').addEventListener('click', function(e) {
    if (process_init) return ;
    if (!Object.keys(VCF).length) return ;

    process_init = true;
    document.getElementsByClassName( 'main' )[0].classList.add('disable');
    Log('Hunting. Please wait!');

    var q = '';
    q += 'maxafs=' + parseFloat($('maxafs').value);
    q += '&coding=' + ($('coding').checked ? 1 : 0);
    q += '&noncalls=' + ($('non-calls').checked ? 1 : 0);

    var data = [];
    for (var chr in VCF) {
        var last = 0, chrbox = chr;
        for (var pos in VCF[chr]) {
            chrbox += '$' + (parseInt(pos - last)).toString(32);
            chrbox += '@' + VCF[chr][pos];
            /* chrbox = $ offet @ zyg @ ref @ alt  */
            last = pos;
        }
        data.push(chrbox);
    }
    q += '&vcf=' + (data.join('!'));

    data = [];
    for (var chr in BED) {
        var last = 0, chrbox = chr;
        for (var i in BED[chr]) {
            chrbox += '$' + (parseInt(BED[chr][i][0] - last)).toString(32);
            chrbox += '@' + (parseInt(BED[chr][i][1] - BED[chr][i][0])).toString(32);
            last = BED[chr][i][0];
        }
        data.push(chrbox);
    }
    q += '&bed=' + (data.join('!'));
    
    Request('/upload', q, function(e){
        init(e);
    });
});

$('prc_demo').addEventListener('click', function(e) {
    process_init = true;
    document.getElementsByClassName( 'main' )[0].classList.add('disable');
    Log('Hunting (demo). Please wait!');
    Request('/upload', 'coding=demo', function(e){
        init(e);
    });
});

$('gset-open').addEventListener('click', function(e) {
    e.preventDefault();
    $('gset').style.display = 'block';
    $('gset-area').value = Object.keys(genes).join('\n');
});

$('gset-close').addEventListener('click', function(e) {
    $('gset').style.display = 'none';
});

$('gset-save').addEventListener('click', function(e) {
    genes = {};
    var txt = $('gset-area').value.replace(/(?:\r\n|\r| |\n)/g, ',');
    txt.replace(/ /g, '').split(',').map(function(gene){
        if (gene != '') genes[gene] = true;
    });

    var cnt = Object.keys(genes).length;
    $('gset-count').innerHTML = cnt > 0 ? (' (' + cnt + ')') : '';
    $('gset').style.display = 'none';
    
    if (cnt == 0)
    {
        key = global_key;
        for (var i in tables)
        {
            tables[i].count = tables[i].initial_count;
            tables[i].load(1);
        }
        return ;
    }

    var q = 'key=' + global_key;
    q += '&genes=' + JSON.stringify(Object.keys(genes));
    Request('/genes', q, function(e){
        var x = JSON.parse(e);
        key = x[0];
        for (var i in tables){
            tables[i].count = x[1][i];
            tables[i].load(1);
        }
    });
});


var draggable = false, current = [500,200];

window.onmousemove = function(e){
    if (!draggable) return;
    $('gset').style.top  = current[1] - draggable[1] + e.clientY + 'px';
    $('gset').style.left = current[0] - draggable[0] + e.clientX + 'px';
};

window.onmouseup = function(e){
    if (!draggable) return;
    current = [
        current[0] - draggable[0] + e.clientX, 
        current[1] - draggable[1] + e.clientY
    ];
    draggable = false;
};

$('header-gset').onmousedown = function(e){
    draggable = [e.clientX, e.clientY];
};

/* Выбран BED -> предобработка, подготовка к отправке */
$('bed').addEventListener('click', function(e) {
    if (process_init) {
        e.preventDefault();
        return false;
    }
});

$('bed').addEventListener('change', function(evt) {
    var filename = evt.target.value.split( '\\' ).pop();
    Log('Loading .bed file: ' + filename);
    var reader = new FileReader();
    reader.onload = function(e) {
        BED = {};
        var lines = e.target.result.split("\n"), count = 0;
        for (var i in lines) {
            var t = lines[i].split("\t"), chr = t[0];
            if (chr.substr(0, 3) == 'chr') chr = chr.substr(3);
            if (chr == "") continue;
            if (!BED[chr]) BED[chr] = [];
            BED[chr].push([parseInt(t[1]), parseInt(t[2])]);
            count ++;
        }
        Log('Done. Total intervals: ' + count);
    };
    reader.readAsText(evt.target.files[0]);
}, false);

/* Выбран VCF -> предобработка, подготовка к отправке */ 
$('vcf').addEventListener('click', function(e) {
    if (process_init) {
        e.preventDefault();
        return false;
    }
});

$('vcf').addEventListener('change', function(evt) {
    if (process_init) return ;
    var filename = evt.target.value.split( '\\' ).pop();
    Log('Loading .vcf file: ' + filename);
    var reader = new FileReader();
    reader.onload = function(e) {
        VCF = {};
        var lines = e.target.result.split("\n"), count = 0;
        var zygosity = {'1/1' : 1, '0/1' : 2, '0/0' : 3, './1' : 4, './.' : 5};
        for (var i in lines) {
            if (lines[i] == "" || lines[i][0] == '#') continue;
            var t = lines[i].split("\t"), chr = t[0], pos = t[1];
            if (chr.substr(0, 3) == 'chr') chr = chr.substr(3);
            if (!VCF[chr]) VCF[chr] = {};
            if (!VCF[chr][pos]) count ++;

            var F = t[8].split(':');
            var V = t[9].split(':');

            VCF[chr][pos]  = zygosity[V[F.indexOf('GT')]] || '5';
            /* + REF + ALT */
            VCF[chr][pos] += '@' + t[3] + '@' + t[4];
            /* Покрытие */
            var cpos = F.indexOf('AD');
            VCF[chr][pos] += '@' + (cpos == -1 ? '0,0' : (V[cpos] || '0,0'));
        }
        Log('Done. Total lines: ' + count);
        if (Object.keys(VCF).length) {
            document.getElementsByClassName( 'pane' )[0].classList.remove('disable');
        }
    };
    reader.readAsText(evt.target.files[0]);
}, false);
