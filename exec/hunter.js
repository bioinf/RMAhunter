const express = require('express');
const parser  = require('body-parser');
const exec    = require('child_process').exec;

const app     = express();
const server  = require('http').createServer(app);
const fs      = require('fs');
const genes   = require('./genes.js').e;

/* -------------------------------------------------------------------------- */

var port = parseInt(process.argv[2]);
server.listen(port);

app.use(express.static('./web'));
app.use(parser.urlencoded({extended : true, limit: '500mb'}));

/* -------------------------------------------------------------------------- */
app.post('/upload', (req, res) => {
    // Демонстрашка
    if (req.body.coding == 'demo'){
        setTimeout(function(){
            res.send(JSON.stringify(['demo', [15,1220,8]]));
        }, 2 * 1000);
        return ;
    }

    var coding = req.body.coding == 1 ? 'Y' : 'N';
    var noncalls = req.body.noncalls == 1 ? 'Y' : 'N';

    var maxafs = parseFloat(req.body.maxafs);
    if (isNaN(maxafs)) maxafs = 0.5;

	// Загрузка пользовательского ввода
	var vcf_str = req.body.vcf || '';
	var bed_str = req.body.bed || '';
	var zygosity = ['X', '1/1', '0/1', '0/0', './1', './.'];

    // Распаковаваем, что пришло
	var vcf = '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample1\n';
	vcf_str.split('!').map(function(chr_box){
		var blocks = chr_box.split('$'), last = 0;
		var chr = blocks[0];
		if (!chr || isNaN(chr)) return false;
		vcf[chr] = {};
		blocks.splice(1).map(function(block){
			var tmp = block.split('@'); // offet @ zyg @ ref @ alt,alt,alt,alt @ coverage
			var d3  = tmp[3].split(',');
			last += parseInt(tmp[0], 32);
			for (var d in d3) {
    			// chr (num), pos, ref, alt, zyg (num)
    			vcf += [chr, last, '', tmp[2], d3[d], '','','','GT:AD', zygosity[parseInt(tmp[1])] + ':' + tmp[4]].join('\t');
    			vcf += '\n';
			}
		});
	});

	var bed = '';
	bed_str.split('!').map(function(chr_box){
		var blocks = chr_box.split('$'), last = 0;
		var chr = blocks[0];
		if (!chr || isNaN(chr)) return false;
		bed[chr] = [];
		blocks.splice(1).map(function(block){
			var tmp = block.split('@');
			last += parseInt(tmp, 32);
    		bed += [chr, last, last + parseInt(tmp[1], 32)].join('\t');
    		bed += '\n';
		});
	});

    // Кладём локально во временные файлы
    var key = 'E' + Math.random().toString(36).substring(2).toUpperCase();
    fs.writeFileSync('/tmp/' + key + '.xvcf', vcf);
    fs.writeFileSync('/tmp/' + key + '.xbed', bed);
    
    // Обсчёт + разбиение на страницы
    var argv = ['./exec/app.sh', key, coding, noncalls, maxafs].join(' ');
    console.log(argv)
    
    exec(argv, function callback(error, stdout, stderr) {
        res.send(JSON.stringify([key, stdout.replace('\n', '').split(' ')]));
    });
});

app.post('/genes', (req, res) => {
	// Фильтруем файл по исходному ключу и выписываем новый в ответ
	var key_src = req.body.key || '';
	var gsx = JSON.parse(req.body.genes || '[]').filter(function(g){
    	return genes[g] && g != '';
	}).map(function(g){
    	return ',' + g + ','; 
	}).join('|');

    var new_key = 'EG' + Math.random().toString(36).substring(2).toUpperCase();
    var argv = ['./exec/filter.sh', key_src, '"' + gsx + '"', new_key].join(' ');

    exec(argv, function callback(error, stdout, stderr) {
        res.send(JSON.stringify([new_key, stdout.replace('\n', '').split(' ')]));
    });
});
