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

app.use(express.static(__dirname));
app.use(parser.urlencoded({extended : true, limit: '500mb'}));

/* -------------------------------------------------------------------------- */
app.post('/upload', (req, res) => {
    // Демонстрашка
    if (req.body.coding == 'demo')
    {
        setTimeout(function(){
            res.send(JSON.stringify(['demo', [208,2378,7,9]]));
        }, 2 * 1000);
        return ;
    }

    var coding = req.body.coding == 1 ? '1' : '0';
    var maxafs = parseFloat(req.body.maxafs);
    if (isNaN(maxafs)) maxafs = 0.5;

	// Загрузка пользовательского ввода
	var vcf_str = req.body.vcf || '';
	var bed_str = req.body.bed || '';

    // Распаковаваем, что пришло
	var vcf = '';
	vcf_str.split('!').map(function(chr_box){
		var blocks = chr_box.split('$'), last = 0;
		var chr = parseInt(blocks[0].replace('X', 23).replace('Y', 24).replace('M', 25))
		if (chr == 0 || isNaN(chr)) return false;
		vcf[chr] = {};
		blocks.splice(1).map(function(block){
			var tmp = block.split('@');
			var d3  = tmp[3].split(',');
			last += parseInt(tmp[0], 32);
			for (var d in d3){
    			vcf += [chr, last, tmp[2], d3[d], parseInt(tmp[1])].join('\t');
    			vcf += '\n';
			}
		});
	});

	var bed = '';
	bed_str.split('!').map(function(chr_box){
		var blocks = chr_box.split('$'), last = 0;
		var chr = parseInt(blocks[0].replace('X', 23).replace('Y', 24).replace('M', 25))
		if (chr == 0 || isNaN(chr)) return false;
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
    var argv = ['../exec/app.sh', key, coding, maxafs].join(' ');
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
    var argv = ['../exec/filter.sh', key_src, '"' + gsx + '"', new_key].join(' ');
    console.log(argv)
    exec(argv, function callback(error, stdout, stderr) {
        res.send(JSON.stringify([new_key, stdout.replace('\n', '').split(' ')]));
    });
});
