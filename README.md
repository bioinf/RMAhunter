# RMA Hunter

[http://rmahunter.bioinf.me/](http://rmahunter.bioinf.me/)  

This is RMA Hunter — a web-based tool to systematically analyze and correct 
reference minor alleles in variant calling data. The tool provides a complete 
list of all potentially interesting RMAs that are not called in the sample 
analyzed or found reference homo-/heterozygous, as well as all 
reference-synonymous variants called in the RMA loci. All variants are 
annotated with correct pathogenicity predictions and reference allele 
frequencies.

To start, please upload your VCF file and (for exome sequencing or target 
sequencing assays) a BED-file containing enrichment intervals. To analyze only 
specific genes of interest, please enter a list of genes (separate with comma 
or a newline) in the box provided.

**Please cite the tool as:**  
Barbitoff Y.A., Bezdvornykh I.V., Serebryakova E.A., Glotov A.S., Glotov O.S. 
and Predeus A.V.  
Systematic correction of reference minor alleles in clinical variant calling. 
(2016)


## Quick Start

### How to install & run local version

~~~
git clone https://github.com/bioinf/RMAhunter.git && cd RMAhunter
gzip -d data/sdf_plus.csv.gz data/sdf.csv.gz
chmod +x exec/*
~~~

~~~
Usage:
  ./exec/hunter.py [input vcf-file] [Optional arguments]

Optional arguments:
  -f  Path to `input.vcf` file
  -v  Show log [N or Y]. Default: Y
  -c  Report coding only [N or Y]. Default: Y
  -m  Allelic frequency cutoff. Default: 0.01
  -o  Output dir name
  -z  Show non-calls [N or Y]. Default: Y

Examples:
  ./exec/hunter.py input.vcf
  ./exec/hunter.py -f input.vcf -c 0 -m 0.05 -o results
~~~


## How to start a server with a web-version

Install Node.JS, npm, forever. Example (ubuntu):

~~~
curl -sL https://deb.nodesource.com/setup_6.x | sudo -E bash -
sudo apt install -y nodejs npm
npm install forever -g
npm install

# Make demo samples
cp data/demo.xvcf /tmp/demo.xvcf
cp data/demo.xbed /tmp/demo.xbed
./exec/app.sh demo Y 0.01

# Starting web-server on port `8915`
nodejs ./exec/hunter.js 8915 # for debug
forever start ./exec/hunter.js 8915 # for production
~~~

If file `data/sdf.csv` has been updated, you need to create a file with a list of genes for the web version:

~~~
echo "exports.e = {" $(
  echo $(
   cat data/sdf.csv | \
    awk -F  "," {'print $6'} | sort | uniq | \
    awk '{print "\""$1"\":true"}'
  ) | sed 's/ /,/g'
) "}" > exec/genes.js
~~~
