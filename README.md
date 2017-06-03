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
git clone https://github.com/bioinf/RMAhunter.git ./
gzip -d build/data/sdf_plus.csv.gz build/data/sdf.csv.gz
chmod +x build/exec/*, build/hunter.py
~~~

~~~
Usage:
  ./build/hunter.py [input vcf-file] [Optional arguments]

Optional arguments:
  -f  Path to `input.vcf` file
  -c  Report coding only [0 or 1]. Default: 1
  -m  Allelic frequency cutoff. Default: 0.01
  -o  Output dir name

Examples:
  ./build/hunter.py input.vcf
  ./build/hunter.py -f input.vcf -c 0 -m 0.05 -o results
~~~


## How to start a server with a web-version

Install Node.JS, npm, forever. Example (ubuntu):

~~~
curl -sL https://deb.nodesource.com/setup_6.x | sudo -E bash -
sudo apt install -y nodejs
sudo apt install npm
npm install forever -g
npm install
~~~

Starting web-server on port `8915`

~~~
nodejs ./build/web/hunter.js 8915 # for debug
forever start ./build/web/hunter.js 8915 # for production
~~~

If file `build/data/sdf.csv` has been updated, you need to create a file with a list of genes for the web version:

~~~
echo "exports.e = {" $(
  echo $(
   cat build/data/sdf.csv | \
    awk -F  "," {'print $6'} | sort | uniq | \
    awk '{print "\""$1"\":true"}'
  ) | sed 's/ /,/g'
) "}" > build/web/genes.js
~~~

### How to build from sources

~~~
g++ -Werror -Wall -std=c++11 src/hunter.cpp -o build/exec/hunter # for debug
g++ -g -std=c++11 src/hunter.cpp -o build/exec/hunter # for production

# Delete temporary files
rm -f build/web/results/*

# Running on test files
cp build/data/test.xvcf /tmp/demo.xvcf
cp build/data/test.xbed /tmp/demo.xbed
./build/exec/app.sh demo 1 0.1
~~~
