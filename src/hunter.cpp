#include <iostream>
#include <fstream>
#include <unistd.h>
#include <stdlib.h>
#include <unordered_map>
#include "hunter.h"

sdfline s;

int main(int argc, char *argv[])
{
    using namespace std;
    
    string xvcf = argv[1]; // .vcf (processed) filename
    string xbed = argv[2]; // .bed (processed) filename
    
    intervals bed(xbed);
    vcfdata vcf(xvcf, 1000); // 1000 это декомпозиция для быстрого поиска
    
    string tpx = argv[3]; // ? only coding
    string out = argv[5]; // output prefix
    string sdf = argv[6]; // supplementary data file
    string sec = argv[7]; // supplementary data file2
    
    double afs = strtod(argv[4], NULL); // max ref afs
    
    ofstream tbl1(out + ".t1");
    ofstream tbl2(out + ".t2");
    ofstream tbl3(out + ".t3");
    ofstream tbl4(out + ".t4");

    unordered_map<string, xarray<sdfline> *> second;
    ifstream ssec(sec);
    while (ssec >> s)
    {
        string key = s.e[0] + ":" + s.e[17] + ":" + s.e[18] + "->" + s.e[19];
        auto search = second.find(key);
        if (search == second.end()) second[key] = new xarray<sdfline>;
        second[key]->append(s);
    }
    
    ifstream sdfs(sdf);
    while (sdfs >> s)
    {
        // Только выбранный класс интересует нас (only coding? = 1)
        // tpx = '1','0'
        if (s.chr == 0) continue;
        if (tpx[0] == '1' && s.e[20] == "NO") continue;

        // Проверка, есть ли точка в интервалах BED
        if (bed.exist && !bed.find(s.chr, s.pos)) continue;

        // Поиск в файле юзера
        // zygosity = ['X', '1/1', '0/1', '0/0', './1', './.'];
        short int zyg = vcf.zyg(s);

        // 4. (2 таблица)
        // Что находится в файле юзера в зиготности 1/1
        // сохраняется в список 2 (оставляем зиготность 1/1)
        if (zyg == 1)
        {
            s.e[21] = "1/1";
            if (s.maxafs > afs) tbl2 << s;
            // 4.1 (3,4 таблица)
            string key = s.e[0] + ":" + s.e[1] + ":" + s.e[2] + "->" + s.e[3];
            auto search = second.find(key);
            if (search != second.end())
            {
                for (int xp = 0; xp < second[key]->count; ++xp)
                {
                    short int exist = vcf.zyg(second[key]->data[xp]);
                    if (exist != -1)
                    {
                        if (second[key]->data[xp].e[20] == "1") tbl3 * second[key]->data[xp];
                        if (second[key]->data[xp].e[20] == "0") tbl4 * second[key]->data[xp];
                    }
                }
            }
        }

        if (s.maxafs > afs) continue;

        // 1. (1 таблица)
        // Что не находится в файле юзера сохраняется в список 1 c генотипом ./.
        if (zyg == -1)
        {
            s.e[21] = "./.";
            tbl1 << s;
        }
        
        // 2. (1 таблица)
        // Что находится в файле юзера в зиготности 0/0
        // сохраняется в список 1 (изменяем зиготность 1/1, помечаем)
        if (zyg == 3)
        {
            s.e[21] = "1/1";
            s.e[22] = "*";
            tbl1 << s;
        }
        
        // 3. (1 таблица)
        // Что находится в файле юзера в зиготности 0/1
        // сохраняется в список 1 (оставляем зиготность 0/1)
        if (zyg == 2)
        {
            s.e[21] = "0/1";
            tbl1 << s;
            // 3.1 (3,4 таблица)
            // ...
            // tblX << s.raw << "\n";
        }
        
    }

    tbl1.close();
    tbl2.close();
    tbl3.close();
    tbl4.close();
    
    return 0;
}

