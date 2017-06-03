#pragma once
#include <stdlib.h>
#include <iomanip>

using namespace std;

template <typename T>
struct xarray
{
    xarray()
    : size(2)
    , count(0)
    , current(0)
    {
        data = new T[size];
    }
    
    ~xarray()
    {
        delete [] data;
    }
    
    void resize()
    {
        T * t = new T[size * 2];
        for (int k = 0; k < size; ++k) t[k] = data[k];
        delete [] data;
        data = t;
        size *= 2;
    }
    
    void append(T item)
    {
        if (size - 1 == count) resize();
        data[++count] = item;
    }
    
    T * data;
    int size;
    int count;
    int current;
};

struct intervals
{
    struct interval
    {
        interval()
        {}
        
        interval(short int chr, int from, int to)
        : chr(chr)
        , from(from)
        , to(to)
        {}
        
        interval& operator=(const interval& other)
        {
            chr = other.chr;
            from = other.from;
            to = other.to;
            return *this;
        }
        
        friend istream& operator >> ( istream& is, interval& e )
        {
            short int chr;
            int from, to;
            is >> chr >> from >> to;
            e = interval(chr, from, to);
            return is;
        }
        
        short int chr;
        int from, to;
    };
    
    intervals(string filename)
    {
        exist = false;
        chrx = new xarray<int> [26];
        ifstream raw(filename);
        interval line;
        while (raw >> line)
        {
            if (line.chr == 0) continue;
            chrx[line.chr].append(line.from);
            chrx[line.chr].append(line.to);
            exist = true;
        }
    }
    
    ~intervals()
    {
        delete [] chrx;
    }
    
    bool find(short int chr, int pos)
    {
        while (pos > chrx[chr].data[chrx[chr].current])
        {
            if (chrx[chr].current == chrx[chr].count - 1) break;
            chrx[chr].current += 1;
        }
        return chrx[chr].current % 2 ? false : true;
    }
    
    bool exist;
    xarray<int> * chrx;
};

struct sdfline
{
    sdfline()
    {}
    
    sdfline(string raw)
    : raw(raw)
    , chr(0)
    , maxafs(99)
    {
        // SRC:
        // Chromosome, Position, Ref, Alt, ID,
        // Gene, Type, Substitution, UniProt_ID, 1000G_AF,
        // ExAC_AF, ESP_AF, COSMIC_count, ClinVar, PROVEAN_score,
        // PROVEAN_prediction, Polyphen_score, Polyphen_prediction, SIFT_score, SIFT_prediction,
        // Coding
        
        // SRC-2:
        // Chromosome, Position, Ref, Alt, rsID,
        // Gene, Annotated_Type, Real_Type, 1000G_AF, ExAC_AF,
        // ESP_AF, PROVEAN_score, PROVEAN_prediction, Polyphen_score, Polyphen_prediction,
        // SIFT_score, SIFT_prediction, Relative,Gained
        
        unsigned short int k = 0;
        e[0] = "";
        for (unsigned short int x = 0; x < raw.length(); ++x)
        {
            if (raw[x] == ',')
            {
                e[++k] = "";
                continue;
            }
            e[k] += raw[x];
        }
        
        
        try {
            if (e[0] == "X") chr = 23;
            if (e[0] == "Y") chr = 24;
            if (e[0] == "Z") chr = 25;
            if (chr  == 0  ) chr = stoi(e[0]);
            pos = stoi(e[1]);
            
            maxafs = atof(e[9].c_str());
            maxafs = max(maxafs, atof(e[10].c_str()));
            maxafs = max(maxafs, atof(e[11].c_str()));
        } catch (...) {}
    }
    
    friend istream& operator >> ( istream& is, sdfline& e )
    {
        string raw;
        is >> raw;
        e = sdfline(raw);
        return is;
    }
    
    friend ostream& operator * ( ostream& out, const sdfline& h )
    {
        out << "0 ";
        for (unsigned short int i = 0; i < 23; ++i)
        {
            if (i == 8) out << ",";
            if (i == 11) out << ",,";
            out << (i == 18 ? "1/1" : h.e[i]) << (i == 22 ? "\n" : ",");
        }
        return out;
    }
    
    
    friend ostream& operator << ( ostream& out, const sdfline& h )
    {
        // first column (sorting)
        out << fixed << h.maxafs;
        out << string(2 - h.e[0].length(), '0') << h.e[0];
        out << string(9 - h.e[1].length(), '0') << h.e[1];
        out << " ";
        
        // OUTPUT:
        for (unsigned short int i = 0; i < 23; ++i) out << h.e[i] << (i == 22 ? "\n" : ",");
        return out;
    }
    
    string raw;
    short int chr;
    int pos;
    double maxafs;
    string e[23];
};

struct vcfdata
{
    struct vcfline
    {
        vcfline()
        {
        }
        
        vcfline(short int chr, int pos, string ref, string alt, short int zyg)
        : chr(chr), pos(pos), ref(ref), alt(alt), zyg(zyg)
        {
        }
        
        vcfline& operator=(const vcfline& other)
        {
            chr = other.chr;
            pos = other.pos;
            ref = other.ref;
            alt = other.alt;
            zyg = other.zyg;
            return *this;
        }
        
        friend istream& operator >> ( istream& is, vcfline& e )
        {
            short int chr, zyg;
            string ref, alt;
            int pos;
            is >> chr >> pos >> ref >> alt >> zyg;
            e = vcfline(chr, pos, ref, alt, zyg);
            return is;
        }
        
        friend ostream& operator << ( ostream& out, const vcfline& e )
        {
            out << e.chr << "\t" << e.pos << "\t" << e.ref << "\t"
            << e.alt << "\t" << e.pos << "\t" << e.zyg << "\n";
            return out;
        }
        
        short int chr;
        int pos;
        string ref;
        string alt;
        short int zyg;
    };
    
    vcfdata(string filename, int sqr)
    : sqr(sqr)
    {
        data = new xarray<vcfline> [sqr];
        ifstream raw(filename);
        vcfline line;
        while (raw >> line)
        {
            if (line.chr == 0) continue;
            data[line.pos % sqr].append(line);
        }
    }
    
    ~vcfdata()
    {
        delete [] data;
    }
    
    int get_p(short int chr, int pos, string ref, string alt)
    {
        int key = pos % sqr;
        for (int p = 0; p < data[key].count; ++p)
        {
            if (data[key].data[p].chr != chr) continue;
            if (data[key].data[p].pos != pos) continue;
            if (data[key].data[p].ref != ref) continue;
            if (data[key].data[p].alt != alt) continue;
            return p;
        }
        return -1;
    }
    
    short int zyg(sdfline& l)
    {
        int p = get_p(l.chr, l.pos, l.e[2], l.e[3]);
        return (p != -1) ? data[l.pos % sqr].data[p].zyg : -1;
    }
    
    int sqr;
    xarray<vcfline> * data;
};



