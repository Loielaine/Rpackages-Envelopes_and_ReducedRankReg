#ifndef _Read_Matrix_615_H_
#define _Read_Matrix_615_H_

#include<iostream>
#include<vector>
#include<string>
#include<cstdlib>
#include<climits>
#include<cmath>
#include <boost/tokenizer.hpp>
#include<fstream>
#include<boost/lexical_cast.hpp>
#include<Eigen/Dense>

using namespace std;
using namespace boost;
using namespace Eigen;


template <class T>
Matrix<T,Dynamic, Dynamic> readFromFile(const char* fileName) {
      vector<T> temp; 
      ifstream ifs(fileName);
      if ( ! ifs.is_open() ) {
            cerr<< "Cannot open file" <<fileName<< endl;
            abort();
       }
       string line;
       char_separator<char> sep(", \t \"");
       typedef tokenizer<char_separator<char> > wsTokenizer;
       temp.clear();
       int nr=0, nc=0, nc_check = 0;
       while(getline(ifs, line) ) {
            if (line[0]=='#') continue;
             wsTokenizer t(line,sep);
             nc_check = 0;
             for(wsTokenizer::iterator i=t.begin(); i !=t.end(); ++i) {
                 temp.push_back(lexical_cast<T>(i->c_str()));
                 //temp.push_back(lexical_cast<T>(*i));
                 if (nr==0) ++nc;
                 ++nc_check;
              }
             
             if (nc != nc_check ) {
              cerr<<"The input file is not rectabgle at line "<<nr<<endl;
              abort();
        }
            ++nr;
        }
        
        Map<Matrix<T, Dynamic, Dynamic,RowMajor> > mf(&temp[0],nr,nc);
        return mf;
}

#endif
