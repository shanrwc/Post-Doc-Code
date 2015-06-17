#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cstdlib>
#include <vector>

using namespace std;


void mergeTextFiles(int argc, const char* argv[]){
  
  //Create output file
  cout<<"Merged file saved to current directory. Extension used '.txt' by default."<<endl;
  fstream myOutfile("mergedFile.txt", std::ios_base::out | std::ios_base::trunc); 
  
  int evt = 1;
  //open input files one by one
  ifstream myInfile;
  string line;
  for (int i = 1; i<argc; i++){

    cout<<"Accessing file "<<argv[i]<<endl;
    myInfile.open(argv[i]);
  
    if (!myInfile) {
      cout<<"FILE NOT FOUND! Skipping File!"<<endl;
      continue;
    }
  
    while (!myInfile.eof()){
      line.clear();
      getline(myInfile, line);
      if(line != "" ) myOutfile<<line<<"\n"; //write to output file
      evt++;
    }
    myInfile.close();    
  }
  
  myOutfile.close();
  cout<<"Merging text files completed."<<endl;
  
}

    

    
#ifndef __CINT__
    
int main ( int argc, const char* argv[] ) {

  mergeTextFiles( argc,  argv) ;
  
  return 0 ;  
}

#endif
