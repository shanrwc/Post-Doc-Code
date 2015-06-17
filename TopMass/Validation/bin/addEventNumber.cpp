#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cstdlib>
#include <vector>

using namespace std;


void addEventNumber(int argc, const char* argv[]){
  
  //Create output file
  cout<<"Adding event number to "<<endl;

 

  int evt = 1 ;


  //open input files one by one
  ifstream myInfile;
  string line;
  for (int i = 1; i<argc; i++){

    int first_evt = evt;
    int old_hyp = 1;

    stringstream outName;
    outName<<argv[i]<<".evt";
    fstream myOutfile(outName.str().c_str(), std::ios_base::out | std::ios_base::trunc); 
    
    cout<<"Accessing file "<<argv[i]<<endl;
    myInfile.open(argv[i]);
  
    if (!myInfile) {
      cout<<"FILE NOT FOUND! Skipping File!"<<endl;
      continue;
    }
  
    while (!myInfile.eof()){
      line.clear();
      getline(myInfile, line);
      
      if(line != "" ) {
	int pos = line.find(".");
	int hyp = atoi(line.substr(0,pos).c_str()) + 0 ;
	//	int pos1 = line.rfind("0");
	//	line.erase(pos1, 1);
	if (hyp == old_hyp){

	  myOutfile<<line<<" "<<evt<<"\n"; //write to output file
	}else{
	  evt=first_evt;
	  myOutfile<<line<<" "<<evt<<"\n"; 
	  old_hyp=hyp;
	}
	evt++;
      }
    }
    myInfile.close();    
    myOutfile.close();
 
  }
}

    

    
#ifndef __CINT__
    
int main ( int argc, const char* argv[] ) {

  addEventNumber( argc,  argv) ;
  
  return 0 ;  
}

#endif
