#include <fstream>
#include <iostream>
#include <vector>
#include <cstring>
#include <cstdio>
#include <cstdlib>
using namespace std;

// functions defined below 'main' program
bool parsecmd(int, char**,bool&, int&, int&, bool&, int&, string&);
void smooth(vector<double> &y);
void background(vector<double> y, vector<double> &bg, int w, int SMOOTH);
void printResults(vector<double>,vector<double>,vector<double>,string,string);
void printUsage(const char*);
void plot(string,string,int);

//
// Read in an ASCII spectrum and
// determine its 'pseudo-continuum'
//
int main( int argc, char **argv ) {

  bool BXCR, PLT;
  int smo, ssm, del;
  string inFile;
  // check input info
  if( !parsecmd(argc,argv,BXCR,smo,ssm,PLT,del,inFile) ) {
    printUsage(argv[0]);
    return 1;
  }
  else { // print info
    cout << endl << "Running with the following options:" << endl
	 << (BXCR ? "\tdo " : "\tdon't ") << "boxcar smooth" << endl
	 << "\tSNIP smoothing = " << smo << endl
	 << "\tsimultaneous smoothing to " << ssm*2+1 << " pixels" << endl
	 << (PLT ? "\t" : "\tnot ") << "plotting output" << endl;
  }
  // read in spectrum
  ifstream input; input.open(inFile.c_str());
  if( !input ) {
    cout << "File: " << inFile << " is not readable." << endl << endl;
    return 2;
  }
  else ;
  string dataDir, normDir;
  input >> dataDir >> normDir;
  cout << "Reading from: " << dataDir
       << "/ and writing to: " << normDir << "/\n";
  double x, y;
  vector<double> wave, flux, bkgr;
  string fname;
  char tmpname[256];
  ifstream instar;
  while( input >> fname ) {
    sprintf(tmpname,"%s/%s",dataDir.c_str(),fname.c_str());
    instar.open(tmpname);
    wave.clear(); flux.clear(); bkgr.clear();

    while( instar >> x >> y ) {
      wave.push_back(x);
      flux.push_back(y);
    }
    if( wave.size() == 0 ) { // check to see if we really have something
      cout << endl << "BAD INPUT: " << fname << " !!!" << endl << endl;
      return 3;
    }
    else if( BXCR ){ // boxcar smooth if wanted
      smooth(flux);
    }
    else ;

    // calculate the background
    background(flux,bkgr,ssm,smo);
    // create output
    printResults(wave,flux,bkgr,normDir,fname);

    if( PLT ) { // plot if wanted
      //system("gnuplot showBG.gps");
      plot(normDir,fname,del);
    }
    else ;

    instar.close();
  }
  input.close();

  return 0;
}

//
// Weighted boxcar smoothing
//
void smooth(vector<double> &y) {
  // smooth to 3 pixels
  int N = y.size();
  double c, *z;
  z = new double[N];
  double w[3] = {1., 3., 1.};
  for(int i=1; i<N-2; i++) {
    c=0.;
    for(int j=0; j<3; j++) c += w[j]*y[i-1+j];
    c /= 5.;
    z[i] = c;
  }
  for(int i=1; i<N-2; i++) y[i] = z[i];
  delete [] z;
}

//
// SNIP background estimation:
// decreasing window with smoothing
//
// REFERENCE:
// Algorithm D in Morhac, M. 2009, NIMA, 600, 478
//
void background(vector<double> y, vector<double> &bg,
		int BC, int SMOOTH) {
  if( bg.size() > 0 ) bg.clear(); else ;
  int N=y.size();
  int m=SMOOTH; // half the 'feature' size
  int w=BC; // simultaneous smoothing
  double a, b, c, *z;
  z = new double[N];
  for(int p=m; p>=1; p--) {
    for(int i=p; i<N-p; i++) {
      a = y[i];
      b = (y[i-p]+y[i+p])/2.;
      c = 0.;
      for(int j=i-w; j<=i+w; j++) c += y[j];
      c /= (2.*w+1.);
      z[i] = ( a>b ? c : b );
    }
    for(int i=p; i<N-p; i++) y[i] = z[i];
  }
  for(int i=0; i<N; i++) bg.push_back(y[i]);
  delete [] z;
}

//
// Print results to file 'norm.dat'
// columns:
//   1. wavelength
//   2. flux
//   3. background
//
void printResults(vector<double> w, vector<double> f, vector<double> b,
		  string dir, string fname) {
  //FILE *fout = fopen("norm.dat","w");
  char tmpname[256]; sprintf(tmpname,"%s/%s",dir.c_str(),fname.c_str());
  FILE *fout = fopen(tmpname,"w");
  for(int i=0; i<w.size(); i++) {
    fprintf(fout,"%10.5f  %12.5f  %12.5f\n",w[i],f[i],b[i]);
  }
  fclose(fout);
  cout << endl << "File: " << tmpname << " has been created." << endl;
}

//
// Plot results using gnuplot
//
void plot( string dir, string fname, int del ) {
  char cmd[256];
  sprintf(cmd,
	  "sed \'s/FILENAME/%s\\/%s/g\' template.gps |",
	  dir.c_str(),fname.c_str());
  sprintf(cmd,
	  "%s sed \'s/STAR/%s/\' |",
	  cmd,fname.c_str());
  if( del > 0 && del <= 10 )
    sprintf(cmd,"%s sed \'s/TIME/%d/\' > plot.gps",cmd,del);
  else if( del == -1 )
    sprintf(cmd,
	    "%s sed \'s/TIME/-1 \\\"Hit Enter to continue\\\"/\' > plot.gps",
	    cmd);
  else
    sprintf(cmd,"%s sed \'s/TIME/3/\' > plot.gps",cmd);
  system(cmd);
  system("gnuplot plot.gps");
  system("rm plot.gps");
}

//
// Figure out what to do
//
bool parsecmd(int argc, char **argv,
	      bool &bxcr, int &smo, int &ssm, bool &plt, int &del,
	      string &fname) {
  if( argc == 2 ) {
    fname = argv[1];
    bxcr = false;
    smo = 33;
    ssm = 1;
    plt = true;
    del = 3;
    return true;
  }
  else if( argc > 2 ) {
    bxcr = false;
    smo = 33;
    ssm = 1;
    plt = true;
    del = 3;
    fname = argv[argc-1];
    for(int i=1; i<argc; i++) {
      if( !strcmp(argv[i],"--boxcar") ) bxcr = true; else ;
      if( !strcmp(argv[i],"--smooth") ) smo = atoi(argv[i+1]); else ;
      if( !strcmp(argv[i],"--sismoo") ) ssm = atoi(argv[i+1]); else ;
      if( !strcmp(argv[i],"--no-plot") ) plt = false; else ;
      if( !strcmp(argv[i],"--delay") ) del = atoi(argv[i+1]); else ;
    }
    if( smo <= 0 || ssm < 0 || del > 10 || del == 0 ) {
      cout << endl << "Something wrong with input values!!!" << endl << endl;
      return false;
    }
    else
      return true;
  }
  else
    return false;
}

//
// Say how this works
//
void printUsage(const char *exe) {
  cout << endl << "Usage:" << endl
       << exe << " [option(s)] [filename]" << endl << endl
       << "Options:" << endl
       << "\t--boxcar    -> boxcar smooth the spectrum" << endl
       << "\t--smooth #  -> smoothing level for the SNIP algorithm" << endl
       << "\t--sismoo #  -> simultaneous smoothing "
       << "(boxcar smoothing during background estimation)" << endl
       << "\t--[no-]plot -> plot spectrum+background in gnuplot" << endl
       << "\t--delay #   -> delay interval for plots "
       << "(-1 = Enter to continue, 1-10 = seconds to wait)" << endl
       << endl
       << "Defaults:" << endl
       << "do NOT boxcar smooth (1 3 1)" << endl
       << "SNIP smoothing = 33" << endl
       << "simultaneously smoothing = 1 (3 pixels)" << endl
       << "show results in gnuplot" << endl
       << "plot delay of 3 secs" << endl << endl;
}
