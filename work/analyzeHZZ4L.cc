//----------------------------------------------------------------------------
// make reduced ntuples from Delphes files for HZZ4L RGS analysis
//----------------------------------------------------------------------------
#include <iostream>
#include <string>
#include <stdlib.h>
#include <getopt.h>
#include "monoHZZ4L.h"
//----------------------------------------------------------------------------
using namespace std;
int main(int argc, char** argv)
{
  if ( argc < 2 )
    {
      cout << "Usage: " << endl;
      cout << "./mkRGSntuple -l lumi -x xection Delphe file" << endl;
      exit(0);
    }


  string inputFile("");
  
  char optchar=0;        // single-character command line option

  // Loop until there are no more entries in argv

  // ab:c: means that this program recognizes the following options
  // -a    without a value
  // -b=X  with value X
  // -c=Y  with value Y
  double lumi=1.0;
  double xsection = 1.0;
  while ( (optchar = getopt(argc, argv, "l:x:")) != -1)
    {
      switch (optchar)
        {
        case 'l':
          {
	    lumi = atof(optarg);
            break;
          }
        case 'x':
          {
	    xsection = atof(optarg);
            break;
          }
        default:
          {
	    cout << "\tgo boil your head!" << endl;
            exit(0);
          }
        }
      if ( optind < argc )
	{
	  inputFile = string(argv[optind]);
	  break;
	}
    }
  cout << "input File:    " << inputFile << endl;
  cout << "cross section: " << xsection << " fb" << endl;
  
  int nevents = -1;
  
  monoHZZ4L::analysis(inputFile, nevents, lumi, xsection);

  return 0;
}
