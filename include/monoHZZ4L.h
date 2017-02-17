#ifndef MONOHZZ4L_H
#define MONOHZZ4L_H

#include <string>

struct monoHZZ4L
{
  static
  void analysis(std::string inputFile,
		int numberEvents=-1,    // all events
		double luminosity=1,    // 1/fb
		double xsection=1);     // 1 fb
};
#endif
