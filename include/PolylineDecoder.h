#ifndef POLYLINEDECODER_H
#define	POLYLINEDECODER_H

#include <string>
#include <vector>

namespace polyline {


class PolylineDecoder {

public:

    static void decodePoly(const std::string& encoded, std::vector<std::pair<double,double> >& decodedPoints, double pointPrecision);

};


} ///namespace polyline

#endif	/* POLYLINEDECODER_H */

