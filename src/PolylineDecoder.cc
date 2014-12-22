#include "include/PolylineDecoder.h"

#include <math.h>

namespace polyline {


void PolylineDecoder::decodePoly(const std::string& encoded, std::vector<std::pair<double,double > >& decodedPoints, double pointPrecision) {

    size_t index = 0, len = encoded.size();
    int lat = 0, lng = 0;

    while (index < len) {
        int b, shift = 0, result = 0;
        do {
            b = encoded[index++] - 63;
            result |= (b & 0x1f) << shift;
            shift += 5;
        } while (b >= 0x20);
        int dlat = ((result & 1) != 0 ? ~(result >> 1) : (result >> 1));
        lat += dlat;

        shift = 0;
        result = 0;
        do {
            b = encoded[index++] - 63;
            result |= (b & 0x1f) << shift;
            shift += 5;
        } while (b >= 0x20);
        int dlng = ((result & 1) != 0 ? ~(result >> 1) : (result >> 1));
        lng += dlng;

        double decodeFactor = pow(10, pointPrecision);
        decodedPoints.push_back(std::pair<double, double>(lat / decodeFactor, lng / decodeFactor));
    }
}

} ///namespace polyline
