#ifndef POLYLINEENCODER_H
#define	POLYLINEENCODER_H

#include <vector>
#include <string>
#include <cmath>
#include <memory>

namespace polyline {

class PolylineEncoder {
public:
    explicit PolylineEncoder(int numLevels=18, 
                             int zoomFactor=2, 
                             double threshold=0.00001, 
                             bool forceEndpoints=true,
                             const int& precision = 5);
    virtual ~PolylineEncoder();

    std::unique_ptr<std::pair<std::string, std::string> > dpEncode(std::vector<std::pair<double, double> >& points);
    std::string encodePoints(std::vector<std::pair<double,double> >& points);//encode only points and not levels; do not simplify points
    int getNumLevels() { return numLevels; }
    int getZoomFactor() { return zoomFactor; }
    
private:
    int numLevels;
    int zoomFactor;
    double threshold;
    bool forceEndpoints;
    double *zoomLevelBreaks;
    int precision;
    
    void _buildZoomLevelBreaks();
    double distance(const std::pair<double,double>& p0, const std::pair<double,double>& p1, const std::pair<double,double>& p2) const;
    inline int floorPrecision(double coordinate, int precision) const;
    inline int roundPrecision(double coordinate, int precision) const;
    inline int quick_pow10(int n) const;
    std::string encodeSignedNumber(int num) const;
    std::string encodeNumber(int num) const;
    int computeLevel(const double absMaxDist) const;
    std::unique_ptr<std::pair<std::string, std::string> > encode(std::vector<std::pair<double,double> >& points, const double dists[], double absMaxDist);
 
    // if we want to allow copying we need to implement these so we copy zoomLevelBreaks
    PolylineEncoder(const PolylineEncoder &rhs);
    PolylineEncoder & operator = (const PolylineEncoder &rhs);
};


int PolylineEncoder::quick_pow10(int n) const{
    static int pow10[10] = {
        1, 10, 100, 1000, 10000, 
        100000, 1000000, 10000000, 100000000, 1000000000
    };

    return pow10[n]; 
}

int PolylineEncoder::floorPrecision(double coordinate, int precision) const{
    return static_cast<int> (floor(coordinate * quick_pow10(precision)));
}

int PolylineEncoder::roundPrecision(double coordinate, int precision) const{
    return static_cast<int> (round(coordinate * quick_pow10(precision)));
}

} ///namespace polyline

#endif	/* POLYLINEENCODER_H */

