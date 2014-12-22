#include "include/PolylineEncoder.h"

#include <stack>
#include <sstream>

namespace polyline {


PolylineEncoder::PolylineEncoder(int numLevels, int zoomFactor, double threshold, bool forceEndpoints, const int& precision)
    : numLevels(numLevels), 
    zoomFactor(zoomFactor), 
    threshold(threshold), 
    forceEndpoints(forceEndpoints),
    zoomLevelBreaks(new double[numLevels]),
    precision(precision)
    {
        for (int i=0; i<numLevels; ++i) {
            zoomLevelBreaks[i] = threshold * pow(static_cast<double> (zoomFactor), numLevels - i - 1);
        }
    }

PolylineEncoder::~PolylineEncoder() {
    delete[] zoomLevelBreaks;
}

std::unique_ptr<std::pair<std::string, std::string> > PolylineEncoder::dpEncode(std::vector<std::pair<double, double> >& points) {
    int i, maxLoc = 0;
	typedef std::pair<int, size_t> pair_int_int;
	std::stack<pair_int_int> stack;
    double *dists = new double[points.size()];
	std::fill(&dists[0], &dists[points.size()], 0.0);
    double maxDist, absMaxDist = 0.0, temp = 0.0;

    // simplify using Douglas-Peucker
    if (points.size() > 2) {
		stack.push(pair_int_int(0, (points.size() - 1)));

        while (stack.size() > 0) {
			pair_int_int current = stack.top();
            stack.pop();
            maxDist = 0;

            for (i = current.first + 1; i < current.second; i++) {
                temp = distance(points[i], points[current.first], points[current.second]);
                if (temp > maxDist) {
                    maxDist = temp;
                    maxLoc = i;
                    if (maxDist > absMaxDist) {
                        absMaxDist = maxDist;
                    }
                }
            }
            if (maxDist > threshold) {
                dists[maxLoc] = maxDist;
				stack.push(pair_int_int(current.first, maxLoc));
				stack.push(pair_int_int(maxLoc, current.second));
            }
        }
    }

	std::unique_ptr<std::pair<std::string, std::string> > r = encode(points, dists, absMaxDist);
    delete[] dists;
    return r;
}

/**
 * distance(p0, p1, p2) computes the distance between the point p0 and the
 * segment [p1,p2]. This could probably be replaced with something that is a
 * bit more numerically stable.
 */
double PolylineEncoder::distance(const std::pair<double, double>& p0, const std::pair<double, double>& p1, const std::pair<double, double>& p2) const {
    double u, out = 0.0;

    if (p1.second == p2.second && p1.first == p2.first) {
        out = sqrt(pow(p2.second - p0.second, 2) + pow(p2.first - p0.first, 2));
    } else {
        u = ((p0.second - p1.second) * (p2.second - p1.second) + (p0.first - p1.first) * (p2.first - p1.first))
                / (pow(p2.second - p1.second, 2) + pow(p2.first - p1.first, 2));

        if (u <= 0) {
            out = sqrt(pow(p0.second - p1.second, 2) + pow(p0.first - p1.first, 2));
        }
        else if (u >= 1) {
            out = sqrt(pow(p0.second - p2.second, 2) + pow(p0.first - p2.first, 2));
        }
        else if (0 < u && u < 1) {
            out = sqrt(pow(p0.second - p1.second - u * (p2.second - p1.second), 2)
                    + pow(p0.first - p1.first - u * (p2.first - p1.first), 2));
        }
    }
    return out;
}

std::string PolylineEncoder::encodeSignedNumber(int num) const {
    int sgn_num = num << 1;
    if (num < 0) {
        sgn_num = ~(sgn_num);
    }
    return (encodeNumber(sgn_num));
}

std::string PolylineEncoder::encodeNumber(int num) const {
	std::ostringstream encodeString;

    while (num >= 0x20) {
        int nextValue = (0x20 | (num & 0x1f)) + 63;
        encodeString << (static_cast<char> (nextValue));
        num >>= 5;
    }

    num += 63;
    encodeString << (static_cast<char> (num));

    return encodeString.str();
}


std::unique_ptr<std::pair<std::string, std::string> > PolylineEncoder::encode(std::vector<std::pair<double, double> >& points, const double dists[], double absMaxDist) {
	std::ostringstream encodedLevels;
	std::ostringstream encodedPoints;

    int plat = 0;
    int plng = 0;

    if (forceEndpoints) {
        encodedLevels << (encodeNumber(numLevels - 1));
    } else {
        encodedLevels << (encodeNumber(numLevels - computeLevel(absMaxDist) - 1));
    }
    
    size_t n_points = points.size();
    for (size_t i=0; i<n_points; i++) {
        if ((i > 0) && (i < n_points-1) && (dists[i] != 0.0)) {
            encodedLevels << (encodeNumber(numLevels - computeLevel(dists[i]) - 1));
        }
        if (dists[i] != 0 || i == 0 || i == points.size() - 1) {
			std::pair<double, double> point = points[i];

            int lateMul = roundPrecision(point.second, precision);
            int lngeMul = roundPrecision(point.first, precision);

            int dlat = lateMul - plat;
            int dlng = lngeMul - plng;

            plat = lateMul;
            plng = lngeMul;

            encodedPoints << encodeSignedNumber(dlat);
            encodedPoints << encodeSignedNumber(dlng);
        }
    }
    
    if (forceEndpoints) {
        encodedLevels << encodeNumber(numLevels - 1);
    } else {
        encodedLevels << encodeNumber(numLevels - computeLevel(absMaxDist) - 1);
    }
    
	std::unique_ptr<std::pair<std::string, std::string> > r(new std::pair<std::string, std::string>);
    r->first = encodedPoints.str();
    r->second = encodedLevels.str();
    return r;
}

/**
 * This computes the appropriate zoom level of a point in terms of it's
 * distance from the relevant segment in the DP algorithm. Could be done in
 * terms of a logarithm, but this approach makes it a bit easier to ensure
 * that the level is not too large.
 */
int PolylineEncoder::computeLevel(const double absMaxDist) const {
    int lev = 0;
    if (absMaxDist > threshold) {
        while (absMaxDist < zoomLevelBreaks[lev]) {
            lev++;
        }
    }
    return lev;
}

//encode only points and not levels; do not simplify points
std::string PolylineEncoder::encodePoints(std::vector<std::pair<double,double> >& points){
    std::ostringstream encodedPoints;

    int plat = 0;
    int plng = 0;
    size_t n_points = points.size();
    for (size_t i=0; i<n_points; ++i) {
            int lateMul = roundPrecision(points[i].second, precision);
            int lngeMul = roundPrecision(points[i].first, precision);

            int dlat = lateMul - plat;
            int dlng = lngeMul - plng;

            plat = lateMul;
            plng = lngeMul;

            encodedPoints << encodeSignedNumber(dlat);
            encodedPoints << encodeSignedNumber(dlng);
    }
    return encodedPoints.str();
}

} ///namespace polyline

