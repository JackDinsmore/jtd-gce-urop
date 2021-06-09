#include <iostream>
#include <string>
#include <vector>
#include <assert.h>
#include <fstream>
#include <sstream>
#include <thread>
#include <future>
#include <boost/math/special_functions/gamma.hpp>

typedef std::pair<double, double> CoordPair;
typedef std::pair<int, int> IndexPair;
typedef std::vector<std::vector<double>> DoubleVector;

//#define LOTS_OF_THREADS // Use on Erebus
/* When LOTS_OF_THREADS is defined, a target of 80-90 threads are made and run concurrently.
   Otherwise, four threads are made and run concurrently. */

#define SENSITIVITY_DIVISOR 10.0

#define pi 3.14159265358979323
#define ONE_PLUS_EPSILON 1.0000000001
#define max(a, b) (((a) > (b)) ? (a) : (b))

#define ALPHA 1.94
#define L_MIN 1.0e29
#define DIST_TO_CENTER 8.5 // kpc
#define CM_PER_KPC_SQUARED 9.523396e+42
#define FLUX_EXCESS 1.2953417255755896e-09  // All flux units are in ergs per second per square centimeter

#define FLUX_HIST_BINS 40
#define FLUX_HIST_LOW 1.0731612e-12
#define FLUX_HIST_HIGH 1.3719972e-09
#define ALPHA_HIST 1.94
#define L_MIN_HIST 1.0e29
#define L_MAX_HIST 1.0e35
#define L0_HIST 0.88e34
#define SIGMA_HIST 0.62
#define L0_PLOEG_HIST 1.6084e+32
#define SIGMA_PLOEG_HIST 0.7003
#define NPTF_L_BREAK_HIST 8.656487610122969e+33
#define NPTF_N_1_HIST -0.66
#define NPTF_N_2_HIST 18.2



std::array<double, FLUX_HIST_BINS> fluxValues, fluxWidths;

const double L_MIN_RANGE[2] = { 1.0e28, 1.0e34 };
const double L_MAX_RANGE[2] = {1.0e34, 1.0e38}; //{ 1.0e34, 1.0e36 };
const double ALPHA_RANGE[2] = { 1.1, 2.5 };
const double L0_RANGE[2] = {1.0e30, 2.0e36}; //{ 1.0e32, 2.0e34 };
const double SIGMA_RANGE[2] = { 0.001, 1 };

const double fermilabPaperPoint[2] = { 1e35, 1e29 };
const double logNormalPaperPoint[2] = { 0.88e34, 0.62 };
const double logNormalPloegPoint[2] = { 1.6084e+32, 0.7003 };

// Values for the computation of the integral of the NFW profile.
#define GAMMA 1.2
#define NFW_SCALE_DIST 20 // kpc
#define INTEGRAL_HALF_WIDTH (DIST_TO_CENTER * 0.9)// / 2 // kpc
#define NFW_INTEGRAL_STEPS 1000
#define PLOT_SIZE 50
#define FERMILAB_ALPHA 1.94

#define ROOT "C:/Users/goods/Dropbox (MIT)/GCE UROP/"
#define SENSITIVITY_PATH (ROOT "sensitivity/sensitivity.txt")
#define PLOEG_PATH (ROOT "luminosity-models-step/ploeg/data/disk.csv")
#define DISPLAY_SIZE 20.0
#define myAbs(a) ((a) > 0 ? (a) : -(a))

double nptfLBreak;
double nptfPremul = 0;
double nptfLMin = 0;

enum class VALUE {
    TOTAL_NUM,
    TOTAL_FLUX,
    SEEN_NUM,
    SEEN_FLUX,
};

enum class LUMINOSITY_FUNCTION {
    POWER_LAW,
    LOG_NORMAL,
    PLOEG,
    NPTF,
    POWER_LAW_ALPHA,
    ERROR,
};

enum class METHOD {
    COUNT,
    HIST,
    ERROR,
};

LUMINOSITY_FUNCTION luminosityFunction = LUMINOSITY_FUNCTION::POWER_LAW;

std::string sciNot(double value) {
    int power = log10(value);
    if (log10(value) < 0 && value != pow(10, power)) {
        power--;
    }
    if (myAbs(power) > 6) {
        return std::to_string(value / pow(10.0, power)) + "e" + std::to_string(power);
    }
    return std::to_string(value);
}

class SensitivityMap {
public:
    SensitivityMap() {
        std::ifstream sensitivityFile;
        sensitivityFile.open(SENSITIVITY_PATH);
        if (sensitivityFile.is_open()) {
            std::string line;
            while (std::getline(sensitivityFile, line)) {
                std::vector<double> v;
                std::stringstream ss(line);
                for (double d; ss >> d;) {
                    v.push_back(d / SENSITIVITY_DIVISOR);
                    if (ss.peek() == ',') {
                        ss.ignore();
                    }
                }
                fluxes.push_back(v);
            }
        }
        else {
            std::cout << "Sensitivity file not found." << std::endl;
            std::cin.get();
        }

        upperLeft[0] = int((1 - DISPLAY_SIZE / 90.0) * fluxes.size() / 2.0);
        upperLeft[1] = int((1 - DISPLAY_SIZE / 180.0) * fluxes[0].size() / 2.0);
        lowerRight[0] = int((1 + DISPLAY_SIZE / 90.0) * fluxes.size() / 2.0);
        lowerRight[1] = int((1 + DISPLAY_SIZE / 180.0) * fluxes[0].size() / 2.0);
        skyShape = { lowerRight[0] - upperLeft[0], lowerRight[1] - upperLeft[1] };
    }

public:
    CoordPair indexToLatLon(IndexPair xy) const {
        double deltaXFrac = 1 - xy.first / (skyShape.first / 2.0);
        double deltaYFrac = xy.second / (skyShape.second / 2.0) - 1;
        return { DISPLAY_SIZE * deltaXFrac * pi / 180.0, DISPLAY_SIZE * deltaYFrac * pi / 180.0 };
    }
    IndexPair latLonToIndex(CoordPair latLon) const {
        double deltaXFrac = latLon.first * 180.0 / pi / DISPLAY_SIZE;
        double deltaYFrac = latLon.second * 180.0 / pi / DISPLAY_SIZE;
        return { int((1 - deltaXFrac) * (skyShape.first / 2.0)), int((1 + deltaYFrac) * (skyShape.second / 2.0)) };
    }
    double getFluxThreshold(CoordPair latLon) const {
        IndexPair xy = latLonToIndex(latLon);
        return fluxes[xy.first + upperLeft[0]][xy.second + upperLeft[1]];
    }

public:
    IndexPair skyShape;
    int upperLeft[2];
    int lowerRight[2];
    DoubleVector fluxes;
};

class PloegData {
public:
    PloegData() {
        std::ifstream ploegFile;
        ploegFile.open(PLOEG_PATH);
        if (ploegFile.is_open()) {
            std::string line;
            while (std::getline(ploegFile, line)) {
                std::vector<double> v;
                std::stringstream ss(line);
                for (double d; ss >> d;) {
                    v.push_back(d); // L, dN/dlogL
                    if (ss.peek() == ',') {
                        ss.ignore();
                    }
                }
                xs.push_back(pow(10, v[0]));
                logxs.push_back(v[0]);
                ys.push_back(log10(exp(1)) / pow(10, v[0]) * v[1]);
            }
        }
        else {
            std::cout << "Ploeg data not found." << std::endl;
            std::cin.get();
        }

        normalization = integrate(pow(10, logxs[0]) * ONE_PLUS_EPSILON, false);
    }
    double operator()(double l) const {// l is the log10 of the luminosity
        if (l <= logxs[0]) { return ys[0]; }
        if (l >= logxs[logxs.size() - 1]) { return ys[ys.size() - 1]; }

        // Binary search for l:
        int leftIndex = 0;
        int rightIndex = logxs.size() - 1;
        for (;;) {
            int midIndex = (rightIndex + leftIndex) / 2;
            if (logxs[midIndex] > l) {
                rightIndex = midIndex;
            }
            else if (logxs[midIndex] < l) {
                leftIndex = midIndex;
            }
            else {
                return ys[midIndex];
            }
            if (rightIndex - leftIndex <= 1) {
                assert(leftIndex != rightIndex);
                return ys[leftIndex] + (ys[rightIndex] - ys[leftIndex]) * (l - logxs[leftIndex]) / (logxs[rightIndex] - logxs[leftIndex]);
            }
        }
    }

public:
    double integrate(double start, bool normalize=true) const {
        int i;
        for (i = 0; i < xs.size() && xs[i] < start; i++);
        assert(i < xs.size());
        assert(i > 0);
        double frac = (start - xs[i - 1]) / (xs[i] - xs[i - 1]);
        double startY = ys[i - 1] + (ys[i] - ys[i - 1]) * frac;
        double integral = 0.5 * (startY + ys[i]) * (xs[i] - start);
        for (; i < xs.size() - 1; i++) {
            integral += 0.5 * (ys[i] + ys[i + 1]) * (xs[i + 1] - xs[i]);
        }
        if (normalize) {
            return integral / normalization;
        }
        else {
            return integral;
        }

        /*double integral = 0;
        for (double logx = log10(start); logx < logxs[logxs.size() - 1]; logx += INTEGRATION_LOG_STEP) {
            // Trapezoidal integration
            double xLow = pow(10, logx);
            double xHigh = pow(10, logx + INTEGRATION_LOG_STEP);
            double yLow = operator()(logx);
            double yHigh = operator()(logx + INTEGRATION_LOG_STEP);
            integral += (yLow + yHigh) / 2 * (xHigh - xLow);
        }
        if (normalize) {
            return integral / normalization;
        }
        else {
            return integral;
        }*/
    }

    double lintegrate(double start, bool normalize = true) const {int i;
        for (i = 0; i < xs.size() && xs[i] < start; i++);
        assert(i < xs.size());
        assert(i > 0);
        double frac = (start - xs[i - 1]) / (xs[i] - xs[i - 1]);
        double startY = ys[i - 1] + (ys[i] - ys[i - 1]) * frac;
        double integral = 1 / 6.0 * (-start * start * (ys[i] + 2 * startY) + start * xs[i] * (startY - ys[i]) + xs[i] * xs[i] * (2 * ys[i] + startY));
        for (; i < xs.size() - 1; i++) {
            integral += 1 / 6.0 * (ys[i] * (-2 * xs[i] * xs[i] + xs[i] * xs[i + 1] + xs[i + 1] * xs[i + 1])
                + ys[i + 1] * (-xs[i] * xs[i] - xs[i] * xs[i + 1] + 2 * xs[i + 1] * xs[i + 1]));
        }
        if (normalize) {
            return integral / normalization;
        }
        else {
            return integral;
        }
        /*double integral = 0;
        for (double logx = log10(start); logx < logxs[logxs.size() - 1]; logx += INTEGRATION_LOG_STEP) {
            // Trapezoidal integration
            double xLow = pow(10, logx);
            double xHigh = pow(10, logx + INTEGRATION_LOG_STEP);
            double yLow = operator()(logx);
            double yHigh = operator()(logx + INTEGRATION_LOG_STEP);
            integral += (xLow * yLow + xHigh * yHigh) / 2 * (xHigh - xLow);
        }
        if (normalize) {
            return integral / normalization;
        }
        else {
            return integral;
        }*/
    }

public:
    std::vector<double> logxs, ys, xs;
    double normalization;
};

template<int N>
void addArrays(std::array<double, N>& result, std::array<double, N> summand) {
    for (int i = 0; i < N; i++) {
        result[i] += summand[i];
    }
}

const SensitivityMap thresholds;
const PloegData ploegData;

void setNPTFPremul(double n1, double n2) {
    nptfPremul = (n1 - n2) / (n1 - n2 - (1 - n2) * pow(nptfLMin / nptfLBreak, 1 - n1));
}


// ================= These functions characterize the luminosity function ============

// The luminosity function must be normalized.

double improvedGamma(double s, double x) {
    assert(x >= 0);
    while (s < 0) {
        return (improvedGamma(s+1, x) - pow(x, s) * exp(-x)) / s;
    }
    return boost::math::tgamma(s, x);
}

double integrate(double start, double arg1, double arg2) {// arg1, arg2; l0, sigma; nBelow, nAbove
    switch (luminosityFunction) {
    case LUMINOSITY_FUNCTION::POWER_LAW:
        return improvedGamma(1 - ALPHA, start / arg2) / improvedGamma(1 - ALPHA, arg1 / arg2);
    case LUMINOSITY_FUNCTION::LOG_NORMAL:
        return 0.5 * (1 - erf((log10(start)-log10(arg1))/(sqrt(2) * arg2)));
    case LUMINOSITY_FUNCTION::PLOEG:
        return ploegData.integrate(start);
    case LUMINOSITY_FUNCTION::NPTF:
        if (start < nptfLBreak) {
            return nptfPremul * (1 - pow(nptfLBreak / start, arg1 - 1) * (arg2 - 1) / (arg2 - arg1));
        }
        else {
            return nptfPremul * pow(nptfLBreak / start, arg2 - 1) * (1 - arg1) / (arg2 - arg1);
        }
    case LUMINOSITY_FUNCTION::POWER_LAW_ALPHA:
        return improvedGamma(1 - arg1, start / arg2) / improvedGamma(1 - arg1, L_MIN / arg2);
    default:
        std::cout << "Running integrate with an invalid luminosity function." << std::endl;
        std::cin.get();
    }
}

double lintegrate(double start, double arg1, double arg2) {// lMin, lMax; l0, sigma
    switch (luminosityFunction) {
    case LUMINOSITY_FUNCTION::POWER_LAW:
        if (start == NULL) { start = arg1; }
        return arg2 * improvedGamma(2 - ALPHA, start / arg2) / improvedGamma(1 - ALPHA, arg1 / arg2);
    case LUMINOSITY_FUNCTION::LOG_NORMAL:
        if (start == NULL) { start = 0; }
        return 0.5 * arg1 * exp(arg2 * arg2 * log(10)*log(10) / 2) * (1 - erf((log10(start) - log10(arg1) - arg2 * arg2 * log(10)) / (sqrt(2) * arg2)));
    case LUMINOSITY_FUNCTION::PLOEG:
        if (start == NULL) { start = pow(10, ploegData.logxs[0]) * ONE_PLUS_EPSILON; }
        return ploegData.lintegrate(start);
    case LUMINOSITY_FUNCTION::NPTF:
        if (start < nptfLBreak) {
            return nptfPremul * nptfLBreak * (1 - arg1) * (1-arg2) * (1 / ((arg1-2) * (arg2-2)) + pow(nptfLBreak / start, arg1 - 2) / ((arg1 - 2) * (arg1 - arg2)));
        }
        else {
        return nptfPremul * nptfLBreak * (1 - arg1) * (1 - arg2) * (pow(nptfLBreak / start, arg2 - 2) / ((arg2 - 2) * (arg1 - arg2)));
        }
    case LUMINOSITY_FUNCTION::POWER_LAW_ALPHA:
        if (start == NULL) { start = L_MIN; }
        return arg2 * improvedGamma(2 - arg1, start / arg2) / improvedGamma(1 - arg1, L_MIN / arg2);
    default:
        std::cout << "Running integrate with an invalid luminosity function." << std::endl;
        std::cin.get();
    }
}


// ================= Declare helper functions =======================

double numSeenFunc(double threshold, double arg1, double arg2) {// Returns(unscaled) number of visible pulsars
    //assert(threshold > arg1);
    return integrate(threshold, arg1, arg2);
}

double fluxSeenFunc(double threshold, double arg1, double arg2) {// Returns(unscaled) amount of luminosity visible
    //assert(threshold > arg1);
    return lintegrate(threshold, arg1, arg2);
}

double totalFluxFunc(double threshold, double arg1, double arg2) {// Returns(unscaled) total luminosity
    return lintegrate(NULL, arg1, arg2);
}


double getValueAtLatLon(CoordPair latLon, double fluxThreshold, VALUE value, double arg1, double arg2) {
    // returns the (unscaled) value for a specific lon and lat

    // Integrate value along the line of sight:
    const double deltaRadialDistance = INTEGRAL_HALF_WIDTH * 2 / NFW_INTEGRAL_STEPS;
    double radialDistance = DIST_TO_CENTER - INTEGRAL_HALF_WIDTH;
    double integral = 0;
    double distFromCenter, nfwSquaredValue, lThreshold;
    const double cosLat = cos(latLon.first);
    const double cosLon = cos(latLon.second);
    while (radialDistance < DIST_TO_CENTER + INTEGRAL_HALF_WIDTH) {
        distFromCenter = sqrt(radialDistance * radialDistance + DIST_TO_CENTER * DIST_TO_CENTER
            - 2 * DIST_TO_CENTER * radialDistance * cosLon * cosLat);
        nfwSquaredValue = pow(distFromCenter / NFW_SCALE_DIST, -GAMMA) * pow(1 + distFromCenter / NFW_SCALE_DIST, -3 + GAMMA);
        nfwSquaredValue = nfwSquaredValue * nfwSquaredValue;

        lThreshold = fluxThreshold * 4 * pi * radialDistance * radialDistance * CM_PER_KPC_SQUARED;

        switch (value) {
        case VALUE::TOTAL_NUM:
            integral += nfwSquaredValue * deltaRadialDistance * radialDistance * radialDistance * cosLat * CM_PER_KPC_SQUARED;
            break;
        case VALUE::SEEN_NUM:
            integral += nfwSquaredValue * deltaRadialDistance * radialDistance * radialDistance * cosLat * numSeenFunc(lThreshold, arg1, arg2) * CM_PER_KPC_SQUARED;
            break;
        case VALUE::TOTAL_FLUX:
            integral += nfwSquaredValue * deltaRadialDistance * cosLat / (4 * pi) * totalFluxFunc(lThreshold, arg1, arg2);
            break;
        case VALUE::SEEN_FLUX:
            integral += nfwSquaredValue * deltaRadialDistance * cosLat / (4 * pi) * fluxSeenFunc(lThreshold, arg1, arg2);
            break;
        }
        // Ignore the angular part since it introduces a constant cofactor, except for the cosine term which accounts for shinking pixels far from the equator.

        radialDistance += deltaRadialDistance;
    }

    return integral;
}

std::array<double, FLUX_HIST_BINS> getHistAtLatLon(CoordPair latLon, double fluxThreshold) {
    // returns the (unscaled) value for a specific lon and lat

    // Integrate value along the line of sight:
    const double deltaRadialDistance = INTEGRAL_HALF_WIDTH * 2 / NFW_INTEGRAL_STEPS;
    double radialDistance = DIST_TO_CENTER - INTEGRAL_HALF_WIDTH;
    std::array<double, FLUX_HIST_BINS> integral;
    integral.fill(0);
    double distFromCenter, nfwSquaredValue, lThreshold;
    const double cosLat = cos(latLon.first);
    const double cosLon = cos(latLon.second);
    while (radialDistance < DIST_TO_CENTER + INTEGRAL_HALF_WIDTH) {
        distFromCenter = sqrt(radialDistance * radialDistance + DIST_TO_CENTER * DIST_TO_CENTER
            - 2 * DIST_TO_CENTER * radialDistance * cosLon * cosLat);
        nfwSquaredValue = pow(distFromCenter / NFW_SCALE_DIST, -GAMMA) * pow(1 + distFromCenter / NFW_SCALE_DIST, -3 + GAMMA);
        nfwSquaredValue = nfwSquaredValue * nfwSquaredValue;

        lThreshold = fluxThreshold * 4 * pi * radialDistance * radialDistance * CM_PER_KPC_SQUARED;

        for (int i = 0; i < FLUX_HIST_BINS; i++) {
            double flux = fluxValues[i];
            double lum = flux * 4 * pi * radialDistance * radialDistance * CM_PER_KPC_SQUARED;
            double lumRange = fluxWidths[i] * 4 * pi * radialDistance * radialDistance * CM_PER_KPC_SQUARED;
            if (lum < lThreshold) {
                continue;
            }
            double prob = 0;
            switch (luminosityFunction) {
            case LUMINOSITY_FUNCTION::POWER_LAW:
            case LUMINOSITY_FUNCTION::POWER_LAW_ALPHA:
                prob = pow(lum, -ALPHA_HIST) * exp(-lum / L_MAX_HIST) / (improvedGamma(1 - ALPHA_HIST, L_MIN_HIST / L_MAX_HIST) * pow(L_MAX_HIST, 1 - ALPHA_HIST));
                if (lum < L_MIN_HIST) {
                    prob = 0;
                }
                break;
            case LUMINOSITY_FUNCTION::LOG_NORMAL:
                prob = log10(exp(1))/(SIGMA_HIST * sqrt(2 * pi) * lum) * exp(-pow(log10(lum) - log10(L0_HIST), 2) / (2 * SIGMA_HIST * SIGMA_HIST));
                break;
            case LUMINOSITY_FUNCTION::PLOEG:
                prob = log10(exp(1)) / (SIGMA_PLOEG_HIST * sqrt(2 * pi) * lum) * exp(-pow(log10(lum) - log10(L0_PLOEG_HIST), 2) / (2 * SIGMA_PLOEG_HIST * SIGMA_PLOEG_HIST));
                break;

            case LUMINOSITY_FUNCTION::NPTF:
                if (lum < NPTF_L_BREAK_HIST) {
                    prob = pow(lum / NPTF_L_BREAK_HIST, -NPTF_N_1_HIST);
                }
                else {
                    prob = pow(lum / NPTF_L_BREAK_HIST, -NPTF_N_2_HIST);
                }
                prob *= (1 - NPTF_N_1_HIST) * (1 - NPTF_N_2_HIST) / (NPTF_L_BREAK_HIST * (NPTF_N_1_HIST - NPTF_N_2_HIST));
                break;
            }
            integral[i] += nfwSquaredValue * radialDistance * radialDistance * deltaRadialDistance * prob * lumRange * cosLat * CM_PER_KPC_SQUARED; // Treating nfw here as a number density
        }
        radialDistance += deltaRadialDistance;
    }

    return integral;
}

void generateValueSkyMap(VALUE value, DoubleVector* skyMap, double arg1, double arg2) {
    // returns an array of the value across the entire sky.
    // The values are scaled relative to each other in the image, but do not produce the correct total luminosity across the entire GCE
    *skyMap = DoubleVector(thresholds.skyShape.first, std::vector<double>(thresholds.skyShape.second, 0));
    for (int x = 0; x < skyMap->size(); x++) {
        //std::cout << x << '/' << skyMap->size() << std::endl;
        for (int y = 0; y < (*skyMap)[x].size(); y++) {
            CoordPair latLon = thresholds.indexToLatLon({ x, y });
            if (myAbs(latLon.first) > 2 * pi / 180) {// Cut out 2 degrees around the equator on each side.
                //std::cout << "Evaluating latLon " << latLon.first << " " << latLon.second << std::endl;
                double val = getValueAtLatLon(latLon, thresholds.getFluxThreshold(latLon), value, arg1, arg2);
                (*skyMap)[x][y] = val;
            }
            else {
                (*skyMap)[x][y] = 0;
            }
        }
    }
}

double getValueAtConfig(VALUE value, double arg1, double arg2) {
    // Return the(unscaled) value at a certain arg1 or arg2
    DoubleVector skyMap;
    generateValueSkyMap(value, &skyMap, arg1, arg2);
    double sum = 0;
    for (int x = 0; x < skyMap.size(); x++) {
        for (int y = 0; y < skyMap[x].size(); y++) {
            sum += skyMap[x][y];// The values in the sky map have already been multiplied by the position distro.
        }
    }
    return sum;
}

void generatePowerLawPlotMap(VALUE value, DoubleVector* plotMap) {
    //Return the(unscaled) plot map
    CoordPair powerSteps = { pow(L_MIN_RANGE[1] / L_MIN_RANGE[0], 1.0 / PLOT_SIZE),
                pow(L_MAX_RANGE[1] / L_MAX_RANGE[0], 1.0 / PLOT_SIZE) };

    *plotMap = DoubleVector(PLOT_SIZE, std::vector<double>(PLOT_SIZE, 0));

#ifndef LOTS_OF_THREADS
    for (int i = 0; i < PLOT_SIZE; i++) {
        double lMin = L_MIN_RANGE[0] * pow(powerSteps.first, i);
        //std::cout << lMin << std::endl;
        for (int j = 0; j < PLOT_SIZE; j++) {
            double lMax = L_MAX_RANGE[0] * pow(powerSteps.second, j);
            (*plotMap)[i][j] = getValueAtConfig(value, lMin, lMax);
        }
    }

#else 
    for (int i = 0; i < PLOT_SIZE; i++) {
        double lMin = L_MIN_RANGE[0] * pow(powerSteps.first, i);
        std::vector<std::future<double>*> futures(PLOT_SIZE, nullptr);
        for (int j = 0; j < PLOT_SIZE; j++) {
            double lMax = L_MAX_RANGE[0] * pow(powerSteps.second, j);
            //futures[j] = new std::future<double>(std::move(std::async(std::launch::async, &getValueAtConfig, value, lMin, lMax)));
        }
        std::cout << PLOT_SIZE << " threads made for value " << (int)value << "." << std::endl;
        for (int j = 0; j < PLOT_SIZE; j++) {
            (*plotMap)[i][j] = 1; //futures[j]->get();
            //delete futures[j];
        }
    }
    std::cout << "Value " << (int)value << " completed." << std::endl;
#endif
}

void generatePowerLawAlphaPlotMap(VALUE value, DoubleVector* plotMap) {
    //Return the(unscaled) plot map
    CoordPair powerSteps = { (ALPHA_RANGE[1] - ALPHA_RANGE[0]) / PLOT_SIZE,
                pow(L_MAX_RANGE[1] / L_MAX_RANGE[0], 1.0 / PLOT_SIZE) };

    *plotMap = DoubleVector(PLOT_SIZE, std::vector<double>(PLOT_SIZE, 0));

#ifndef LOTS_OF_THREADS
    for (int i = 0; i < PLOT_SIZE; i++) {
        double alpha = i * powerSteps.first + ALPHA_RANGE[0];
        //std::cout << alpha << std::endl;
        for (int j = 0; j < PLOT_SIZE; j++) {
            double lMax = L_MAX_RANGE[0] * pow(powerSteps.second, j);
            (*plotMap)[i][j] = getValueAtConfig(value, alpha, lMax);
        }
    }

#else 
    std::cout << "Unimplemented" << std::endk;
    throw "Unimplemented";
#endif
}

void generateLogNormalPlotMap(VALUE value, DoubleVector* plotMap) {
    //Return the(unscaled) plot map
    const CoordPair powerSteps = { pow(L0_RANGE[1] / L0_RANGE[0], 1.0 / PLOT_SIZE),
                    (SIGMA_RANGE[1] - SIGMA_RANGE[0]) / PLOT_SIZE };
    *plotMap = DoubleVector(PLOT_SIZE, std::vector<double>(PLOT_SIZE, 0));

#ifndef LOTS_OF_THREADS
    for (int i = 0; i < PLOT_SIZE; i++) {
        double l0 = L0_RANGE[0] * pow(powerSteps.first, i);
        //std::cout << i << " / " << plotShape.first << std::endl;
        for (int j = 0; j < PLOT_SIZE; j++) {
            double sigma = SIGMA_RANGE[0] + j * powerSteps.second;
            (*plotMap)[i][j] = getValueAtConfig(value, l0, sigma);
        }
    }

#else 
    for (int i = 0; i < PLOT_SIZE; i++) {
        double l0 = L0_RANGE[0] * pow(powerSteps.first, i);
        std::vector<std::future<double>*> futures(PLOT_SIZE, nullptr);
        for (int j = 0; j < PLOT_SIZE; j++) {
            double sigma = SIGMA_RANGE[0] + powerSteps.second * j;
            //futures[j] = new std::future<double>(std::move(std::async(std::launch::async, &getValueAtConfig, value, l0, sigma)));
        }
        std::cout << PLOT_SIZE << " threads made for value " << (int)value << "at l0 " << l0 <<"." << std::endl;
        for (int j = 0; j < PLOT_SIZE; j++) {
            (*plotMap)[i][j] = 1;//futures[j]->get();
            //delete futures[j];
        }
    }
    std::cout << "Value " << (int)value << " completed." << std::endl;
#endif
}

// ========================= Generate data =========================
int powerLaw() {
    std::cout << "Using a power law luminosity function" << std::endl << std::endl;
    double paperScale = FLUX_EXCESS / getValueAtConfig(VALUE::TOTAL_FLUX, fermilabPaperPoint[1], fermilabPaperPoint[0]);
    std::string writeText = "Paper values: (lMin = " + sciNot(fermilabPaperPoint[1]) + " ergs/s, lMax = " + sciNot(fermilabPaperPoint[0]) + " ergs/s)\n" +
        "\tTotal number of pulsars: " + sciNot(getValueAtConfig(VALUE::TOTAL_NUM, fermilabPaperPoint[1], fermilabPaperPoint[0]) * paperScale) +
        "\n\tNumber of visible pulsars: " + sciNot(getValueAtConfig(VALUE::SEEN_NUM, fermilabPaperPoint[1], fermilabPaperPoint[0]) * paperScale) +
        "\n\tFraction of seen luminosity: " + sciNot(getValueAtConfig(VALUE::SEEN_FLUX, fermilabPaperPoint[1], fermilabPaperPoint[0]) * paperScale / FLUX_EXCESS) + "\n";
    std::cout << paperScale << std::endl;
    std::cout << writeText << std::endl;
    std::ofstream recordFile;
    recordFile.open(ROOT "luminosity-models-position/data-" + std::to_string((int)SENSITIVITY_DIVISOR) + "x/power-law/record.txt");
    recordFile << writeText;

    // Generate unscaled data
    DoubleVector totalNum, numSeen, fluxSeen, totalFlux;

#ifndef LOTS_OF_THREADS
    std::thread totalNumThread(generatePowerLawPlotMap, VALUE::TOTAL_NUM, &totalNum);
    std::thread numSeenThread(generatePowerLawPlotMap, VALUE::SEEN_NUM, &numSeen);
    std::thread lumSeenThread(generatePowerLawPlotMap, VALUE::SEEN_FLUX, &fluxSeen);
    std::thread totalLumThread(generatePowerLawPlotMap, VALUE::TOTAL_FLUX, &totalFlux);
    totalNumThread.join();
    numSeenThread.join();
    lumSeenThread.join();
    totalLumThread.join();
#else
    // Do two at a time, because 50 * 2 = 100, and Erebus uses 80 threads at a time./
    // Of course, the number of threads is weakly less than 100.
    std::thread totalNumThread(generatePowerLawPlotMap, VALUE::TOTAL_NUM, &totalNum);
    std::thread numSeenThread(generatePowerLawPlotMap, VALUE::SEEN_NUM, &numSeen);
    totalNumThread.join();
    numSeenThread.join();
    std::thread lumSeenThread(generatePowerLawPlotMap, VALUE::SEEN_FLUX, &fluxSeen);
    std::thread totalLumThread(generatePowerLawPlotMap, VALUE::TOTAL_FLUX, &totalFlux);
    lumSeenThread.join();
    totalLumThread.join();
#endif

    // Write data to an output file
    std::ofstream totalNumFile, numSeenFile, lumSeenFile;
    totalNumFile.open(ROOT "luminosity-models-position/data-" + std::to_string((int)SENSITIVITY_DIVISOR) + "x/power-law/total-num.txt");
    numSeenFile.open(ROOT "luminosity-models-position/data-" + std::to_string((int)SENSITIVITY_DIVISOR) + "x/power-law/num-seen.txt");
    lumSeenFile.open(ROOT "luminosity-models-position/data-" + std::to_string((int)SENSITIVITY_DIVISOR) + "x/power-law/lum-seen.txt");
    for (int x = 0; x < totalNum.size(); x++) {
        for (int y = 0; y < totalNum[x].size(); y++) {
            double scale = FLUX_EXCESS / totalFlux[x][y];
            totalNumFile << totalNum[x][y] * scale << (y == totalNum[x].size() - 1 ? "" : ", ");
            numSeenFile << numSeen[x][y] * scale << (y == totalNum[x].size() - 1 ? "" : ", ");
            lumSeenFile << fluxSeen[x][y] * scale << (y == totalNum[x].size() - 1 ? "" : ", ");
        }
        totalNumFile << std::endl;
        numSeenFile << std::endl;
        lumSeenFile << std::endl;
    }

    std::cout << "Done." << std::endl;
    system("python \"" ROOT "luminosity-models-position/graph-power-law.py\"");
    std::cin.get();
    return 1;
}

int powerLawAlpha() {
    std::cout << "Using a power law alpha luminosity function" << std::endl << std::endl;
    double paperScale = FLUX_EXCESS / getValueAtConfig(VALUE::TOTAL_FLUX, FERMILAB_ALPHA, fermilabPaperPoint[0]);
    std::string writeText = "Paper values: (alpha = " + std::to_string((int)FERMILAB_ALPHA) + ", lMax = " + sciNot(fermilabPaperPoint[0]) + " ergs/s)\n" +
        "\tTotal number of pulsars: " + sciNot(getValueAtConfig(VALUE::TOTAL_NUM, FERMILAB_ALPHA, fermilabPaperPoint[0]) * paperScale) +
        "\n\tNumber of visible pulsars: " + sciNot(getValueAtConfig(VALUE::SEEN_NUM, FERMILAB_ALPHA, fermilabPaperPoint[0]) * paperScale) +
        "\n\tFraction of seen luminosity: " + sciNot(getValueAtConfig(VALUE::SEEN_FLUX, FERMILAB_ALPHA, fermilabPaperPoint[0]) * paperScale / FLUX_EXCESS) + "\n";

    std::cout << writeText << std::endl;
    std::ofstream recordFile;
    recordFile.open(ROOT "luminosity-models-position/data-" + std::to_string((int)SENSITIVITY_DIVISOR) + "x/power-law-alpha/record.txt");
    recordFile << writeText;

    // Generate unscaled data
    DoubleVector totalNum, numSeen, fluxSeen, totalFlux;

#ifndef LOTS_OF_THREADS
    std::thread totalNumThread(generatePowerLawAlphaPlotMap, VALUE::TOTAL_NUM, &totalNum);
    std::thread numSeenThread(generatePowerLawAlphaPlotMap, VALUE::SEEN_NUM, &numSeen);
    std::thread lumSeenThread(generatePowerLawAlphaPlotMap, VALUE::SEEN_FLUX, &fluxSeen);
    std::thread totalLumThread(generatePowerLawAlphaPlotMap, VALUE::TOTAL_FLUX, &totalFlux);
    totalNumThread.join();
    numSeenThread.join();
    lumSeenThread.join();
    totalLumThread.join();
#else
    // Do two at a time, because 50 * 2 = 100, and Erebus uses 80 threads at a time./
    // Of course, the number of threads is weakly less than 100.
    std::thread totalNumThread(generatePowerLawPlotMap, VALUE::TOTAL_NUM, &totalNum);
    std::thread numSeenThread(generatePowerLawPlotMap, VALUE::SEEN_NUM, &numSeen);
    totalNumThread.join();
    numSeenThread.join();
    std::thread lumSeenThread(generatePowerLawPlotMap, VALUE::SEEN_FLUX, &fluxSeen);
    std::thread totalLumThread(generatePowerLawPlotMap, VALUE::TOTAL_FLUX, &totalFlux);
    lumSeenThread.join();
    totalLumThread.join();
#endif

    // Write data to an output file
    std::ofstream totalNumFile, numSeenFile, lumSeenFile;
    totalNumFile.open(ROOT "luminosity-models-position/data-" + std::to_string((int)SENSITIVITY_DIVISOR) + "x/power-law-alpha/total-num.txt");
    numSeenFile.open(ROOT "luminosity-models-position/data-" + std::to_string((int)SENSITIVITY_DIVISOR) + "x/power-law-alpha/num-seen.txt");
    lumSeenFile.open(ROOT "luminosity-models-position/data-" + std::to_string((int)SENSITIVITY_DIVISOR) + "x/power-law-alpha/lum-seen.txt");
    for (int x = 0; x < totalNum.size(); x++) {
        for (int y = 0; y < totalNum[x].size(); y++) {
            double scale = FLUX_EXCESS / totalFlux[x][y];
            totalNumFile << totalNum[x][y] * scale << (y == totalNum[x].size() - 1 ? "" : ", ");
            numSeenFile << numSeen[x][y] * scale << (y == totalNum[x].size() - 1 ? "" : ", ");
            lumSeenFile << fluxSeen[x][y] * scale << (y == totalNum[x].size() - 1 ? "" : ", ");
        }
        totalNumFile << std::endl;
        numSeenFile << std::endl;
        lumSeenFile << std::endl;
    }

    std::cout << "Done." << std::endl;
    system("python \"" ROOT "luminosity-models-position/graph-power-law-alpha.py\"");
    std::cin.get();
    return 1;
}

int logNormal() {
    std::cout << "Using a log normal luminosity function" << std::endl << std::endl;
    double paperScale = FLUX_EXCESS / getValueAtConfig(VALUE::TOTAL_FLUX, logNormalPaperPoint[0], logNormalPaperPoint[1]);;
    std::string paperText = "Paper values: (l0 = " + sciNot(logNormalPaperPoint[0]) + " ergs/s, sigma = " + sciNot(logNormalPaperPoint[1]) + ")\n"
        + "\tTotal number of pulsars: " + sciNot(getValueAtConfig(VALUE::TOTAL_NUM, logNormalPaperPoint[0], logNormalPaperPoint[1]) * paperScale)
        + "\n\tNumber of visible pulsars: " + sciNot(getValueAtConfig(VALUE::SEEN_NUM, logNormalPaperPoint[0], logNormalPaperPoint[1]) * paperScale)
        + "\n\tFraction of seen luminosity: " + sciNot(getValueAtConfig(VALUE::SEEN_FLUX, logNormalPaperPoint[0], logNormalPaperPoint[1]) * paperScale / FLUX_EXCESS) + "\n";

    double ploegScale = FLUX_EXCESS / getValueAtConfig(VALUE::TOTAL_FLUX, logNormalPloegPoint[0], logNormalPloegPoint[1]);
    std::string ploegText = "Ploeg values: (l0 = " + sciNot(logNormalPloegPoint[0]) + " ergs/s, sigma = " + sciNot(logNormalPloegPoint[1]) + ")\n"
        + "\tTotal number of pulsars: " + sciNot(getValueAtConfig(VALUE::TOTAL_NUM, logNormalPloegPoint[0], logNormalPloegPoint[1]) * ploegScale)
        + "\n\tNumber of visible pulsars: " + sciNot(getValueAtConfig(VALUE::SEEN_NUM, logNormalPloegPoint[0], logNormalPloegPoint[1]) * ploegScale)
        + "\n\tFraction of seen luminosity: " + sciNot(getValueAtConfig(VALUE::SEEN_FLUX, logNormalPloegPoint[0], logNormalPloegPoint[1]) * ploegScale / FLUX_EXCESS) + "\n";

    std::cout << paperText << std::endl;
    std::cout << ploegText << std::endl;
    std::ofstream recordFile;
    recordFile.open(ROOT "luminosity-models-position/data-" + std::to_string((int)SENSITIVITY_DIVISOR) + "x/log-normal/record.txt");
    recordFile << paperText << std::endl;
    recordFile << ploegText << std::endl;

    // Generate unscaled data
    DoubleVector totalNum, numSeen, fluxSeen, totalFlux;

#ifndef LOTS_OF_THREADS
    std::thread totalNumThread(generateLogNormalPlotMap, VALUE::TOTAL_NUM, &totalNum);
    std::thread numSeenThread(generateLogNormalPlotMap, VALUE::SEEN_NUM, &numSeen);
    std::thread lumSeenThread(generateLogNormalPlotMap, VALUE::SEEN_FLUX, &fluxSeen);
    std::thread totalLumThread(generateLogNormalPlotMap, VALUE::TOTAL_FLUX, &totalFlux);
    totalNumThread.join();
    numSeenThread.join();
    lumSeenThread.join();
    totalLumThread.join();
#else
    // Do two at a time, because 50 * 2 = 100, and Erebus uses 80 threads at a time./
    // Of course, the number of threads is weakly less than 98.
    std::thread totalNumThread(generateLogNormalPlotMap, VALUE::TOTAL_NUM, &totalNum);
    std::thread numSeenThread(generateLogNormalPlotMap, VALUE::SEEN_NUM, &numSeen);
    totalNumThread.join();
    numSeenThread.join();
    std::thread lumSeenThread(generateLogNormalPlotMap, VALUE::SEEN_FLUX, &fluxSeen);
    std::thread totalLumThread(generateLogNormalPlotMap, VALUE::TOTAL_FLUX, &totalFlux);
    lumSeenThread.join();
    totalLumThread.join();
#endif

    // Write data to an output file
    std::ofstream totalNumFile, numSeenFile, lumSeenFile;
    totalNumFile.open(ROOT "luminosity-models-position/data-" + std::to_string((int)SENSITIVITY_DIVISOR) + "x/log-normal/total-num.txt");
    numSeenFile.open(ROOT "luminosity-models-position/data-" + std::to_string((int)SENSITIVITY_DIVISOR) + "x/log-normal/num-seen.txt");
    lumSeenFile.open(ROOT "luminosity-models-position/data-" + std::to_string((int)SENSITIVITY_DIVISOR) + "x/log-normal/lum-seen.txt");
    for (int x = 0; x < totalNum.size(); x++) {
        for (int y = 0; y < totalNum[x].size(); y++) {
            double scale = FLUX_EXCESS / totalFlux[x][y];
            totalNumFile << totalNum[x][y] * scale << (y == totalNum[x].size() - 1 ? "" : ", ");
            numSeenFile << numSeen[x][y] * scale << (y == totalNum[x].size() - 1 ? "" : ", ");
            lumSeenFile << fluxSeen[x][y] * scale << (y == totalNum[x].size() - 1 ? "" : ", ");
        }
        totalNumFile << std::endl;
        numSeenFile << std::endl;
        lumSeenFile << std::endl;
    }

    std::cout << "Done." << std::endl;
    system("python \"" ROOT "luminosity-models-position/graph-log-normal.py\"");
    std::cin.get();
    return 1;
}

int ploeg() {
    std::cout << "Using Ploeg et al.'s luminosity function" << std::endl << std::endl;

    std::future<double> totalFlux = std::async(std::launch::async, &getValueAtConfig, VALUE::TOTAL_FLUX, 0, 0);
    std::future<double> totalNum = std::async(std::launch::async, &getValueAtConfig, VALUE::TOTAL_NUM, 0, 0);
    std::future<double> seenFlux = std::async(std::launch::async, &getValueAtConfig, VALUE::SEEN_FLUX, 0, 0);
    std::future<double> seenNum = std::async(std::launch::async, &getValueAtConfig, VALUE::SEEN_NUM, 0, 0);

    double ploegScale = FLUX_EXCESS / totalFlux.get();
    double totalNumScaled = totalNum.get() * ploegScale;
    double numVisibleScaled = seenNum.get() * ploegScale;
    double fracFluxScaled = seenFlux.get() * ploegScale / FLUX_EXCESS;

    std::string ploegText = "Ploeg luminosity function \n"
        "\tTotal number of pulsars: " + sciNot(totalNumScaled)
        + "\n\tNumber of visible pulsars: " + sciNot(numVisibleScaled)
        + "\n\tFraction of seen luminosity: " + sciNot(fracFluxScaled) + "\n";

    std::cout << ploegText << std::endl;
    std::ofstream recordFile;
    recordFile.open(ROOT "luminosity-models-position/data-" + std::to_string((int)SENSITIVITY_DIVISOR) + "x/ploeg/record.txt");
    recordFile << ploegText << std::endl;

    std::cin.get();
    return 1;
}

int nptf() {
    std::cout << "Using the NPTF luminosity function" << std::endl << std::endl;

    nptfLBreak = 8.656487610122969e+33;
    setNPTFPremul(-0.66, 18.2);

    std::future<double> totalFlux = std::async(std::launch::async, getValueAtConfig, VALUE::TOTAL_FLUX, -0.66, 18.2);
    std::future<double> totalNum = std::async(std::launch::async, getValueAtConfig, VALUE::TOTAL_NUM, -0.66, 18.2);
    std::future<double> seenFlux = std::async(std::launch::async, getValueAtConfig, VALUE::SEEN_FLUX, -0.66, 18.2);
    std::future<double> seenNum = std::async(std::launch::async, getValueAtConfig, VALUE::SEEN_NUM, -0.66, 18.2);    

    double nptfScale = FLUX_EXCESS / totalFlux.get();
    double totalNumScaled = totalNum.get() * nptfScale;
    double numVisibleScaled = seenNum.get() * nptfScale;
    double fracFluxScaled = seenFlux.get() * nptfScale / FLUX_EXCESS;

    std::string nptfText = "NPTF NFW PS luminosity function (nBelow = -0.66, nAbove = 18.2, fluxBreak = " + sciNot(nptfLBreak) + " ergs/s/cm^2 \n"
        "\tTotal number of pulsars: " + sciNot(totalNumScaled)
        + "\n\tNumber of visible pulsars: " + sciNot(numVisibleScaled)
        + "\n\tFraction of seen luminosity: " + sciNot(fracFluxScaled) + "\n\n";


    nptfLBreak = 3.344552031183874e+35;
    nptfLMin = 1e29;
    setNPTFPremul(1.40, 17.5);

    totalFlux = std::async(std::launch::async, &getValueAtConfig, VALUE::TOTAL_FLUX, 1.40, 17.5);
    totalNum = std::async(std::launch::async, &getValueAtConfig, VALUE::TOTAL_NUM, 1.40, 17.5);
    seenFlux = std::async(std::launch::async, &getValueAtConfig, VALUE::SEEN_FLUX, 1.40, 17.5);
    seenNum = std::async(std::launch::async, &getValueAtConfig, VALUE::SEEN_NUM, 1.40, 17.5);

    nptfScale = FLUX_EXCESS / totalFlux.get();
    totalNumScaled = totalNum.get() * nptfScale;
    numVisibleScaled = seenNum.get() * nptfScale;
    fracFluxScaled = seenFlux.get() * nptfScale / FLUX_EXCESS;

    nptfText += "NPTF Disk PS luminosity function (nBelow = 1.40, nAbove = 17.5, fluxBreak = " + sciNot(nptfLBreak) + " ergs/s/cm^2 \n"
        "\tTotal number of pulsars: " + sciNot(totalNumScaled)
        + "\n\tNumber of visible pulsars: " + sciNot(numVisibleScaled)
        + "\n\tFraction of seen luminosity: " + sciNot(fracFluxScaled) + "\n";

    std::cout << nptfText << std::endl; std::ofstream recordFile;
    recordFile.open(ROOT "luminosity-models-position/data-" + std::to_string((int)SENSITIVITY_DIVISOR) + "x/nptf/record.txt");
    recordFile << nptfText << std::endl;

    std::cin.get();
    return 1;
}


int hist() {
    nptfLBreak = NPTF_L_BREAK_HIST;
    nptfLMin = 1e29;
    setNPTFPremul(NPTF_N_1_HIST, NPTF_N_2_HIST);


    switch (luminosityFunction) {
    case LUMINOSITY_FUNCTION::POWER_LAW:
    case LUMINOSITY_FUNCTION::POWER_LAW_ALPHA:
        std::cout << "Using a power law luminosity function" << std::endl << std::endl;
        break;
    case LUMINOSITY_FUNCTION::LOG_NORMAL:
        std::cout << "Using a log normal luminosity function" << std::endl << std::endl;
        break;
    case LUMINOSITY_FUNCTION::PLOEG:
        std::cout << "Using a log normal luminosity function (Ploeg)" << std::endl << std::endl;
        break;
    case LUMINOSITY_FUNCTION::NPTF:
        std::cout << "Using an NPTF luminosity function" << std::endl << std::endl;
        break;
    }


    for (double i = 0; i < FLUX_HIST_BINS; i++) {
        double low = pow(10, (log10(FLUX_HIST_HIGH) - log10(FLUX_HIST_LOW)) * (i / FLUX_HIST_BINS) + log10(FLUX_HIST_LOW));
        double high = pow(10, (log10(FLUX_HIST_HIGH) - log10(FLUX_HIST_LOW)) * ((i + 1) / FLUX_HIST_BINS) + log10(FLUX_HIST_LOW));
        fluxValues[i] = (low + high) / 2;
        fluxWidths[i] = high - low;
    }

    std::array<double, FLUX_HIST_BINS> probs;
    probs.fill(0);

    for (int x = 0; x < thresholds.skyShape.first; x++) {
        //std::cout << x << '/' << skyMap->size() << std::endl;
        for (int y = 0; y < thresholds.skyShape.second; y++) {
            CoordPair latLon = thresholds.indexToLatLon({ x, y });
            if (myAbs(latLon.first) > 2 * pi / 180) {// Cut out 2 degrees around the equator on each side.
                std::array<double, FLUX_HIST_BINS> probs_here = getHistAtLatLon(latLon, thresholds.getFluxThreshold(latLon));
                addArrays(probs, probs_here);
            }
        }
    }


    double scale;
    switch (luminosityFunction) {
    case LUMINOSITY_FUNCTION::POWER_LAW:
    case LUMINOSITY_FUNCTION::POWER_LAW_ALPHA:
        scale = FLUX_EXCESS / getValueAtConfig(VALUE::TOTAL_FLUX, L_MIN_HIST, L_MAX_HIST);
        break;
    case LUMINOSITY_FUNCTION::LOG_NORMAL:
        scale = FLUX_EXCESS / getValueAtConfig(VALUE::TOTAL_FLUX, L0_HIST, SIGMA_HIST);
        break;
    case LUMINOSITY_FUNCTION::PLOEG:
        scale = FLUX_EXCESS / getValueAtConfig(VALUE::TOTAL_FLUX, L0_PLOEG_HIST, SIGMA_PLOEG_HIST);
        break;
    case LUMINOSITY_FUNCTION::NPTF:
        scale = FLUX_EXCESS / getValueAtConfig(VALUE::TOTAL_FLUX, NPTF_N_1_HIST, NPTF_N_2_HIST);
        break;
    };

    std::cout << scale << std::endl;


    for (double& prob : probs) {
        prob *= scale;
    }

    // Write data to an output file
    std::ofstream histFile;
    switch (luminosityFunction) {
    case LUMINOSITY_FUNCTION::POWER_LAW:
    case LUMINOSITY_FUNCTION::POWER_LAW_ALPHA:
        histFile.open(ROOT "luminosity-models-position/data-" + std::to_string((int)SENSITIVITY_DIVISOR) + "x/power-law/flux-hist.txt");
        histFile << "alpha = " << ALPHA_HIST << ", L_min = " << L_MIN_HIST << ", L_max = " << L_MAX_HIST << std::endl;
        break;
    case LUMINOSITY_FUNCTION::LOG_NORMAL:
        histFile.open(ROOT "luminosity-models-position/data-" + std::to_string((int)SENSITIVITY_DIVISOR) + "x/log-normal/flux-hist.txt");
        histFile << "L0 = " << L0_HIST << ", sigma = " << SIGMA_HIST << std::endl;
        break;
    case LUMINOSITY_FUNCTION::PLOEG:
        histFile.open(ROOT "luminosity-models-position/data-" + std::to_string((int)SENSITIVITY_DIVISOR) + "x/ploeg/flux-hist.txt");
        histFile << "L0 = " << L0_PLOEG_HIST << ", sigma = " << SIGMA_PLOEG_HIST << std::endl;
        break;
    case LUMINOSITY_FUNCTION::NPTF:
        histFile.open(ROOT "luminosity-models-position/data-" + std::to_string((int)SENSITIVITY_DIVISOR) + "x/nptf/flux-hist.txt");
        histFile << "L_break = " << NPTF_L_BREAK_HIST << ", n1 = " << NPTF_N_1_HIST << ", n2 = " << NPTF_N_2_HIST << std::endl;
        break;
    }
    for (double value : probs) {
        histFile << value << ',';
    }
    histFile << std::endl;

    for (double value : fluxValues) {
        histFile << value << ',';
    }
    histFile << std::endl;

    std::cout << "Done." << std::endl;
    std::cin.get();
    return 1;
}


int main(int argc, char** argv) {
    std::cout << "Sensitivity " << SENSITIVITY_DIVISOR << std::endl;
    METHOD method = METHOD::HIST;
    if (argc == 3) {
        if (strcmp(argv[2], "powerlaw") == 0) {
            luminosityFunction = LUMINOSITY_FUNCTION::POWER_LAW;
        }
        else if (strcmp(argv[2], "lognormal") == 0) {
            luminosityFunction = LUMINOSITY_FUNCTION::LOG_NORMAL;
        }
        else if (strcmp(argv[2], "ploeg") == 0) {
            luminosityFunction = LUMINOSITY_FUNCTION::PLOEG;
        }
        else if (strcmp(argv[2], "nptf") == 0) {
            luminosityFunction = LUMINOSITY_FUNCTION::NPTF;
        }
        else if (strcmp(argv[2], "alpha") == 0) {
            luminosityFunction = LUMINOSITY_FUNCTION::POWER_LAW_ALPHA;
        }
        else {
            luminosityFunction = LUMINOSITY_FUNCTION::ERROR;
        }

        if (strcmp(argv[1], "count") == 0) {
            method = METHOD::COUNT;
        }
        else if (strcmp(argv[1], "hist") == 0) {
            method = METHOD::HIST;
        }
        else {
            method = METHOD::ERROR;
        }
    }

    if (method == METHOD::COUNT) {
        std::cout << "Count" << std::endl;
        switch (luminosityFunction) {
        case LUMINOSITY_FUNCTION::POWER_LAW:
            return powerLaw();
        case LUMINOSITY_FUNCTION::LOG_NORMAL:
            return logNormal();
        case LUMINOSITY_FUNCTION::PLOEG:
            return ploeg();
        case LUMINOSITY_FUNCTION::NPTF:
            return nptf();
        case LUMINOSITY_FUNCTION::POWER_LAW_ALPHA:
            return powerLawAlpha();
        case LUMINOSITY_FUNCTION::ERROR:
            std::cout << "The argument \"" + std::string(argv[2]) + "\" is not supported. Options are \"powerlaw\", \"lognormal\", \"ploeg\", \"alpha\" and \"nptf\"." << std::endl;
            return 0;
        }
    }

    else if (method == METHOD::HIST) {
        std::cout << "Hist" << std::endl;
        switch (luminosityFunction) {
        case LUMINOSITY_FUNCTION::ERROR:
            std::cout << "The argument \"" + std::string(argv[2]) + "\" is not supported. Options are \"powerlaw\", \"lognormal\", \"ploeg\", \"alpha\" and \"nptf\"." << std::endl;
            return 0;
        default:
            return hist();
        }
    }

    else {
        std::cout << "The argument \"" + std::string(argv[1]) + "\" is not supported. Options are \"count\" and \"hist\"." << std::endl;
        return 0;
    }
    return 1;
}