#include <iostream>
#include <string>
#include <vector>
#include <assert.h>
#include <fstream>
#include <sstream>
#include <thread>
#include <boost/math/special_functions/gamma.hpp>

typedef std::pair<double, double> CoordPair;
typedef std::pair<int, int> IndexPair;
typedef std::vector<std::vector<double>> DoubleVector;

#define POWER_STEP 1.1 // 1 is the minimum
#define pi 3.14159265358979323

#define ALPHA 1.94
#define L_EXCESS 6.37e36  // All units are in ergs per second
#define DIST_TO_CENTER 8.5 // kpc
#define CM_PER_KPC 3.086e21
double L_MIN_RANGE[2] = { 1.0e28, 1.0e34 };
double L_MAX_RANGE[2] = { 1.0e34, 1.0e36 };

// Values for the computation of the integral of the NFW profile.
#define GAMMA 1.2
#define NFW_SCALE_DIST 20 // kpc
#define INTEGRAL_HALF_WIDTH DIST_TO_CENTER// / 2 // kpc
#define INTEGRAL_STEPS 100000
#define LOAD_DISTROS true

#define SENSITIVITY_PATH L"C:/Users/goods/Dropbox (MIT)/GCE UROP/sensitivity/sensitivity.txt"
#define DISPLAY_SIZE 20.0

#define MULTIPROCESSING

double paperPoint[2] = { 1e35, 1e29 };
IndexPair plotShape = { int(log(L_MIN_RANGE[1] / L_MIN_RANGE[0]) / log(POWER_STEP)),
	int(log(L_MAX_RANGE[1] / L_MAX_RANGE[0]) / log(POWER_STEP)) };

class SensitivityMap {
public:
	SensitivityMap() {
		std::ifstream sensitivityFile;
		sensitivityFile.open(SENSITIVITY_PATH);
		std::string line;
		while (std::getline(sensitivityFile, line)) {
			std::vector<double> v;
			std::stringstream ss(line);
			for (double i; ss >> i;) {
				v.push_back(fluxToThresholdLuminosity(i));
				if (ss.peek() == ',') {
					ss.ignore();
				}
			}
			luminosities.push_back(v);
		}

		upperLeft[0] = int((1 - DISPLAY_SIZE / 90.0) * luminosities.size() / 2.0);
		upperLeft[1] = int((1 - DISPLAY_SIZE / 180.0) * luminosities[0].size() / 2.0);
		lowerRight[0] = int((1 + DISPLAY_SIZE / 90.0) * luminosities.size() / 2.0);
		lowerRight[1] = int((1 + DISPLAY_SIZE / 180.0) * luminosities[0].size() / 2.0);
		skyShape = { lowerRight[0] - upperLeft[0], lowerRight[1] - upperLeft[1] };
	}

public:
	double fluxToThresholdLuminosity(double flux) const {
		return flux * (4 * pi * pow(DIST_TO_CENTER * CM_PER_KPC, 2));
	}
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
	double getThreshold(CoordPair latLon) const {
		IndexPair xy = latLonToIndex(latLon);
		return luminosities[xy.first + upperLeft[0]][xy.second + upperLeft[1]];
	}

public:
	IndexPair skyShape;
	int upperLeft[2];
	int lowerRight[2];
	DoubleVector luminosities;
};

SensitivityMap thresholds;


// ================= These functions characterize the luminosity function ============

double improvedGamma(double s, double x) {
	assert(x >= 0);
	while (s < 0) {
		return (improvedGamma(s+1, x) - pow(x, s) * exp(-x)) / s;
	}
	return boost::math::tgamma(s, x);
}

double integrate(double start, double lMin, double lMax) {
	return pow(lMax, 1 - ALPHA) * improvedGamma(1 - ALPHA, start / lMax);
}

double lintegrate(double start, double lMin, double lMax) {
	return pow(lMax, 2 - ALPHA) * improvedGamma(2 - ALPHA, start / lMax);
}


// ================= Declare helper functions =======================

void generateNFWMaps(DoubleVector* numberDistro, DoubleVector* lumDistro) {
	*numberDistro = DoubleVector(thresholds.skyShape.first, std::vector<double>(thresholds.skyShape.second, 0));
	*lumDistro = DoubleVector(thresholds.skyShape.first, std::vector<double>(thresholds.skyShape.second, 0));
	if (LOAD_DISTROS) {
		std::ifstream numberFile;
		std::ifstream lumFile;
		numberFile.open("number-distro.txt");
		lumFile.open("lum-distro.txt");

		std::string line;
		int x = 0;
		while (std::getline(numberFile, line)) {
			std::stringstream ss(line);
			int y = 0;
			for (double d; ss >> d;) {
				(*numberDistro)[x][y] = d;
				if (ss.peek() == ',') {
					ss.ignore();
				}
				y++;
			}
			x++;
		}

		x = 0;
		while (std::getline(lumFile, line)) {
			std::stringstream ss(line);
			int y = 0;
			for (double d; ss >> d;) {
				(*lumDistro)[x][y] = d;
				if (ss.peek() == ',') {
					ss.ignore();
				}
				y++;
			}
			x++;
		}
	}
	else {
		std::string numberText, lumText;
		double deltaRadialDistance = INTEGRAL_HALF_WIDTH * 2 / INTEGRAL_STEPS;
		for (int x = 0; x < numberDistro->size(); x++) {
			std::cout << x << "/" << numberDistro->size() << std::endl;
			for (int y = 0; y < (*numberDistro)[x].size(); y++) {
				CoordPair latLon = thresholds.indexToLatLon({ x, y });
				double radialDistance = DIST_TO_CENTER - INTEGRAL_HALF_WIDTH;
				double numberIntegral = 0;
				double lumIntegral = 0;
				while (radialDistance < DIST_TO_CENTER + INTEGRAL_HALF_WIDTH) {
					double distFromCenter = sqrt(pow(radialDistance * cos(latLon.second) * cos(latLon.first) - DIST_TO_CENTER, 2) +
						pow(radialDistance * sin(latLon.second) * cos(latLon.first), 2) + pow(radialDistance * sin(latLon.first), 2));
					double volumeElement = deltaRadialDistance * radialDistance * radialDistance;
					double nfwValue = pow(pow(distFromCenter / NFW_SCALE_DIST, -GAMMA) * pow(1 + distFromCenter / NFW_SCALE_DIST, -3 + GAMMA), 2);
					numberIntegral += nfwValue * deltaRadialDistance * radialDistance * radialDistance; // Ignore the angular part since it introduces a constant cofactor.
					lumIntegral += nfwValue * deltaRadialDistance * DIST_TO_CENTER * DIST_TO_CENTER;// Ignore the angular part since it introduces a constant cofactor.
					radialDistance += deltaRadialDistance;
				}
				(*numberDistro)[x][y] = numberIntegral;
				(*lumDistro)[x][y] = lumIntegral;
				numberText += std::to_string(numberIntegral) + (y == (*numberDistro)[x].size() - 1 ? "" : ",");
				lumText += std::to_string(lumIntegral) + (y == (*numberDistro)[x].size() - 1 ? "" : ",");
			}
			numberText += (x == numberDistro->size() - 1 ? "" : "\n");
			lumText += (x == numberDistro->size() - 1 ? "" : "\n");
		}
		std::ofstream numberFile;
		std::ofstream lumFile;
		numberFile.open("number-distro.txt");
		lumFile.open("lum-distro.txt");
		numberFile << numberText;
		lumFile << lumText;
	}
}

double totalNumFunc(double threshold, double lMin, double lMax) {
	return integrate(lMin, lMin, lMax);
}

double numSeenFunc(double threshold, double lMin, double lMax) {// Returns(unscaled) number of visible pulsars
	assert(threshold > lMin);
	return integrate(threshold, lMin, lMax);
}

double lumSeenFunc(double threshold, double lMin, double lMax) {// Returns(unscaled) amount of luminosity visible
	assert(threshold > lMin);
	return lintegrate(threshold, lMin, lMax);
}

double totalLumFunc(double threshold, double lMin, double lMax) {// Returns(unscaled) total luminosity
	return lintegrate(lMin, lMin, lMax);
}


double getValueAtPos(CoordPair latLon, double (*valueFunc)(double, double, double), double lMin, double lMax) {
	// returns the(unscaled) value for a specific lon and lat
	double threshold = thresholds.getThreshold(latLon);
	return valueFunc(threshold, lMin, lMax);
}

void generateValueSkyMap(double (*valueFunc)(double, double, double), DoubleVector* skyMap, DoubleVector* positionDistro, double lMin, double lMax) {
	// returns an array of numTotal, numSeen, fracSeen across the entire sky.
	// The values are scaled relative to each other in the image, but do not produce the correct total luminosity across the entire GCE
	*skyMap = DoubleVector(thresholds.skyShape.first, std::vector<double>(thresholds.skyShape.second, 0));
	for (int x = 0; x < skyMap->size(); x++) {
		for (int y = 0; y < (*skyMap)[x].size(); y++) {
			CoordPair latLon = thresholds.indexToLatLon({ x, y });
			double val = (*positionDistro)[x][y] * getValueAtPos(latLon, valueFunc, lMin, lMax);
			(*skyMap)[x][y] = val;
		}
	}
}

double getValueAtConfig(double (*valueFunc)(double, double, double), DoubleVector* positionDistro, double lMin, double lMax) {
	// Return the(unscaled) value at a certain lMin or lMax
	DoubleVector skyMap;
	generateValueSkyMap(valueFunc, &skyMap, positionDistro, lMin, lMax);
	double sum = 0;
	for (int x = 0; x < skyMap.size(); x++) {
		for (int y = 0; y < skyMap[x].size(); y++) {
			sum += skyMap[x][y];
		}
	}
	return sum;
}

void generatePlotMap(double (*valueFunc)(double, double, double), DoubleVector* plotMap, DoubleVector* positionDistro) {
	//Return the(unscaled) plot map
	*plotMap = DoubleVector(plotShape.first, std::vector<double>(plotShape.second, 0));
	for (int i = 0; i < plotShape.first; i++) {
		double lMin = L_MIN_RANGE[0] * pow(POWER_STEP, i);
		//std::cout << i << " / " << plotShape.first << std::endl;
		for (int j = 0; j < plotShape.second; j++) {
			double lMax = L_MAX_RANGE[0] * pow(POWER_STEP, j);
			(*plotMap)[i][j] = getValueAtConfig(valueFunc, positionDistro, lMin, lMax);
		}
	}
}

// ========================= Generate data =========================
int main() {
	DoubleVector numberDistro, lumDistro;
	generateNFWMaps(&numberDistro, &lumDistro);
	std::cout << "Distros generated" << std::endl;
	double paperScale = L_EXCESS / getValueAtConfig(totalLumFunc, &lumDistro, paperPoint[1], paperPoint[0]);
	std::cout << "Paper values: (lMin = " << paperPoint[1] << " ergs/s, lMax = " << paperPoint[0] << " ergs/s)" << std::endl
		<< "\tTotal number of pulsars: " << getValueAtConfig(totalNumFunc, &numberDistro, paperPoint[1], paperPoint[0]) * paperScale << std::endl
		<< "\tNumber of visible pulsars: " << getValueAtConfig(numSeenFunc, &numberDistro, paperPoint[1], paperPoint[0]) * paperScale << std::endl
		<< "\tFraction of seen luminosity: " << getValueAtConfig(lumSeenFunc, &lumDistro, paperPoint[1], paperPoint[0]) * paperScale / L_EXCESS << std::endl
		<< "\tTotal luminosity: " << getValueAtConfig(totalLumFunc, &lumDistro, paperPoint[1], paperPoint[0]) * paperScale << " ergs/s" << std::endl;

	// Generate unscaled data
	DoubleVector totalNum, numSeen, lumSeen, totalLum;

#ifdef MULTIPROCESSING
	std::thread totalNumThread(generatePlotMap, totalNumFunc, &totalNum, &numberDistro);
	std::thread numSeenThread(generatePlotMap, numSeenFunc, &numSeen, &numberDistro);
	std::thread lumSeenThread(generatePlotMap, lumSeenFunc, &lumSeen, &lumDistro);
	std::thread totalLumThread(generatePlotMap, totalLumFunc, &totalLum, &lumDistro);
	totalNumThread.join();
	numSeenThread.join();
	lumSeenThread.join();
	totalLumThread.join();
#else
	generatePlotMap(totalNumFunc, &totalNum, &numberDistro);
	std::cout << "Total number map completed" << std::endl;
	generatePlotMap(numSeenFunc, &numSeen, &numberDistro);
	std::cout << "Seen number map completed" << std::endl;
	generatePlotMap(lumSeenFunc, &lumSeen, &lumDistro);
	std::cout << "Seen luminosity map completed" << std::endl;
	generatePlotMap(totalLumFunc, &totalLum, &lumDistro);
	std::cout << "Total luminosity map completed" << std::endl;
#endif


	// Write data to an output file
	std::ofstream totalNumFile, numSeenFile, lumSeenFile;
	totalNumFile.open("C:/Users/goods/Dropbox (MIT)/GCE UROP/luminosity-models-position/power-law/data/total-num.txt");
	numSeenFile.open("C:/Users/goods/Dropbox (MIT)/GCE UROP/luminosity-models-position/power-law/data/num-seen.txt");
	lumSeenFile.open("C:/Users/goods/Dropbox (MIT)/GCE UROP/luminosity-models-position/power-law/data/lum-seen.txt");
	for (int x = 0; x < totalNum.size(); x++) {
		for (int y = 0; y < totalNum[x].size(); y++) {
			double scale = L_EXCESS / totalLum[x][y];
			totalNumFile << totalNum[x][y] * scale << (y == totalNum[x].size() - 1 ? "" : ", ");
			numSeenFile << numSeen[x][y] * scale << (y == totalNum[x].size() - 1 ? "" : ", ");
			lumSeenFile << lumSeen[x][y] * scale << (y == totalNum[x].size() - 1 ? "" : ", ");
		}
		totalNumFile << (x == totalNum.size() - 1 ? "" : "\n");
		numSeenFile << (x == totalNum.size() - 1 ? "" : "\n");
		lumSeenFile << (x == totalNum.size() - 1 ? "" : "\n");
	}
}