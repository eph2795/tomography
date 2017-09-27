#ifndef THRESHOLD_H
#define THRESHOLD_H

#include <vector>

#include <stxxl/vector>

#include "predefined.h"

class Threshold {
public:
    Threshold(int nPhases_, const std::vector<uchar> &Low_, const std::vector<uchar> &High_) :
            nPhases(nPhases_), specifiedLowThresholds(Low_), specifiedHighThresholds(High_) { }

    Threshold() {}

    virtual std::vector<uchar>LowThresholds() const {
        return specifiedLowThresholds;
    }

    virtual std::vector<uchar> HighThresholds() const {
        return specifiedHighThresholds;
    }

    virtual int Low() const {
        return specifiedHighThresholds[0];
    }

    virtual int High() const {
        return specifiedLowThresholds[specifiedHighThresholds.size() - 1];
    }

    virtual int PhasesCount() const {
        return nPhases;
    }

private:
    int nPhases;
    std::vector<uchar> specifiedLowThresholds, specifiedHighThresholds;
};

#endif
