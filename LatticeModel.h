#ifndef LATTICEMODEL_H
#define LATTICEMODEL_H

#include <vector>

#include <stxxl/vector>

#include "threshold.h"
#include "predefined.h"

class LatticeModel {
public:
    LatticeModel(int mWidth_, int mHeight_, int mDepth_, const Threshold &threshold_) :
            mWidth(mWidth_), mHeight(mHeight_), mDepth(mDepth_),
            means(threshold_.PhasesCount()), vars(threshold_.PhasesCount()),
            threshold(threshold_) {
        mLabels = threshold_.PhasesCount();
        prev.resize(mWidth * mHeight);
        cur.resize(mWidth * mHeight);
        next.resize(mWidth * mHeight);
        prevPhases.resize(mWidth * mHeight);
        curPhases.resize(mWidth * mHeight);
        nextPhases.resize(mWidth * mHeight);
    }

    LatticeModel(const LatticeModel&) {}

    virtual ~LatticeModel() {}

    void ConditionalImage(const stxxl::vector<uchar> &I, stxxl::vector<uchar> &C) const {
        C.resize(mWidth * mHeight * mDepth);

        std::vector<uchar> lows = threshold.LowThresholds();
        std::vector<uchar> highs = threshold.HighThresholds();

        for (std::pair<LatticeConstIterator_stxxl, LatticeIterator_stxxl> iters(I.begin(), C.begin());
             iters.first != I.end(); iters.first++, iters.second++) {

            *(iters.second) = mLabels;
            uchar v = *(iters.first);
            double min_delta = INF;

            for (size_t i = 0; i < mLabels; i++) {
                if ((v >= lows[i]) && (v < highs[i])) {
                    *(iters.second) = i;
                    break;
                }

                double cur_delta = abs(*(iters.first) - means[i]);
                if (vars[i] != 0) {
                    cur_delta /= vars[i];
                }
                if (cur_delta < min_delta) {
                    min_delta = cur_delta;
                    *(iters.second) = i;
                }

            }
		}
	}

    inline int SelectPhase(uchar c) const {
        size_t nPhase = 0;
        double min_delta = INF;
        for (size_t i = 0; i < mLabels; i++) {
            double cur_delta = abs(c - means[i]);
            if (vars[i] != 0) {
                cur_delta /= vars[i];
            }
            if (cur_delta < min_delta) {
                min_delta = cur_delta;
				nPhase = i;
			}
		}
		return nPhase;
	}


    void StatsNL(const stxxl::vector<uchar> &img)	{
        std::vector<uchar> Lows = threshold.LowThresholds();
        std::vector<uchar> Highs = threshold.HighThresholds();

        std::vector<long long> S(mLabels, 0);
        std::fill(means.begin(), means.end(), 0);
        std::fill(vars.begin(), vars.end(), 0);

        for (LatticeConstIterator_stxxl it = img.begin(); it != img.end(); it++) {
            for (size_t i = 0; i < mLabels; i++) {
                if ((*it >= Lows[i]) && (*it < Highs[i])) {
                    means[i] += *it;
					S[i]++;
				};
			}
		}

        std::transform(means.begin(), means.end(), S.begin(), means.begin(),
                       [](double &mean, long long count) { return count > 0 ? mean / count : 0; });

        for (LatticeConstIterator_stxxl it = img.begin(); it != img.end(); it++) {
            uchar v = *it;
            for (size_t i = 0; i < mLabels; i++) {
                if ((v >= Lows[i]) && (v < Highs[i])) {
                    double d = (v - means[i]);
                    vars[i] += d * d;
                }
			}
		}

        std::transform(vars.begin(), vars.end(), S.begin(), vars.begin(),
                       [](double &var, long long count) { return count > 0 ? sqrt(var / count) : 0; });
	}

    void GetStatsNL(std::vector<double> &means_, std::vector<double> &vars_) const {
        means_.resize(mLabels);
        vars_.resize(mLabels);
        means_.assign(means.begin(), means.end());
        vars_.assign(means.begin(), means.end());
    }

protected:

    void setCur(size_t d_,
                const stxxl::vector<uchar> &grayscale,
                stxxl::vector<uchar> &phases) {
        d = d_;
        if (d != 0) {
            std::copy(grayscale.begin() + (d_ - 1) * mWidth * mHeight,
                      grayscale.begin() + d_ * mWidth * mHeight,
                      prev.begin());
            std::copy(phases.begin() + (d_ - 1) * mWidth * mHeight,
                      phases.begin() + d_ * mWidth * mHeight,
                      prevPhases.begin());
        }
        std::copy(grayscale.begin() + d_ * mWidth * mHeight,
                  grayscale.begin() + (d_ + 1) * mWidth * mHeight,
                  cur.begin());
        std::copy(phases.begin() + d_ * mWidth * mHeight,
                  phases.begin() + (d_ + 1) * mWidth * mHeight,
                  curPhases.begin());
        if (d != mDepth - 1) {
            std::copy(grayscale.begin() + (d_ + 1) * mWidth * mHeight,
                      grayscale.begin() + (d_ + 2) * mWidth * mHeight,
                      next.begin());
            std::copy(phases.begin() + (d_ + 1) * mWidth * mHeight,
                      phases.begin() + (d_ + 2) * mWidth * mHeight,
                      nextPhases.begin());
        }
    }

    void inc(const stxxl::vector<uchar> &grayscale,
             stxxl::vector<uchar> &phases) {
        std::copy(curPhases.begin(),
                  curPhases.end(),
                  phases.begin() + mWidth * mHeight * d);
        prev = cur;
        prevPhases = curPhases;
        cur = next;
        curPhases = nextPhases;
        d += 1;
        if (d + 1 < mDepth) {
            std::copy(grayscale.begin() + (d + 1) * mWidth * mHeight,
                      grayscale.begin() + (d + 2) * mWidth * mHeight,
                      next.begin());
            std::copy(phases.begin() + (d + 1) * mWidth * mHeight,
                      phases.begin() + (d + 2) * mWidth * mHeight,
                      nextPhases.begin());
        }
    }

    typedef std::vector<uchar>::iterator LatticeIterator;
    typedef std::vector<uchar>::const_iterator LatticeConstIterator;
    typedef stxxl::vector<uchar>::iterator LatticeIterator_stxxl;
    typedef stxxl::vector<uchar>::const_iterator LatticeConstIterator_stxxl;


    inline LatticeIterator RightNB(LatticeIterator it) {
        return it + 1;
    }

    inline LatticeIterator LeftNB(LatticeIterator it) {
        return it - 1;
    }

    inline LatticeIterator TopNB(LatticeIterator it) {
        return it - mWidth;
    }
    inline LatticeIterator BotNB(LatticeIterator it) {
        return it + mWidth;
    }

    inline LatticeIterator BhdNB(LatticeIterator it) {
        size_t shift = it - curPhases.begin();
        return nextPhases.begin() + shift;
    }

    inline LatticeIterator FrntNB(LatticeIterator it) {
        size_t shift = it - curPhases.begin();
        return prevPhases.begin() + shift;
    }

    inline LatticeConstIterator RightNB(LatticeConstIterator it) {
        return it + 1;
    }

    inline LatticeConstIterator LeftNB(LatticeConstIterator it) {
        return it - 1;
    }

    inline LatticeConstIterator TopNB(LatticeConstIterator it) {
        return it - mWidth;

    }
    inline LatticeConstIterator BotNB(LatticeConstIterator it) {
        return it + mWidth;
    }

    inline LatticeConstIterator BhdNB(LatticeConstIterator it) {
        size_t shift = it - cur.cbegin();
        return next.begin() + shift;
    }

    inline LatticeConstIterator FrntNB(LatticeConstIterator it) {
        size_t shift = it - cur.cbegin();
        return prev.begin() + shift;
    }

    size_t mLabels;
    size_t mWidth, mHeight, mDepth;
    std::vector<double> means;
    std::vector<double> vars;

    size_t d;
    std::vector<uchar> prev, cur, next;
    std::vector<uchar> prevPhases, curPhases, nextPhases;

    Threshold threshold;
};

#endif
