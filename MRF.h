#ifndef MRF_H
#define MRF_H

#include <math.h>
#include <queue>
#include <string>
#include <random>

#include "threshold.h"
#include "LatticeModel.h"
#include "predefined.h"

struct MRFSettings {
    MRFSettings(double Beta_=0.9, double FreezingSpeed_=0.98, double TStart_=4.0, MRFMethods Method_=MRF_SA) :
            Beta(Beta_), TStart(TStart_), FreezingSpeed(FreezingSpeed_), Method(Method_) {}
    double Beta; //weight of neighbour voxels
    double TStart; //start temperature
    double FreezingSpeed; //speed of freezing in simulated annealing
    MRFMethods Method;
};


class MRF : public LatticeModel {
public:
    MRF(int mWidth_, int mHeight_, int mDepth_, const Threshold &threshold_, MRFMethods Method_=MRF_SA) :
            LatticeModel(mWidth_, mHeight_, mDepth_, threshold_), Method(Method_) {
        //mean = new double[mLabels];
        //var = new double[mLabels];
        mBeta = 0.9;
        mDeltaT = 0.9;
        mT0 = 4;
        mEnergyThresh = 0.01;
    }

    MRF(const MRFSettings &pref_, int mWidth_, int mHeight_, int mDepth_, const Threshold &threshold_) :
            LatticeModel(mWidth_, mHeight_, mDepth_, threshold_),
            mBeta(pref_.Beta),  mDeltaT(pref_.FreezingSpeed), mT0(pref_.TStart), Method(pref_.Method){
        //mean = new double[mLabels];
        //var = new double[mLabels];
        mEnergyThresh = 0.01;
    }

    ~MRF() {
		//	delete[] mean;
		//	delete[] var;
    }
		
    typedef std::vector<uchar>::iterator LatticeIterator;
    typedef std::vector<uchar>::const_iterator LatticeConstIterator;
    typedef stxxl::vector<uchar>::iterator LatticeIterator_stxxl;
    typedef stxxl::vector<uchar>::const_iterator LatticeConstIterator_stxxl;

    void SimulatedAnnealing3D(const stxxl::vector<uchar> &img, stxxl::vector<uchar> &Conditional) {
        int r;

        double summa_deltaE;
        double Temperature = mT0;
        std::uniform_int_distribution<int> distribution(1, mLabels - 1);
        std::uniform_int_distribution<int> translation_distr(0, 1);

        mIterCount = 0;
        do {
            summa_deltaE = 0.0;

            setCur(0, img, Conditional);
            typedef double(MRF::*pEnergy)(LatticeConstIterator, LatticeIterator, int);

            for (size_t k = 0; k < mDepth; k++) {
                LatticeIterator cit = curPhases.begin();// +mWidth + 1;
                LatticeConstIterator it = cur.cbegin();// +mWidth + 1;
                pEnergy pE = &MRF::Energy3D;
                if (mDepth == 1)
                    pE = &MRF::Energy;
                else if (k == 0)
                    pE = &MRF::Energy3DTop;
                else if (k == mDepth - 1)
                    pE = &MRF::Energy3DBot;

                for (size_t i = 1; i < mHeight - 1; i++) {
                    for ( size_t j = 1; j < mWidth - 1; j++) {
                        size_t shift_it = i * mWidth + j;
                        it =  cur.cbegin() + shift_it;
                        cit = curPhases.begin() + shift_it;

                        r = (*cit + distribution(generator)) % mLabels;
                        double delta_E = -((this->*pE)(it, cit, *cit) - (this->*pE)(it, cit, r));

//                        if (std::exp(delta_E / Temperature) > 0.5) {
//                            std::cout << std::exp(delta_E / Temperature) << std::endl;
//                        }

                        if (translation_distr(generator) <= -delta_E / Temperature) {
                            summa_deltaE += fabs(delta_E);
                            *cit = r;
                        }
                    }
                }
                inc(img, Conditional);
            }
            Temperature -= mDeltaT;
            ++mIterCount;
        } while ((summa_deltaE / (mWidth - 2) / (mHeight - 2) / (mDepth) > mEnergyThresh) && (Temperature > EPS_THINY));
        std:: cout << "MRF itetartions: "  << mIterCount
                   << ", accumulated deltaE: " << summa_deltaE / (mWidth - 2) / (mHeight - 2) / (mDepth)
                   << ", final Temperature: " << Temperature << std::endl;
    }

//    void Gibbs(const stxxl::vector<uchar> &img, stxxl::vector<uchar> &Conditional)	{
//        size_t s;
//        double sumE;
//        double z;

//        double *Ek = new double[mLabels];
//        double Temperature = mT0;
//        mIterCount = 0;
//        size_t turnsCount = 0;

//        std::vector<uchar> cur(mWidth * mHeight);
//        std::vector<uchar> prev(mWidth * mHeight);
//        std::vector<uchar> next(mWidth * mHeight);
//        do {
//            turnsCount = 0;
//            LatticeIterator cit = Conditional.begin();// +mWidth + 1;
//            LatticeConstIterator it = img.begin();// +mWidth + 1;
//            typedef double(MRF::*pEnergy)(LatticeConstIterator it, LatticeIterator cit, int label);

//            for (size_t k = 0; k < mDepth; k++) {
//                pEnergy pE = &MRF::Energy3D;
//                if (mDepth == 1) {
//                    pE = &MRF::Energy;
//                } else if (k == 0) {
//                    pE = &MRF::Energy3DTop;
//                } else if (k == mDepth - 1) {
//                    pE = &MRF::Energy3DBot;
//                }
//                for (size_t i = 1; i < mHeight - 1; i++) {
//                    for (size_t j = 1; j < mWidth - 1; j++) {
//                        size_t shift_it = k * mWidth * mHeight + i * mWidth + j;
//                        it = img.begin() + shift_it;
//                        cit = Conditional.begin() + shift_it;

//                        sumE = 0.0;
//                        for (s = 0; s < mLabels; s++) {
//                            Ek[s] = exp(-(this->*pE)(it, cit, s) / Temperature);
//                            sumE += Ek[s];
//                        }
//                        z = 0.0;
//                        for (s = 0; s < mLabels; ++s) {
//                            z += Ek[s] / sumE;
//                            std::bernoulli_distribution isTransition(z);
//                            if (isTransition(generator)) {
//                                turnsCount += (*cit != s);
//                                *cit = s;
//                                break;
//                            }
//                        }
//                    }
//                }
//            }
//            Temperature *= mDeltaT;
//            ++mIterCount;
//        } while (((double)turnsCount / (mWidth - 2) / (mHeight - 2) / mDepth > mEnergyThresh) && (Temperature > EPS_THINY));
//        delete[] Ek;
//    }

    void Perform(const stxxl::vector<uchar> &img, stxxl::vector<uchar> &Conditional,
                 bool doStats=true, bool doConditional=true) {
        if (doStats) {
            StatsNL(img);
        }
        if (doConditional) {
            ConditionalImage(img, Conditional);
        }

        if (Method == MRF_SA) {
            SimulatedAnnealing3D(img, Conditional);
        }

//        if (Method == MRF_GIBBS) {
//            Gibbs(img, Conditional);
//        }
    }

private:
    //	double* var;
    //	double* mean;
    double mBeta;
    double mDeltaT;
    double mT0;
    size_t mIterCount;
    double mEnergyThresh;

    MRFMethods Method;

    std::default_random_engine generator;


    double Energy(LatticeConstIterator it, LatticeIterator  cit, int label)	{
        /*double NbSum = 0;
            *BotNB(cit) == label ? NbSum++ : NbSum--;
            *TopNB(cit) == label ? NbSum++ : NbSum--;
            *RightNB(cit) == label ? NbSum++ : NbSum--;
            *LeftNB(cit) == label ? NbSum++ : NbSum--;

            double v = CellOwnEnergy(it, label) + mBeta*NbSum;*/
        return CellOwnEnergy(it, label)
            + mBeta * (
                ((*LeftNB(cit) == label)
                 + (*RightNB(cit) == label)
                 + (*TopNB(cit) == label)
                 + (*BotNB(cit) == label)) * 2 - 4
                );
    }

    double Energy3D(LatticeConstIterator it, LatticeIterator  cit, int label) {
        return CellOwnEnergy(it, label)
            + mBeta * (
                ((*LeftNB(cit) == label)
                 + (*RightNB(cit) == label)
                 + (*TopNB(cit) == label)
                 + (*BotNB(cit) == label)
                 + (*FrntNB(cit) == label)
                 + (*BhdNB(cit) == label)) * 2 - 6
                );
    }

    double Energy3DTop(LatticeConstIterator it, LatticeIterator  cit, int label) {
        return CellOwnEnergy(it, label)
            + mBeta * (
                ((*LeftNB(cit) == label)
                 + (*RightNB(cit) == label)
                 + (*TopNB(cit) == label)
                 + (*BotNB(cit) == label)
                 + (*BhdNB(cit) == label)) * 2 - 5
                );
    }

    double Energy3DBot(LatticeConstIterator it, LatticeIterator cit, int label) {
        return CellOwnEnergy(it, label)
            + mBeta * (
                ((*LeftNB(cit) == label)
                 + (*RightNB(cit) == label)
                 + (*TopNB(cit) == label)
                 + (*BotNB(cit) == label)
                 + (*FrntNB(cit) == label)) * 2 - 5
                );
    }

    double CellOwnEnergy(LatticeConstIterator it, int label) {
        return log(sqrt(2.0 * M_PI * vars[label])) + pow(*it - means[label], 2) / (2.0 * vars[label]);
    }
};


//struct Point3D {
//    Point3D(size_t x_, size_t y_, size_t z_) :x(x_), y(y_), z(z_) {}
//    size_t x, y, z;
//};


//class SRG :public LatticeModel {
//public:
//    SRG(int mWidth_, int mHeight_, int mDepth_, const Threshold &threshold_, MRFMethods Method_) :
//            LatticeModel(mWidth_, mHeight_, mDepth_, threshold_), Method(Method_) {
//        //mean = new double[mLabels];
//        //var = new double[mLabels];
//    }

//    ~SRG() {
//        //	delete[] mean;
//        //	delete[] var;
//    }

//    void Perform(const stxxl::vector<uchar> &img, stxxl::vector<uchar> &Conditional) {
//        StatsNL(img);
//        ConditionalImage(img, Conditional);

//		LatticeConstIterator it = img.begin();
	
//		std::queue<Point3D>SSL;

//		LatticeIterator cit;
//        for (size_t k = 0; k < mDepth - 1; k++) {
//            for (size_t i = 1; i < mHeight - 1; i++) {
//                for (size_t j = 1; j < mWidth - 1; j++) {
//                    size_t shift = mHeight * mWidth * k + mWidth * i + j;
//                    cit = Conditional.begin() + shift;
//                    it = img.begin() + shift;
//                    if ((*cit == mLabels) && (*BotNB(cit) != mLabels || *TopNB(cit) != mLabels ||
//                         *LeftNB(cit) != mLabels || *RightNB(cit) != mLabels ||
//                        (k < mDepth - 1 && *BhdNB(cit) != mLabels) || (k > 0 && *FrntNB(cit) != mLabels))) {
//                        SSL.push(Point3D{ j, i, k });
//                    }
//                }
//            }
//		}

//		cit = Conditional.begin();
//        std::vector<size_t> nbh(mLabels);
		
//        while (!SSL.empty()) {
//			Point3D p = SSL.front();
//			SSL.pop();
//			size_t I = p.y;
//			size_t J = p.x;
//			size_t K = p.z;

//            size_t shift = mWidth * mHeight * K + mWidth * I + J;
//            it = img.begin() + shift;
//            cit = Conditional.begin() + shift;

//            for (size_t i = 0; i < mLabels; i++) {
//                nbh[i] = (I < mHeight - 1 && *BotNB(cit) == i) + (I > 0 && *TopNB(cit) == i) +
//                         (J > 0 && *LeftNB(cit) == i) + (J < mWidth - 1 && *RightNB(cit) == i) +
//                         (K < mDepth - 1 && *BhdNB(cit) == i) + (K > 0 && *FrntNB(cit) == i);
//            }
//            size_t nb_max = 0;
//            size_t nb_max_idx = 0;

//            for (size_t i = 0; i < mLabels; i++) {
//                if (nbh[i] > nb_max) {
//					nb_max_idx = i;
//					nb_max = nbh[i];
//				}
//			}

//            if (nb_max >= 1) {
//				*cit = nb_max_idx;
//            } else {
//                *cit = SelectPhase(*cit);
//			}

//            if (I < mHeight - 1  && *BotNB(cit) == mLabels) {
//				SSL.push(Point3D{ J, I + 1, K });
//                *BotNB(cit) = mLabels + 1;
//			}
//            if (I > 0  && *TopNB(cit) == mLabels) {
//				SSL.push(Point3D{ J, I - 1, K });
//                *TopNB(cit) = mLabels + 1;
//			}

//            if (J > 0 && *LeftNB(cit) == mLabels) {
//				SSL.push(Point3D{ J - 1, I, K });
//                *LeftNB(cit) = mLabels + 1;
//			}

//            if (J < mWidth - 1 && *RightNB(cit) == mLabels) {
//				SSL.push(Point3D{ J + 1, I, K });
//                *RightNB(cit) = mLabels + 1;
//			}

//            if (K > 0 && *FrntNB(cit) == mLabels) {
//                SSL.push(Point3D{ J, I, K - 1 });
//                *FrntNB(cit) = mLabels + 1;
//			}

//            if (K < mDepth - 1 && * BhdNB(cit) == mLabels) {
//                SSL.push(Point3D{ J , I, K + 1 });
//                *BhdNB(cit) = mLabels + 1;
//			}
//		}
//	}
//private:
//	//double* var;
//	//double* mean;
//	int mIterCount;
//    MRFMethods Method;
//};


#endif
