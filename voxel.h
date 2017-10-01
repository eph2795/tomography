#ifndef VOXEL_H
#define VOXEL_H

#include <array>
#include <set>
#include <vector>
#include <algorithm>
#include <map>
#include <utility>
#include <tuple>

#include <stxxl/vector>

#include "predefined.h"


struct Voxel {
    size_t W, H, D;
    size_t WH, WD, HD;

    uchar solidValue;

    stxxl::vector<uchar> grayscaleStack;
    stxxl::vector<uchar> phasesStack;
    stxxl::vector<long long> componentsStack;

    std::vector<long long> top, bot, left, right, back, front;

    long long n_labels;
    long long total_volume;
    std::vector<long long > labels;                     //need to create stxxl containers in future
    std::map<long long, long long> clusterSizesDistr;   //
    std::array<long long, 256> grayscaleHistogram;      //

    std::default_random_engine generator;
    std::bernoulli_distribution doSwap;

    Voxel(size_t W_, size_t H_, size_t D_, uchar solivValue_);
//    : W(W_), H(H_), D(D_), labels(1, 0) {
//        grayscaleHistogram.fill(0);
//        WH = W_ * H_;
//        WD = W_ * D_;
//        HD = H_ * D_;
//    }

    long long getShift(size_t w, size_t h, size_t d);

    void getWSurface(size_t w, std::vector<long long> &surface);

    void getHSurface(size_t h, std::vector<long long> &surface);

    void getDSurface(size_t d, std::vector<long long> &surface);

    long long uf_find(long long  label);

    long long uf_union(long long one, long long another);
    long long uf_make_set();

    void calculatePercolation(const std::vector<long long> &oneLabels,
                              const std::vector<long long> &anotherLabels,
                              std::map<long long, double> &) const;

    long long clusterize();
    void clusterizeCentral();
    void unionOnEdges();
    long long getTotalVolume();
};


#endif
