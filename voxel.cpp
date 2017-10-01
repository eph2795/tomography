#include <array>
#include <set>
#include <vector>
#include <algorithm>
#include <map>
#include <utility>
#include <random>
#include <tuple>

#include <stxxl/vector>

#include "voxel.h"
#include "predefined.h"


Voxel::Voxel(size_t W_, size_t H_, size_t D_, uchar solidValue_)
    : W(W_), H(H_), D(D_), solidValue(solidValue_),
      n_labels(1), total_volume(0), labels(1, 0), doSwap(0.5),
      top(W_* H_), bot(W_ * H_),
      back(H_ * D_), front(H_ * D_),
      left(W_ * D_), right(W_ * D_) {
    grayscaleHistogram.fill(0);
    clusterSizesDistr[0] = 0;
    WH = W_ * H_;
    WD = W_ * D_;
    HD = H_ * D_;
}


long long Voxel::getShift(size_t w, size_t h, size_t d){
    return d * WH + h * W + w;
}


void Voxel::getWSurface(size_t w, std::vector<long long> &surface) {
    surface.resize(HD);
//    for (size_t d = 0; d < D; d++) {
//        std::copy(componentsStack.begin() + getShift(w, 0, d),
//                  componentsStack.begin() + getShift(w + 1, 0, d),
//                  surface.begin() + d * H);
//    }
    for (size_t d = 0; d < D; d++) {
        for (size_t h = 0; h < H; h++) {
            surface[d * H + h] = componentsStack[getShift(w, h, d)];
        }
    }
}


void Voxel::getHSurface(size_t h, std::vector<long long> &surface) {
    surface.resize(WD);
    for (size_t d = 0; d < D; d++) {
        for (size_t w = 0; w < W; w++) {
            surface[d * W + w] = componentsStack[getShift(w, h, d)];
        }
    }
}


void Voxel::getDSurface(size_t d, std::vector<long long> &surface) {
    surface.resize(WH);
//    std::copy(componentsStack.begin() + getShift(0, 0, d),
//              componentsStack.begin() + getShift(0, 0, d + 1),
//              surface.begin());
    for (size_t h = 0; h < H; h++) {
        for (size_t w = 0; w < W; w++) {
            surface[h * W + w] = componentsStack[getShift(w, h, d)];
        }
    }
}


long long Voxel::uf_find(long long label) {
    long long rootLabel = label;
    while (labels[rootLabel] != rootLabel) {
        rootLabel = labels[rootLabel];
    }

    while (labels[label] != label) {
        long long parentLabel = labels[label];
        labels[label] = rootLabel;
        clusterSizesDistr[rootLabel] += clusterSizesDistr[label];
        clusterSizesDistr[label] = 0;
        label = parentLabel;
    }
    return rootLabel;
}


long long Voxel::uf_union(long long one, long long another) {
    long long oneRoot = uf_find(one), anotherRoot = uf_find(another);
    return labels[oneRoot] = anotherRoot;
}


long long Voxel::uf_make_set() {
    labels.push_back(n_labels++);
    clusterSizesDistr.insert(std::pair<long long, long long>(n_labels - 1, 1));
    return n_labels - 1;
}


long long Voxel::clusterize() {
    std::cout << "Begin clustering!" << std::endl;
    clusterizeCentral();

    std::cout << "End clustering!" << std::endl;
    std::cout << "Number of clusters: " << n_labels << std::endl;

    std::map<long long, long long> labels_mapping;
    n_labels = 1;
    labels_mapping[0] = 0;

    std::vector<long long> current(WH);
    for (size_t d = 0; d < D; d++) {
        std::copy(componentsStack.begin() + getShift(0, 0, d),
                  componentsStack.begin() + getShift(0, 0, d + 1),
                  current.begin());

        for (size_t h = 0; h < H; h++) {
            for (size_t w = 0; w < W; w++) {
                long long root = uf_find(current[h * W + w]);

                if (labels_mapping.find(root) == labels_mapping.end()) {
                    labels_mapping[root] = n_labels++;
                }
                current[h * W + w] = labels_mapping[root];

                if (w == 0) {
                    left[d * W + w] = current[h * W + w];
                }
                if (w == W - 1) {
                    right[d * W + w] = current[h * W + w];
                }

            }
            std::copy(current.begin() + getShift(0, 0, 0),
                      current.begin() + getShift(0, 1, 0),
                      back.begin() + getShift(0, h, 0));
            std::copy(current.begin() + getShift(0, H - 1, 0),
                      current.begin() + getShift(0, H, 0),
                      front.begin() + getShift(0, h, 0));
        }
        if (d == 0) {
            top = current;
        }
        if (d == D - 1) {
            bot = current;
        }
        std::copy(current.begin(), current.end(), componentsStack.begin() + getShift(0, 0, d));
    }

    std::map<long long, long long> new_dist;
    for (auto it = labels_mapping.begin(); it != labels_mapping.end(); it++) {
        new_dist[it->second] =clusterSizesDistr[it->first];
    }
    clusterSizesDistr = new_dist;

    labels.resize(n_labels);
    std::iota(labels.begin(), labels.end(), 0);

    return n_labels;
}


bool is_valid(size_t w, size_t h, size_t d, int dw, int dh, int dd,
              size_t W, size_t H, size_t D,
              const std::vector<uchar> &current_phases,
              const std::vector<uchar> &prev_phases) {
    if ((w == 0) && (dw == -1) || (w == W - 1) && (dw == 1)) {
        return false;
    }
    if ((h == 0) && (dh == -1) || (h == H - 1) && (dh == 1)) {
        return false;
    }
    if ((d == 0) && (dd == -1) || (d == D - 1) && (dd == 1)) {
        return false;
    }
    size_t neighbor;
    if (dd == 0) {
        neighbor = current_phases[(h + dh) * W + (w + dw)];
    } else {
        neighbor = prev_phases[(h + dh) * W + (w + dw)];
    }
    return current_phases[h * W + w] == neighbor;
}

void get_valid_neighbors(size_t w, size_t h, size_t d, size_t W, size_t H, size_t D,
                         const std::vector<uchar> &current_phases,
                         const std::vector<uchar> &prev_phases,
                         const std::vector<long long> &current_components,
                         const std::vector<long long> &prev_components,
                         std::vector<long long> &valid_neighbors) {
    valid_neighbors.resize(0);
    if (is_valid(w, h, d, -1, 0, 0, W, H, D, current_phases, prev_phases)) {
        valid_neighbors.push_back(current_components[h * W + (w - 1)]);
    }
    if (is_valid(w, h, d, 0, -1, 0, W, H, D, current_phases, prev_phases)) {
        valid_neighbors.push_back(current_components[(h - 1) * W + w]);
    }
    if (is_valid(w, h, d, 0, 0, -1, W, H, D, current_phases, prev_phases)) {
        valid_neighbors.push_back(prev_components[h * W + w]);
    }
}

void Voxel::clusterizeCentral() {
    componentsStack.resize(W * H * D);
    std::cout << "Memory for clustering allocated!" << std::endl;

    std::vector<uchar> current_phases(WH), prev_phases(WH);
    std::vector<long long> current_components(WH), prev_components(WH);
    for (size_t d = 0; d < D; d++) {
        std::copy(phasesStack.begin() + getShift(0, 0, d),
                  phasesStack.begin() + getShift(0, 0, d + 1), current_phases.begin());
        std::copy(componentsStack.begin() + getShift(0, 0, d),
                  componentsStack.begin() + getShift(0, 0, d + 1), current_components.begin());
        for (size_t h = 0; h < H; h++) {
            for (size_t w = 0; w < W; w++) {
                long long cur_shift = h * W + w;
                uchar cur = current_phases[cur_shift];
                if (cur == solidValue) {
                    std::vector <long long> valid_neighbors;
                    get_valid_neighbors(w, h, d, W, H, D, current_phases, prev_phases,
                                        current_components, prev_components,
                                        valid_neighbors);

                    if (valid_neighbors.size() != 0) {
                        long long parent = valid_neighbors[0];
                        current_components[cur_shift] = parent;
                        clusterSizesDistr[parent] += 1;

                        for (size_t i = 1; i < valid_neighbors.size(); i++) {
                            long long child = valid_neighbors[i];
                            if (doSwap(generator)) {
                                std::swap(child, parent);
                            }
                            uf_union(child, parent);
                        }
                    } else {
                        long long parent = uf_make_set();
                        current_components[cur_shift] = parent;
                    }
                } else {
                    current_components[cur_shift] = 0;
                    clusterSizesDistr[0] += 1;
                }
            }
        }
        std::copy(current_components.begin(), current_components.end(),
                  componentsStack.begin() + getShift(0, 0, d));
        prev_components = current_components;
        prev_phases = current_phases;
    }
}


long long Voxel::getTotalVolume() {
    total_volume = 0;
    for (auto it = clusterSizesDistr.begin(); it != clusterSizesDistr.end(); it++) {
        total_volume += it->second;
    }
    return total_volume;
}


void Voxel::calculatePercolation(const std::vector<long long> &oneLabels,
                                 const std::vector<long long> &anotherLabels,
                                 std::map<long long, double> &clustersStatistic) const {
    std::set<long long> one_labels_set(oneLabels.cbegin(), oneLabels.cend());
    std::set<long long> another_labels_set(anotherLabels.cbegin(), anotherLabels.cend());

    long long percolation_volume = 0;
    long long clusters_num = 0;

    for (auto it = one_labels_set.cbegin(); it != one_labels_set.cend(); it++) {
        if ((another_labels_set.find(*it) != another_labels_set.cend()) && (*it != 0)) {
//            std::cout << " " << *it << " " << clusters_num << std::endl;
            clusters_num += 1;
            percolation_volume += clusterSizesDistr.at(*it);
            clustersStatistic[*it] = 0;
        }
    }

    for (auto it = clustersStatistic.begin(); it != clustersStatistic.end(); it++) {
        it->second = 1.0 * clusterSizesDistr.at(it->first) / percolation_volume;
    }
}
