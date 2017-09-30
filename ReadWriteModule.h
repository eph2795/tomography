#ifndef READWRITEMODULE_H
#define READWRITEMODULE_H

#include <string>
#include <array>

#include <stxxl/vector>

#include "predefined.h"
#include "voxel.h"

unsigned char ToBinary255(unsigned char);

int ReadImageStackFromBinary(const std::string &, size_t, size_t, size_t,
                             stxxl::vector<uchar> &, bool);

int ReadImageStackFromDirectory(const std::string &,
                                bool,
                                size_t &,
                                size_t &,
                                size_t &,
                                stxxl::vector<uchar> &,
                                std::array<long long, 256>&);

int WriteImageStackToBinary(const std::string &,
                            size_t,
                            size_t,
                            size_t,
                            const stxxl::vector<uchar> &,
                            bool);

int WriteGrayscaleToDirectory(const std::string &, size_t, size_t, size_t,
                              const stxxl::vector<uchar> &);

int WriteBinaryToDirectory(const std::string &, size_t, size_t, size_t,
                           const stxxl::vector<uchar> &);

int WriteComponentsToDirectory(const std::string &, size_t, size_t, size_t,
                               const stxxl::vector<long long> &,
                               const std::map<long long, long long> &);

void WriteCharactiristicsToCsv(const std::string &,
                               const std::string &,
                               const std::string &,
                               const std::string &,
                               const std::string &,
                               const std::string &,
                               Voxel &,
                               bool, bool, bool, bool,
                               double);

void handleBatch(const std::string &,
                 const std::string &);

#endif
