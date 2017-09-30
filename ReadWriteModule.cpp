#include <string>
#include <iostream>
#include <math.h>
#include <array>
#include <random>

#include "stxxl/vector"
#include <boost/filesystem.hpp>
#include "boost/filesystem/operations.hpp"
#include "boost/filesystem/path.hpp"
#include "boost/progress.hpp"
#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/contrib/contrib.hpp>

#include "ReadWriteModule.h"
#include "predefined.h"
#include "threshold.h"
#include "MRF.h"

bool match_pattern(const std::string &s, const std::string &p) {
    if ((s.size() >= p.size())
        && (s.substr(s.size() - p.size(), p.size()) == p)) {
        return true;
    }
    return false;
}

uchar ToBinary255(uchar a)
{
	return a != 0 ? 255 : 0;
}


int ReadImageStackFromBinary(const std::string &FileName, size_t W, size_t H, size_t D,
                             stxxl::vector<uchar> &Voxel, bool DoTransform) {
    std::ifstream I(FileName.c_str(), std::ifstream::binary);
    if (I.bad()) {
        std::cout << "Opening failed: " << FileName << std::endl;
        return -1;
    }

    size_t N = W * H;
    char *buff = new char[N];
    Voxel.resize(W * H * D);

    size_t pos = 0;
    for (size_t i = 0; i < D; i++) {
        I.seekg(pos, std::ios::beg);
		I.read(buff, N);
        if (DoTransform) {
            std::transform(buff, buff + N, Voxel.begin() + pos, ToBinary255);
        } else {
            std::copy(buff, buff + N, Voxel.begin() + pos);
        }
        pos += N;
    }

	if (!I) {
		std::cout << "Reading failed: " << FileName << std::endl;
        pos = -1;
	}
	I.close();
	delete []buff;

	return pos;
}

int WriteImageStackToBinary(const std::string &FileName,
                            size_t W,
                            size_t H,
                            size_t D,
                            const stxxl::vector<uchar> &Voxel,
                            bool DoTransform) {
    std::cout << " WRITING " << FileName
              << ", size: " << Voxel.size()
              << ", " << W << " " << H << " " << D << std::endl;

    std::ofstream O(FileName.c_str(), std::ofstream::binary);
    if (O.bad()) {
        std::cout << "Opening failed: " << FileName << std::endl;
        return -1;
    }

    size_t N = W * H;
    char *buff = new char[N];

    size_t pos = 0;
    for (size_t i = 0; i < D; i++) {
        O.seekp(pos, std::ios::beg);
        if (DoTransform) {
            std::transform(Voxel.begin() + pos, Voxel.begin() + pos + N, buff, ToBinary255);
        } else {
            std::copy(Voxel.begin() + pos, Voxel.begin() + pos + N, buff);
        }
        O.write(buff, N);
        pos += N;
    }

    O.close();
    delete []buff;

    return pos;
}

bool isImage(const std::string &fname) {
    return match_pattern(fname, "bmp") && !match_pattern(fname, "_spr.bmp");
}

void getImageStackSizes(const std::vector<boost::filesystem::path> &directory_list,
                        size_t &W, size_t &H, size_t &D) {
    std::string fname = "@";
    for (const boost::filesystem::path &item : directory_list) {
        fname = item.c_str();
        if (isImage(fname)) {
            D += 1;
            if ((W == 0) && (H == 0)) {
                cv::Mat image = cv::imread(fname, CV_LOAD_IMAGE_GRAYSCALE);
                H = image.rows;
                W = image.cols;
            }
        }
    }
    if (fname == "@") {
        std::cout << "No image in folder!" << std::endl;
    }
}

int ReadImageStackFromDirectory(const std::string &PathToDir,
                                bool isBinary,
                                size_t &W, size_t &H, size_t &D,
                                stxxl::vector<uchar> &grayscaleStack,
                                std::array<long long, 256> &grayscaleHistogram) {
    size_t pos = 0;
    boost::filesystem::path path_to_dir(PathToDir);
    grayscaleHistogram.fill(0);
    try  {
        if (boost::filesystem::exists(path_to_dir) && boost::filesystem::is_directory(path_to_dir)) {
            std::vector<boost::filesystem::path> directory_list;
            std::copy(boost::filesystem::directory_iterator(path_to_dir),
                      boost::filesystem::directory_iterator(), std::back_inserter(directory_list));
            std::sort(directory_list.begin(), directory_list.end());
            if (((W != 0) || (H != 0) || (D != 0)) && (directory_list.size() != D)) {
                std::cout << " folder contains wrong number of images: "
                          << directory_list.size() << std::endl;
                return -1;
            }

            if ((W == 0) && (H == 0) && (D == 0)) {
                getImageStackSizes(directory_list, W, H, D);
                if ((W == 0) && (H == 0) && (D == 0)) {
                    std::cout << "Cant read stack sizes: " << PathToDir << std::endl;
                    return -1;
                }
                std::cout << "Readen sizes: " << W << ", " << H << ", " << D << std::endl;
            }

            grayscaleStack.resize(W * H * D);
            size_t N = W * H;

            for (const boost::filesystem::path &item : directory_list) {
                if (!isImage(item.c_str())) {
                    continue;
                }
                cv::Mat image;
                image = cv::imread(item.c_str(), CV_LOAD_IMAGE_GRAYSCALE);

                if ((image.rows != H) || (image.cols != W) || (image.dims != 2)) {
                    std::cout << " folder contains image with wrong sizes: "
                              << image.rows << ", " << image.cols << ", " << image.dims << std::endl;
                    return -1;
                } else {
                    if (isBinary) {
                        std::transform(image.data, image.data + N, grayscaleStack.begin() + pos,
                                       [&grayscaleHistogram](uchar c)
                                       { return grayscaleHistogram[c]++, (uchar)(c != 0); });
                    } else {
                        std::transform(image.data, image.data + N, grayscaleStack.begin() + pos,
                                       [&grayscaleHistogram](uchar c)
                                       { return grayscaleHistogram[c]++, c; });
                    }
                    pos += N;
                }
            }
        } else {
            std::cout << path_to_dir << " does not exist or exists, but is not a directory"
                      << std::endl;
            return -1;
        }
    } catch (const boost::filesystem::filesystem_error& ex) {
        std::cout << ex.what() << std::endl;
    }
    return pos;
}


int WriteBinaryToDirectory(const std::string &PathToDir, size_t W, size_t H, size_t D,
                           const stxxl::vector<uchar> &grayscaleStack) {

    if (grayscaleStack.size() != W * H * D) {
        std::cout << " wrong voxel sizes" << std::endl;
        return -1;
    }

    boost::filesystem::path path_to_dir(PathToDir);
    size_t N = W * H;
    size_t pos = 0;
    uchar *buf = new uchar[N];

    for (size_t i = 0; i < D; i++) {
        std::transform(grayscaleStack.begin() + pos, grayscaleStack.begin() + pos + N, buf,
                       ToBinary255);
        cv::Mat image = cv::Mat(H, W, CV_8UC1, buf);

        boost::filesystem::path current_path = path_to_dir / ("im" + std::to_string(i) + ".bmp");
        cv::imwrite(current_path.c_str(), image);
        pos += N;
    }
    return pos;
}


int WriteGrayscaleToDirectory(const std::string &PathToDir, size_t W, size_t H, size_t D,
                              const stxxl::vector<uchar> &grayscaleStack) {
    if (grayscaleStack.size() != W * H * D) {
        std::cout << " wrong voxel sizes" << std::endl;
        return -1;
    }

    boost::filesystem::path path_to_dir(PathToDir);
    size_t N = W * H;
    size_t pos = 0;
    uchar *buf = new uchar[N];

    auto max_label_it = std::max_element(grayscaleStack.begin(), grayscaleStack.end());
    auto transform_function =[max_label_it](uchar c) { return uchar(1.0 * c / *max_label_it * 255); };
    for (size_t i = 0; i < D; i++) {
        std::transform(grayscaleStack.begin() + pos, grayscaleStack.begin() + pos + N, buf,
                       transform_function);
        cv::Mat image = cv::Mat(H, W, CV_8UC1, buf);

        boost::filesystem::path current_path = path_to_dir / ("im" + std::to_string(i) + ".bmp");
        cv::imwrite(current_path.c_str(), image);
        pos += N;
    }
    return pos;
}


int WriteComponentsToDirectory(const std::string &PathToDir, size_t W, size_t H, size_t D,
                               const stxxl::vector<long long> &components,
                               const std::map<long long, long long> &sizesDistr) {
    if (components.size() != W * H * D) {
        std::cout << " wrong voxel sizes" << std::endl;
        return -1;
    }

    std::vector<size_t> colors_ids(COMPONENTS_NUMBER + 1);
    std::iota(colors_ids.begin(), colors_ids.end(), 0);
    std::random_shuffle(colors_ids.begin() + 1, colors_ids.end());

    std::vector<std::pair<long long, long long>>
            component_to_color(sizesDistr.cbegin(), sizesDistr.cend());
    std::sort(component_to_color.begin(), component_to_color.end(),
              [](const std::pair<long long, long long> &one,
                 const std::pair<long long, long long> &another)
                { return one.second > another.second; });

    for (size_t i = 0; i < component_to_color.size(); i++) {
        if (component_to_color[i].first == 0) {
            component_to_color[i].second = 0;
        } else {
            if (i < COMPONENTS_NUMBER) {
                component_to_color[i].second = colors_ids[i + 1];
            } else {
                component_to_color[i].second = colors_ids[COMPONENTS_NUMBER];
            }
        }
    }

    std::map<long long, size_t> color_map;
    std::copy(component_to_color.begin(), component_to_color.end(),
              std::inserter(color_map, color_map.begin()));

//    for(auto& p : component_to_color) {
//        color_map[p.first] = p.second;
//    }

    boost::filesystem::path path_to_dir(PathToDir);
    size_t N = W * H;
    size_t pos = 0;
    uchar *buf = new uchar[3 * N];

    for (size_t i = 0; i < D; i++) {
        for (auto it = components.begin() + pos;
                  it != components.begin() + pos + N; it++) {
            size_t j = it - components.begin() - pos;
            buf[3 * j] = COLORS[color_map[*it]][0];
            buf[3 * j + 1] = COLORS[color_map[*it]][1];
            buf[3 * j + 2] = COLORS[color_map[*it]][2];
        }
        cv::Mat image = cv::Mat(H, W, CV_8UC3, buf);

        boost::filesystem::path current_path = path_to_dir / ("im" + std::to_string(i) + ".bmp");
        cv::imwrite(current_path.c_str(), image);
        pos += N;
    }

    return pos;
}


void WriteCharactiristicsToCsv(const std::string &pathToDir,
                               const std::string &pathToCSV,
                               const std::string &suffix,
                               const std::string &percXY,
                               const std::string &percXZ,
                               const std::string &percYZ,
                               Voxel &voxel,
                               bool isAbsolute, bool isXY, bool isXZ, bool isYZ,
                               double volumeThreshold) {
    voxel.getTotalVolume();

    std::vector<long long> top, bot, left, right, back, front;
    voxel.getDSurface(0, top);
    voxel.getDSurface(voxel.D - 1, bot);
    voxel.getHSurface(0, left);
    voxel.getHSurface(voxel.H - 1, right);
    voxel.getWSurface(0, back);
    voxel.getWSurface(voxel.W - 1, front);

    std::cout << "Computing percolations... " << std::endl;

    std::map<long long, double> statisticsXY, statisticsXZ, statisticsYZ;
    voxel.calculatePercolation(bot, top, statisticsXY);
    voxel.calculatePercolation(left, right, statisticsXZ);
    voxel.calculatePercolation(back, front, statisticsYZ);

    std::ofstream output(pathToCSV + "/statistics_" + suffix + ".csv",
                         std::ofstream::in | std::ofstream::binary | std::ofstream::app);
    output << std::fixed << std::setprecision(6);

    boost::filesystem::path path_to_dir(pathToDir);
    output << path_to_dir.c_str() << ", \n";

    output << "ID, Relative volume";
    if (isAbsolute) {
        output << ", Absolute volume";
    }
    if (isXY) {
        output << ", Percolation XY, Relative XY volume";
    }
    if (isXZ) {
        output << ", Percolation XZ, Relative XZ volume";
    }
    if (isYZ) {
        output << ", Percolation YZ, Relative YZ volume";
    }
    output << "\n";

    std::vector<long long> order(voxel.n_labels);
    std::iota(order.begin(), order.end(), 0);
    std::sort(order.begin(), order.end(),
              [&voxel](long long one, long long another)
              { return voxel.clusterSizesDistr[one] > voxel.clusterSizesDistr[another]; });

    bool is_xy = false
       , is_xz = false
       , is_yz = false;
    for (long long  i : order) {
        if (i == 0) {
            continue;
        }

        long long cur_volume = voxel.clusterSizesDistr[i];
        double cur_rel_volume = 1.0 * voxel.clusterSizesDistr[i] / voxel.W / voxel.H / voxel.D;
        if (cur_rel_volume < volumeThreshold) {
            break;
        }

        output << i << ", " << cur_rel_volume;
        if (isAbsolute) {
            output << ", " << cur_volume;
        }

        if (isXY) {
            if (statisticsXY.find(i) != statisticsXY.end()) {
                output << ", true," << statisticsXY[i];
                if (!is_xy) {
                    WriteImageStackToBinary(percXY, voxel.W, voxel.H, voxel.D, voxel.grayscaleStack, true);
                    is_xy = true;
                }
            } else {
                output << ", false, -";
            }
        }
        if (isXZ) {
            if (statisticsXZ.find(i) != statisticsXZ.end()) {
                if (!is_xz) {
                    WriteImageStackToBinary(percXZ, voxel.W, voxel.H, voxel.D, voxel.grayscaleStack, true);
                    is_xz = true;
                }
                output << ", true, " << statisticsXZ[i];
            } else {
                output << ", false, -";
            }
        }
        if (isYZ) {
            if (statisticsYZ.find(i) != statisticsYZ.end()) {
                if (!is_yz) {
                    WriteImageStackToBinary(percYZ, voxel.W, voxel.H, voxel.D, voxel.grayscaleStack, true);
                    is_yz = true;
                }
                output << ", true, " << statisticsYZ[i];
            } else {
                output << ", false, -";
            }
        }
        output << "\n";
    }
    output.close();
}


void readConfig(const std::string &pathToConfig,
                MRFSettings &settings,
                Threshold &threshold,
                bool &isBinaryData,
                bool &doMRF,
                uchar &solidValue,
                bool &isAbsolute, bool &isXY, bool &isXZ, bool &isYZ,
                bool &produceBinary,
                bool &produceComponents) {
    std::ifstream config(pathToConfig, std::ios::in);

    size_t nPhases = 2;
    std::vector<unsigned char> Low, High;

    double mBeta, mDeltaT, mT0;

    if (config.is_open()) {
        std::string line;

        std::getline(config, line);
        std::istringstream phases(line);
        phases >> nPhases;

        Low.resize(nPhases);
        High.resize(nPhases);

        std::getline(config, line);
        std::istringstream is_binary_data(line);
        is_binary_data >> isBinaryData;
        if (!isBinaryData) {
            std::getline(config, line);
            std::istringstream lows(line);
            size_t v;
            for (size_t i = 0; i < nPhases; i++) {
                lows >> v;
                Low[i] = v;
            }

            std::getline(config, line);
            std::istringstream highs(line);
            for (size_t i = 0; i < nPhases; i++) {
                highs >> v;
                High[i] = v;
            }
        }

        std::getline(config, line);
        std::istringstream do_mrf(line);
        do_mrf >> doMRF;
        if (doMRF) {
            std::getline(config, line);
            std::istringstream params(line);
            params >> mBeta >> mDeltaT >> mT0;
        }

        std::getline(config, line);
        std::istringstream solid(line);
        size_t value;
        solid >> value;
        solidValue = value;

        std::getline(config, line);
        std::istringstream flags(line);
        flags >> isAbsolute >> isXY >> isXZ >> isYZ;

        std::getline(config, line);
        std::istringstream produce_binary(line);
        produce_binary >> produceBinary;

        std::getline(config, line);
        std::istringstream produce_components(line);
        produce_components >> produceComponents;

        settings.Beta = mBeta;
        settings.TStart = mT0;
        settings.FreezingSpeed = mDeltaT;
        settings.Method = MRF_SA;

        std::cout << "Number of phases: " << nPhases << std::endl;
        std::cout << "Data is divided into phases: " << isBinaryData << std::endl;
        if (!isBinaryData) {
            std::cout << "Phases intervals: ";
            for (size_t i = 0; i < nPhases; i++) {
                std::cout << "[" << (size_t)Low[i] << ", " << (size_t)High[i] << "] ";
            }
            std::cout << std::endl;
        }

        std::cout << "Do MRF: " << doMRF << std::endl;
        if (doMRF) {
            std::cout << "MRF params: ";
            std::cout << "mBeta: " << mBeta << ", mDeltaT: " << mDeltaT << ",  mT0: " << mT0 << std::endl;
        }

        std::cout << "Number of solid phase: " << (size_t)solidValue << std::endl;
        std::cout << "Print absolute values: " << (bool)isAbsolute << std::endl;
        std::cout << "Print XY percolation:  " << (bool)isXY << std::endl;
        std::cout << "Print XZ percolation:  " << (bool)isXZ << std::endl;
        std::cout << "Print YZ percolation:  " << (bool)isYZ << std::endl;

        std::cout << "Produce binarization: " << (bool)produceBinary << std::endl;
        std::cout << "Produce clusterization components: " << (bool)produceComponents << std::endl;

        config.close();
    } else {
        std::cout << "Bad file!" << std::endl;
    }
    threshold = Threshold(nPhases, Low, High);
}


void handleDir(const std::string &targetPath,
               const std::string &pathToDir,
               const MRFSettings &settings,
               const Threshold &threshold,
               bool is_binary_data,
               bool do_mrf,
               uchar solid_value,
               bool is_absolute, bool is_xy, bool is_xz, bool is_yz,
               bool produce_binary,
               bool produce_components) {
    boost::filesystem::path path_to_dir(pathToDir);
    size_t W = 0, H = 0, D = 0;
    stxxl::vector<uchar> grayscaleStack;
    std::array<long long, 256> grayscaleHistogram;

    int readen_bytes = ReadImageStackFromDirectory(pathToDir,
                                                   is_binary_data,
                                                   W, H, D,
                                                   grayscaleStack,
                                                   grayscaleHistogram);
    if (D == 0) {
        return;
    }
    std::cout << "Current folder: " << pathToDir << std::endl;

    std::cout << "Readen bytes: " << readen_bytes << std::endl;

    Voxel voxel(W, H, D, solid_value);
    voxel.grayscaleStack = grayscaleStack;
    voxel.grayscaleHistogram = grayscaleHistogram;
    grayscaleStack.resize(0);

    MRF mrf(settings, voxel.W, voxel.H, voxel.D, threshold);

    std::clock_t begin = clock();

    if (!is_binary_data) {
        mrf.StatsNL(voxel.grayscaleStack);
        mrf.ConditionalImage(voxel.grayscaleStack, voxel.phasesStack);

        if (do_mrf) {
            std::cout << "Starting MRF..." << std::endl;
            mrf.SimulatedAnnealing3D(voxel.grayscaleStack, voxel.phasesStack);
        }
    } else {
       voxel.phasesStack = voxel.grayscaleStack;
    }

    std::string target_name;
    boost::filesystem::path target_path;

    if (produce_binary) {
        target_name = pathToDir + "_threshold";
        target_path = boost::filesystem::path(target_name);
        if (!boost::filesystem::exists(target_path) ||
            !boost::filesystem::is_directory(target_path)) {
            boost::filesystem::create_directory(target_name);
        }
        for (auto it = boost::filesystem::directory_iterator(target_path);
                  it != boost::filesystem::directory_iterator(); it ++) {
            boost::filesystem::remove_all(it->path());
        }
        WriteBinaryToDirectory(target_name,
                               voxel.W, voxel.H, voxel.D, voxel.phasesStack);
    }

    std::clock_t end = clock();
    double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
    begin = end;

    std::cout << "Thresholding and binary writing: " << elapsed_secs << std::endl;

    voxel.clusterize();

    end = clock();
    elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
    begin = end;

    std::cout << "Full clusterization: " << elapsed_secs << std::endl;

    if (produce_components) {
        target_name = pathToDir + "_result";
        target_path = boost::filesystem::path(target_name);
        if (!boost::filesystem::exists(target_path) ||
            !boost::filesystem::is_directory(target_path)) {
            boost::filesystem::create_directory(target_name);
        }
        for (auto it = boost::filesystem::directory_iterator(target_path);
                  it != boost::filesystem::directory_iterator(); it ++) {
            boost::filesystem::remove_all(it->path());
        }

        WriteComponentsToDirectory(target_name,
                                   voxel.W, voxel.H, voxel.D, voxel.componentsStack,
                                   voxel.clusterSizesDistr);
    }

    WriteCharactiristicsToCsv(pathToDir,
                              targetPath,
                              path_to_dir.filename().c_str(),
                              targetPath + "/percXY/" + path_to_dir.filename().c_str(),
                              targetPath + "/percXZ/" + path_to_dir.filename().c_str(),
                              targetPath + "/percYZ/" + path_to_dir.filename().c_str(),
                              voxel, is_absolute, is_xy, is_xz, is_yz,
                              1e-6);

    end = clock();
    elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;

   std::cout << "Colored writing and percolation calc: "
             << elapsed_secs << std::endl;
   std::cout << "Leaving folder " << pathToDir << std::endl << std::endl;
}

bool isDataDirectory(const std::string &directory) {
    bool result = false;
    if (((directory.size() < std::string("result").size()) || !match_pattern(directory, "result"))
        && ((directory.size() < std::string("threshold").size()) || !match_pattern(directory, "threshold"))
        && ((directory.size() < std::string("percXY").size()) || !match_pattern(directory, "percXY"))
        && ((directory.size() < std::string("percXZ").size()) || !match_pattern(directory, "percXZ"))
        && ((directory.size() < std::string("percYZ").size()) || !match_pattern(directory, "percYZ"))) {
        result = true;
    }

    return result;
}


bool produceTruncatedStack(const std::string &pathToDir,
                           size_t biasW, size_t biasH, size_t biasD,
                           size_t newW, size_t newH, size_t newD,
                           std::string &newName) {
    size_t W = 0, H = 0, D = 0;
    stxxl::vector<uchar> grayscaleStack;

    std::array<long long, 256> grayscaleHistogram;

    int readen_bytes = ReadImageStackFromDirectory(pathToDir,
                                                   false,
                                                   W, H, D,
                                                   grayscaleStack,
                                                   grayscaleHistogram);
    if ((D == 0) || (D < newD) || (W < newW) || (H < newH)) {
        return false;
    }

    size_t N = W * H;
    size_t newN = newW * newH;
    stxxl::vector<uchar> newStack(newW * newH * newD);
    for (size_t d = biasD; d < biasD + newD; d++) {
        for (size_t h = biasH; h < biasH + newH; h++) {
            size_t shift = d * N + h * W + biasW;
            size_t shift_new = (d - biasD) * newN + (h - biasH) * newW;
            std::copy(grayscaleStack.begin() + shift, grayscaleStack.begin() + shift + newW, newStack.begin() + shift_new);
        }
    }

    boost::filesystem::path path_to_dir(pathToDir);
    std::string filename = path_to_dir.filename().c_str();
    boost::filesystem::path parent_path = path_to_dir.parent_path();
    newName = (parent_path
        / (filename
           + "_W" + std::to_string(newW)
           + "_H" + std::to_string(newH)
           + "_D" + std::to_string(newD))
    ).c_str();

    if (!boost::filesystem::exists(newName) ||
        !boost::filesystem::is_directory(newName)) {
        boost::filesystem::create_directory(newName);
    }
    WriteGrayscaleToDirectory(newName, newW, newH, newD, newStack);

    return true;
}


void DFS(const std::string currentDir,
         std::vector<std::string> &dirTree,
         bool produceTruncated) {
    boost::filesystem::path path_to_dir(currentDir);

    if (produceTruncated){
        std::string new_name;
        bool is_produced = produceTruncatedStack(currentDir,
                                                 0, 0, 0,
                                                 500, 500, 500,
                                                 new_name);

        std::cout << "Directory: " << currentDir
                  << ", produced: " << is_produced << std::endl;

        if (is_produced) {
            dirTree.push_back(new_name);
        }
    }

    std::vector<boost::filesystem::path> directory_list;
    std::copy(boost::filesystem::directory_iterator(path_to_dir),
              boost::filesystem::directory_iterator(), std::back_inserter(directory_list));

    for (const boost::filesystem::path &item : directory_list) {
        if (boost::filesystem::is_directory(item) && isDataDirectory(item.c_str())) {
            dirTree.push_back(item.c_str());
            std::cout << "     " << item.c_str() << std::endl;

            DFS(item.c_str(), dirTree, produceTruncated);
        }
    }
}


void handleBatch(const std::string &targetPath,
                 const std::string &pathToDir) {
    boost::filesystem::path path_to_dir(pathToDir);
    uchar solid_value;
    bool is_binary_data, do_mrf, is_absolute, is_xy, is_xz, is_yz, produce_binary, produce_components;
    try  {
        if (boost::filesystem::exists(path_to_dir) && boost::filesystem::is_directory(path_to_dir)) {
            Threshold threshold;
            MRFSettings settings;
            readConfig((path_to_dir / "config.txt").c_str(),
                       settings, threshold,
                       is_binary_data, do_mrf,
                       solid_value, is_absolute, is_xy, is_xz, is_yz,
                       produce_binary,
                       produce_components);

            std::vector<boost::filesystem::path> directory_list;
            std::copy(boost::filesystem::directory_iterator(path_to_dir),
                      boost::filesystem::directory_iterator(), std::back_inserter(directory_list));
            std::sort(directory_list.begin(), directory_list.end());

            std::vector<std::string> dirTree(0);

            DFS(pathToDir, dirTree, true);

            for (const std::string &item : dirTree) {
                handleDir(targetPath,
                          item,
                          settings,
                          threshold,
                          is_binary_data,
                          do_mrf,
                          solid_value,
                          is_absolute, is_xy, is_xz, is_yz,
                          produce_binary,
                          produce_components);
            }
        }
    } catch (const boost::filesystem::filesystem_error& ex) {
        std::cout << ex.what() << std::endl;
    }
}
