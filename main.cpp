#include "mainwindow.h"
#include <QApplication>

#include <iostream>
#include <string>
#include <sstream>
#include <ctime>

#include <stxxl/vector>
#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>

#include "ReadWriteModule.h"
#include "LatticeModel.h"
#include "MRF.h"
#include "threshold.h"
#include "voxel.h"

#include <boost/filesystem.hpp>
#include "boost/filesystem/operations.hpp"
#include "boost/filesystem/path.hpp"
#include "boost/progress.hpp"

int main(int argc, char *argv[])
{
    std::clock_t begin = clock();

    if (argc < 2) {
        std::cout << "No config path!" << std::endl;
    }
//    std::string dir("/home/efim/work/tomography_data/source");

    std::string dir(argv[1]);

    std::cout << "Config directory: " << dir << std::endl;

    std::vector<std::string> perc = {"/percXY", "/percXZ", "/percYZ"};
    for (const std::string &suff: perc) {
        std::string name = dir + suff;
        boost::filesystem::path path = boost::filesystem::path(name);
        if (!boost::filesystem::exists(path) ||
            !boost::filesystem::is_directory(path)) {
            boost::filesystem::create_directory(name);
        }
    }

    handleBatch(dir, dir);

    std::clock_t end = clock();
    double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
    std::cout << "Total runtime: " << elapsed_secs << std::endl;
//    QApplication a(argc, argv);
//    MainWindow w(voxel);
//    w.show();

//    return a.exec();
    return 0;
}
