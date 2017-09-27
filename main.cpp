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


int main(int argc, char *argv[])
{
    if (argc < 2) {
        std::cout << "No config path!" << std::endl;
    }
//    std::string dir("/home/efim/work/tomography_data/source");
    std::string dir(argv[1]);

    std::cout << "Config directory: " << dir << std::endl;

    handleBatch(dir);

//    QApplication a(argc, argv);
//    MainWindow w(voxel);
//    w.show();

//    return a.exec();
    return 0;
}
