#-------------------------------------------------
#
# Project created by QtCreator 2017-07-08T17:45:50
#
#-------------------------------------------------

QT       += core gui

greaterThan(QT_MAJOR_VERSION, 4): QT += widgets printsupport

TARGET = tomography
TEMPLATE = app

# The following define makes your compiler emit warnings if you use
# any feature of Qt which as been marked as deprecated (the exact warnings
# depend on your compiler). Please consult the documentation of the
# deprecated API in order to know how to port your code away from it.
DEFINES += QT_DEPRECATED_WARNINGS

# You can also make your code fail to compile if you use deprecated APIs.
# In order to do so, uncomment the following line.
# You can also select to disable deprecated APIs only up to a certain version of Qt.
#DEFINES += QT_DISABLE_DEPRECATED_BEFORE=0x060000    # disables all the APIs deprecated before Qt 6.0.0


INCLUDEPATH += /usr/include/opencv


LIBS += \
    -L /usr/lib/cmake -l stxxl \
    -lboost_date_time \
    -lboost_filesystem \
    -lboost_program_options \
    -lboost_regex \
    -lboost_signals \
    -lboost_system \
    -lopencv_core \
    -lopencv_imgproc \
    -lopencv_highgui \
    -lopencv_ml \
    -lopencv_video \
    -lopencv_features2d \
    -lopencv_calib3d \
    -lopencv_objdetect \
    -lopencv_contrib \
    -lopencv_legacy \
    -lopencv_flann \


SOURCES += \
        main.cpp \
        mainwindow.cpp \
        qcustomplot.cpp \
        voxel.cpp \
        ReadWriteModule.cpp \


HEADERS += \
        mainwindow.h \
        LatticeModel.h \
        MRF.h \
        threshold.h \
        voxel.h \
        predefined.h \
        qcustomplot.h \
        ReadWriteModule.h \

FORMS += \
        mainwindow.ui \
