#include "mainwindow.h"
#include "ui_mainwindow.h"

#include "voxel.h"
#include "predefined.h"


void make_hist(size_t nbins, size_t bin_size,
               const QVector<double> &x,
               const QVector<double> &y,
               QCustomPlot *plot) {
    size_t y_max = *std::max_element(y.begin(), y.end());
    size_t y_min = *std::min_element(y.begin(), y.end());
    size_t x_max = *std::max_element(x.begin(), x.end());
    size_t x_min = *std::min_element(x.begin(), x.end());

    QCPBars *bars = new QCPBars(plot->xAxis, plot->yAxis);
    plot->yAxis->setRange(0.9 * y_min, 1.1 * y_max);

    QVector<double> ticks;
    QVector<QString> labels;

    for (size_t i = 1; i < nbins + 1; i++) {
        ticks << i;
        labels << QString::fromStdString(std::to_string(x_min + (i - 1) * bin_size) +
                                         "-" +
                                         std::to_string(i * bin_size));
    }

    QSharedPointer<QCPAxisTickerText> textTicker(new QCPAxisTickerText);
    textTicker->addTicks(ticks, labels);
    plot->xAxis->setTicker(textTicker);
    plot->xAxis->setTickLabelRotation(60);
    plot->xAxis->setSubTicks(false);
    plot->xAxis->setRange(0, nbins + 2);

    bars->setData(ticks, y);
    plot->legend->setVisible(true);
    return;
}

MainWindow::MainWindow(const Voxel &voxel, QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow)
{
    ui->setupUi(this);

    {
    size_t x_max = 255;
    size_t x_min = 0;
    size_t nbins = 16;
    size_t bin_size = (x_max - x_min + 1) / nbins + ((x_max - x_min + 1) % nbins != 0) ;

    QVector<double> x(nbins), y(nbins, 0.0);
    for (size_t i = x_min; i < x_max; i++) {
        y[(i - x_min) / bin_size] += voxel.grayscaleHistogram[i];
    }

    double y_max = 0, y_min =voxel.grayscaleHistogram[0];
    for (size_t i = 0; i < nbins; i++) {
        x[i] = i;
        if (y[i] > y_max) {
            y_max = y[i];
        }
        if (y[i] < y_min) {
            y_min = y[i];
        }
    }

    make_hist(nbins, bin_size, x, y, ui->thresholdHistogramPlot);
    }

    {
    std::map<long long, long long> data = voxel.clusterSizesDistr;
    long long x_max = data[0];
    long long x_min = data[0];
    for (auto it = data.begin(); it != data.end(); it++) {
        if (x_max < it->second) {
            x_max = it->second;
        }
        if (x_min > it->second) {
            x_min = it->second;
        }
    }

    size_t nbins = 10;
    size_t bin_size = (x_max - x_min + 1) / nbins + ((x_max - x_min + 1) % nbins != 0) ;

    QVector<double> x(nbins), y(nbins, 0.0);
    for (auto it = data.begin(); it != data.end(); it++) {
        if (it->first != 0) {
            y[(it->second - x_min) / bin_size] += it->second;
        }
    }

    double y_max = data[0], y_min = data[0];
    for (size_t i = 0; i < nbins; i++) {
        x[i] = i;
        if (y[i] > y_max) {
            y_max = y[i];
        }
        if (y[i] < y_min) {
            y_min = y[i];
        }
    }

    make_hist(nbins, bin_size, x, y, ui->clustersHistogramPlot);
    }
    // create graph and assign data to it:
//    QCPBars *bars = new QCPBars(ui->thresholdHistogramPlot->xAxis, ui->thresholdHistogramPlot->yAxis);
//    //bars ->setWidth(9/(double)x4.size());

//    ui->thresholdHistogramPlot->yAxis->setRange(0, 1.1 * y_max);
////    ui->thresholdHistogramPlot->yAxis->setScaleType(QCPAxis::stLogarithmic);

//    QVector<double> ticks;
//    QVector<QString> labels;
//    for (size_t i = 1; i < nbins + 1; i++) {
//        ticks << i;
//        labels << QString::fromStdString(std::to_string((i - 1) * bin_size) + "-" + std::to_string(std::min(i * bin_size, (size_t)256)));
//    }

//    QSharedPointer<QCPAxisTickerText> textTicker(new QCPAxisTickerText);
//    textTicker->addTicks(ticks, labels);
//    ui->thresholdHistogramPlot->xAxis->setTicker(textTicker);
//    ui->thresholdHistogramPlot->xAxis->setTickLabelRotation(60);
//    ui->thresholdHistogramPlot->xAxis->setSubTicks(false);
//    ui->thresholdHistogramPlot->xAxis->setRange(0, nbins + 2);

//    bars->setData(ticks, y);
//    ui->thresholdHistogramPlot->legend->setVisible(true);
}

MainWindow::~MainWindow()
{
    delete ui;
}
