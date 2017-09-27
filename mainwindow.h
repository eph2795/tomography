#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>

#include "predefined.h"
#include "voxel.h"

namespace Ui {
class MainWindow;
}

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    explicit MainWindow(const Voxel &, QWidget *parent = 0);
    ~MainWindow();

private:
    Ui::MainWindow *ui;
};

#endif // MAINWINDOW_H
