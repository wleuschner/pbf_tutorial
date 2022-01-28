#ifndef CANVAS_H
#define CANVAS_H

#include <QWidget>
#include <QTimer>
#include <glm/glm.hpp>
#include "abstractsolver.h"

class Canvas : public QWidget
{
    Q_OBJECT
public:
    explicit Canvas(QWidget *parent = nullptr);

protected:
    void resizeEvent(QResizeEvent *event) override;
    void paintEvent(QPaintEvent *event) override;
    void mousePressEvent(QMouseEvent *event) override;
    void mouseMoveEvent(QMouseEvent *event) override;
    void mouseReleaseEvent(QMouseEvent *event) override;
public slots:
    void toggleSimulation();
    void stepSimulation();
    void togglePressureView(bool value);
    void changeParticleSize(double value);
    void changeSearchRadius(double value);
    void changeRestingDensity(double value);
    void changeArtificialDensity(double value);
    void changeViscosity(double value);
    void changeTimestep(double value);
    void changeIterations(int niter);
    void changeSpatialStruct(int index);
    void setGPU(bool gpu);
private slots:
    void simulate();
private:
    bool showPressure;
    float particleSize;
    QTimer updateTimer;

    AbstractSolver* solver;
    glm::vec2 line[2];
};

#endif // CANVAS_H
