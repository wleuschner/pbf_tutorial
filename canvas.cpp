#include "canvas.h"
#include <QPainter>
#include <QMouseEvent>
#include "solver.h"
#include "solver_gpu.h"

Canvas::Canvas(QWidget *parent) : QWidget(parent)
{
    solver = new Solver(width(),height());
    particleSize = 10.0f;
    showPressure = false;

    connect(&updateTimer,SIGNAL(timeout()),this,SLOT(simulate()));
    updateTimer.setInterval(1000.0/60.0);
    updateTimer.setSingleShot(false);
}

void Canvas::toggleSimulation()
{
    if(updateTimer.isActive())
    {
        updateTimer.stop();
    }
    else
    {
        updateTimer.start();
    }
}

void Canvas::stepSimulation()
{
    solver->solve();
    update();
}

void Canvas::togglePressureView(bool value)
{
    showPressure = value;
    update();
}

void Canvas::changeParticleSize(double value)
{
    particleSize = value;
    update();
}

void Canvas::changeSearchRadius(double value)
{
    solver->setSearchRadius(value);
}

void Canvas::changeRestingDensity(double value)
{
    solver->setRestingDensity(value);
}

void Canvas::changeArtificialDensity(double value)
{
    solver->setArtificialDensity(value);
}

void Canvas::changeViscosity(double value)
{
    solver->setViscosity(value);
}

void Canvas::changeTimestep(double value)
{
    solver->setTimeStep(value);
}

void Canvas::changeIterations(int niter)
{
    solver->setIterations(niter);
}

void Canvas::changeSpatialStruct(int index)
{
    Solver::SpatialStructType type;
    switch (index)
    {
        case 0:
        {
            type = Solver::LINEAR_SEARCH;
            break;
        }
        case 1:
        {
            type = Solver::RADIX_SORT;
            break;
        }
    }
    solver->setSpatialStruct(type);
}

void Canvas::setGPU(bool gpu)
{
    std::vector<Particle> particles = solver->getParticles();
    std::vector<Boundary> boundaries = solver->getBoundaries();
    delete solver;
    if(gpu)
    {

        solver = new SolverGPU(width(),height());
    }
    else
    {
        solver = new Solver(width(),height());
    }
    solver->changeDomain(0.0,height()-100,width()-100,0.0f);
    solver->setParticles(particles);
    solver->setBoundaries(boundaries);
}

void Canvas::simulate()
{
    solver->solve();
    update();
}

void Canvas::resizeEvent(QResizeEvent *event)
{
    solver->changeDomain(0.0,height()-100,width()-100,0.0f);
}

void Canvas::paintEvent(QPaintEvent *event)
{
    QPainter painter(this);
    QBrush brush(QColor(0,0,255));
    brush.setStyle(Qt::BrushStyle::SolidPattern);
    painter.setPen(QColor(0,0,0));
    painter.setBrush(brush);

    painter.drawLine(line[0].x,line[0].y,line[1].x,line[1].y);

    const std::vector<Particle>& particles = solver->getParticles();
    const std::vector<Boundary>& boundaries = solver->getBoundaries();
    for(unsigned int i=0;i<boundaries.size();i++)
    {
        const Boundary& b = boundaries[i];
        painter.drawLine(b.line[0].x,height()-b.line[0].y,b.line[1].x,height()-b.line[1].y);
    }

    //Particle Draw Loop
    QColor color = QColor(0,0,255);
    brush = QBrush(color);
    brush.setStyle(Qt::BrushStyle::SolidPattern);
    painter.setPen(color);
    painter.setBrush(brush);

    for(unsigned int i=0;i<particles.size();i++)
    {
        const Particle& p = particles[i];
        if(showPressure)
        {
            color = (0,0,std::min(p.pressure * 1000000.0f,255.0f));
            brush = QBrush(color);
            brush.setStyle(Qt::BrushStyle::SolidPattern);
            painter.setPen(color);
            painter.setBrush(brush);
        }

        painter.drawEllipse(QPoint(p.pos.x,height()-p.pos.y),particleSize,particleSize);
    }
}

void Canvas::mousePressEvent(QMouseEvent *event)
{
    switch (event->button())
    {
        case Qt::MouseButton::LeftButton:
        {
            break;
        }
        case Qt::MouseButton::RightButton:
        {
            line[0] = glm::vec2(event->pos().x(), event->pos().y());
            line[1] = glm::vec2(event->pos().x(), event->pos().y());

            break;
        }
    }
}

void Canvas::mouseMoveEvent(QMouseEvent *event)
{
    switch(event->buttons())
    {
        case Qt::MouseButton::LeftButton:
        {
            Particle p;
            p.pos = glm::vec2(event->pos().x(), height()-event->pos().y());
            solver->addParticle(p);
            break;
        }
        case Qt::MouseButton::RightButton:
        {
            line[1] = glm::vec2(event->pos().x(), event->pos().y());
            break;
        }
    }
    update();

}

void Canvas::mouseReleaseEvent(QMouseEvent *event)
{
    switch (event->button())
    {
        case Qt::MouseButton::LeftButton:
        {
            break;
        }
        case Qt::MouseButton::RightButton:
        {
            Boundary boundary;
            boundary.line[0] = line[0];
            boundary.line[0].y = height() - boundary.line[0].y;
            boundary.line[1] = line[1];
            boundary.line[1].y = height() - boundary.line[1].y;
            if(line[0].x>line[1].x)
            {
                glm::vec2 temp = line[0];
                line[0] = line[1];
                line[1] = temp;
            }
            glm::vec2 line_vec = (line[1]-line[0]);
            glm::vec2 normal = glm::vec2(line_vec.y,line_vec.x);
            boundary.normal = glm::normalize(normal);
            solver->addBoundary(boundary);

            line[0] = glm::vec2(0.0f);
            line[1] = glm::vec2(0.0f);

            break;
        }
    }
}
