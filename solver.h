#ifndef SOLVER_H
#define SOLVER_H
#include <vector>
#include "particle.h"
#include "boundary.h"
#include "abstractsolver.h"
#include "abstractspatialstruct.h"
#include "abstractkernel.h"

class Solver : public AbstractSolver
{
public:
    Solver(unsigned int width, unsigned int height);
    ~Solver();

    void solve() override;
private:
    float computeDensity(const Particle& particle, const std::vector<unsigned int>& neighbors);
    float computePressureGradientSum(const Particle& particle, const std::vector<unsigned int> &neighbors);

    std::vector<glm::vec2> viscAcc;
    AbstractKernel* densityKernel;
    AbstractKernel* gradientKernel;
    AbstractKernel* viscKernel;
};

#endif // SOLVER_H
