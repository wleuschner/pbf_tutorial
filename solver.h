#ifndef SOLVER_H
#define SOLVER_H
#include <vector>
#include "particle.h"
#include "boundary.h"
#include "abstractspatialstruct.h"
#include "abstractkernel.h"
#include <mutex>

class Solver
{
public:
    enum SpatialStructType
    {
        LINEAR_SEARCH,
        RADIX_SORT
    };
    Solver(unsigned int width, unsigned int height);
    ~Solver();

    void solve();

    void changeDomain(float right, float top, float left, float bottom);
    void addBoundary(Boundary& boundary);
    void addParticle(Particle& particle);

    void setIterations(float niter);
    void setTimeStep(float step);
    void setSearchRadius(float radius);
    void setRestingDensity(float restingDensity);
    void setArtificialDensity(float artificialDensity);
    void setViscosity(float visc);
    void setSpatialStruct(SpatialStructType type);

    const std::vector<Boundary>& getBoundaries();
    const std::vector<Particle>& getParticles();
private:
    float computeDensity(const Particle& particle, const std::vector<unsigned int>& neighbors);
    float computePressureGradientSum(const Particle& particle, const std::vector<unsigned int> &neighbors);

    float domain_width = 0;
    float domain_height = 0;

    const glm::vec2 gravity = glm::vec2(0.0,-9.81);

    Boundary domain_boundaries[4];

    AbstractSpatialStruct* spatial_struct;
    std::vector<Particle> particles;
    std::vector<Boundary> boundaries;

    AbstractKernel* densityKernel;
    AbstractKernel* gradientKernel;
    AbstractKernel* viscKernel;

    int iterations;
    float timestep;
    float artificial_density;
    float search_radius;
    float resting_density;
    float viscosity;

    std::mutex mut;
};

#endif // SOLVER_H
