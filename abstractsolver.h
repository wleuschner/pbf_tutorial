#ifndef ABSTRACTSOLVER_H
#define ABSTRACTSOLVER_H
#include "particle.h"
#include "boundary.h"
#include "abstractspatialstruct.h"
#include <vector>

class AbstractSolver
{
public:
    enum SpatialStructType
    {
        LINEAR_SEARCH,
        RADIX_SORT
    };
    AbstractSolver(unsigned int width, unsigned int height);
    virtual ~AbstractSolver();

    virtual void solve() = 0;

    void addBoundary(Boundary& boundary);
    void addParticle(Particle& particle);

    void changeDomain(float right, float top, float left, float bottom);

    void setIterations(float niter);
    void setTimeStep(float step);
    void setSearchRadius(float radius);
    void setRestingDensity(float restingDensity);
    void setArtificialDensity(float artificialDensity);
    void setViscosity(float visc);
    void setSpatialStruct(SpatialStructType type);

    void setParticles(std::vector<Particle> &particles);
    void setBoundaries(std::vector<Boundary> &particles);


    const std::vector<Boundary>& getBoundaries();
    const std::vector<Particle>& getParticles();
protected:
    const glm::vec2 gravity = glm::vec2(0.0,-9.81);

    AbstractSpatialStruct* spatial_struct;
    float domain_width;
    float domain_height;

    std::vector<Particle> particles;
    std::vector<Boundary> boundaries;

    int iterations;
    float timestep;
    float artificial_density;
    float search_radius;
    float resting_density;
    float viscosity;
};

#endif
