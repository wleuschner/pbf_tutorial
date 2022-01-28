#include "abstractsolver.h"
#include "linearsearch.h"
#include "radixsort.h"

AbstractSolver::AbstractSolver(unsigned int width, unsigned int height)
{
    domain_width = width;
    domain_height = height;
    boundaries.resize(4);

    iterations = 10;
    timestep = 0.1f;
    search_radius = 10.0f;
    resting_density = 10.0f;
    artificial_density = 0.0001f;
    viscosity = 0.0f;

    spatial_struct = new LinearSearch(particles);
}

AbstractSolver::~AbstractSolver()
{
    delete spatial_struct;
}

void AbstractSolver::addBoundary(Boundary& boundary)
{
    boundaries.push_back(boundary);
}

void AbstractSolver::addParticle(Particle& particle)
{
    particle.vel = glm::vec2(0.0f,0.0f);
    particles.push_back(particle);
}

void AbstractSolver::changeDomain(float right, float top, float left, float bottom)
{
    domain_width = (-right+left);
    domain_height = (-bottom+top);

    boundaries[0].line[0] = glm::vec2(left,bottom);
    boundaries[0].line[1] = glm::vec2(left,top);
    boundaries[0].normal = glm::vec2(-1.0,0.0);

    boundaries[1].line[0] = glm::vec2(left,top);
    boundaries[1].line[1] = glm::vec2(right,top);
    boundaries[1].normal = glm::vec2(0.0,-1.0);


    boundaries[2].line[0] = glm::vec2(right,top);
    boundaries[2].line[1] = glm::vec2(right,bottom);
    boundaries[2].normal = glm::vec2(-1.0,0.0);


    boundaries[3].line[0] = glm::vec2(right,bottom);
    boundaries[3].line[1] = glm::vec2(left,bottom);
    boundaries[3].normal = glm::vec2(0.0,-1.0);

}

void AbstractSolver::setIterations(float niter)
{
    iterations = niter;
}

void AbstractSolver::setTimeStep(float step)
{
    timestep = step;
}

void AbstractSolver::setSearchRadius(float radius)
{
    search_radius = radius;
}

void AbstractSolver::setRestingDensity(float restingDensity)
{
    resting_density = restingDensity;
}

void AbstractSolver::setArtificialDensity(float artificialDensity)
{
    artificial_density = artificialDensity;
}

void AbstractSolver::setViscosity(float visc)
{
    viscosity = visc;
}

void AbstractSolver::setSpatialStruct(SpatialStructType type)
{
    delete spatial_struct;
    if(type == LINEAR_SEARCH)
    {
        spatial_struct = new LinearSearch(particles);
    }
    else if(type == RADIX_SORT)
    {
        spatial_struct = new RadixSort(particles);
    }
}

void AbstractSolver::setParticles(std::vector<Particle> &parts)
{
    particles.clear();
    particles.insert(particles.begin(), parts.begin(), parts.end());
}

void AbstractSolver::setBoundaries(std::vector<Boundary> &bounds)
{
    boundaries.clear();
    boundaries.insert(boundaries.begin(), bounds.begin(), bounds.end());
}

const std::vector<Boundary>& AbstractSolver::getBoundaries()
{
    return boundaries;
}

const std::vector<Particle>& AbstractSolver::getParticles()
{
    return particles;
}
