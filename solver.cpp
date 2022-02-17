#include "solver.h"
#include "linearsearch.h"
#include "radixsort.h"
#include "poly6kernel.h"
#include "spikykernel.h"
#include "visckernel.h"
#include <cstdio>
#include <glm/vec2.hpp>

float sign(float a)
{
    return a>0.0f ? 1.0f : -1.0f;
}

glm::vec3 intersect_line_line(const glm::vec2& a, const glm::vec2& b,
                              const glm::vec2& c, const glm::vec2& d)
{
    float t = ((a.x-c.x)*(c.y-d.y)-(a.y-c.y)*(c.x-d.x))/((a.x-b.x)*(c.y-d.y)-(a.y-b.y)*(c.x-d.x));
    return glm::vec3(a+t*(b-a),t);
}

Solver::Solver(unsigned int width, unsigned int height) : AbstractSolver(width,height)
{
    densityKernel = new Poly6Kernel();
    gradientKernel = new SpikyKernel();
    viscKernel = new ViscKernel();
}

Solver::~Solver()
{
    delete densityKernel;
    delete gradientKernel;
    delete viscKernel;
}

void Solver::solve()
{
    spatial_struct->rebuild(domain_width,domain_height,search_radius);

    #pragma omp parallel for
    for(unsigned int i=0;i<particles.size();i++)
    {
        Particle& p = particles[i];
        p.vel += timestep*gravity;
        p.proj_pos = p.pos+timestep*p.vel;
    }

    for(unsigned int i=0;i<iterations;i++)
    {
        //Calculate Lambda
        #pragma omp parallel for
        for(unsigned int p=0;p<particles.size();p++)
        {
            Particle& particle = particles[p];
            std::vector<unsigned int> neighbors = spatial_struct->findNeighbors(p, search_radius);
            float density = computeDensity(particle, neighbors);
            float densityConstraint = density/resting_density;
            float gradSum = computePressureGradientSum(particle, neighbors);
            particle.pressure = density/resting_density;
            particle.lambda = -densityConstraint/(gradSum + artificial_density);
        }
        //Collision response
        #pragma omp parallel for
        for(unsigned int p=0;p<particles.size();p++)
        {
            Particle& particle = particles[p];
            std::vector<unsigned int> neighbors = spatial_struct->findNeighbors(p, search_radius);
            particle.delta = glm::vec2(0.0f,0.0f);
            for(unsigned int i=0;i<neighbors.size();i++)
            {
                const Particle& neighbor = particles[neighbors[i]];
                glm::vec2 distVec = particle.proj_pos - neighbor.pos;
                float sCorr = -0.1f * std::pow((densityKernel->evaluate(glm::length(distVec),search_radius)/densityKernel->evaluate(0.3*search_radius,search_radius)),4.0f);
                particle.delta += (1.0f/resting_density) * (particle.lambda + neighbor.lambda + sCorr) * gradientKernel->gradient(glm::distance(particle.proj_pos,neighbor.pos),search_radius,glm::normalize(particle.proj_pos - neighbor.pos));
            }

            for(unsigned int b=0;b<boundaries.size();b++)
            {
                Boundary& bound = boundaries[b];
                float d1 = glm::dot(bound.normal,(bound.line[0]-particle.pos));
                float d2 = glm::dot(bound.normal,(bound.line[0]-(particle.proj_pos + particle.delta)));

                //Projected particle is on other side than curren position
                //Thus check for collision response
                if((sign(d1)!=sign(d2)) || (sign(d1)==0.0f && sign(d2)!=0.0f))
                {

                    float si = sign(d2);
                    glm::vec3 pos = intersect_line_line(bound.line[0],
                                                        bound.line[1],
                                                        particle.pos,
                                                        particle.proj_pos + particle.delta);
                    //If particle is not on line segment, ignore response
                    if(pos.z<0.0 || pos.z>1.0) continue;

                    particle.delta = -particle.proj_pos + (glm::vec2(pos) + 0.1f * si*bound.normal);
                }
            }
        }

        //Update projected particle position for next solver iteration
        #pragma omp parallel for
        for(unsigned int p=0;p<particles.size();p++)
        {
            Particle& particle = particles[p];
            particle.proj_pos = particle.proj_pos + particle.delta;
        }
    }

    //Update particle positions and velocities based on projected positions
    std::vector<glm::vec2> viscAcc(particles.size());
    #pragma omp parallel for
    for(unsigned int i=0;i<particles.size();i++)
    {
        std::vector<unsigned int> neighbors = spatial_struct->findNeighbors(i,search_radius);

        Particle& p = particles[i];
        p.vel = (1.0f/timestep)*(p.proj_pos-p.pos);
        p.pos = p.proj_pos;
    }

    //Compute viscosity influence
    for(unsigned int i=0;i<particles.size();i++)
    {
        std::vector<unsigned int> neighbors = spatial_struct->findNeighbors(i,search_radius);
        Particle& p = particles[i];

        viscAcc[i] = glm::vec2(0.0f,0.0f);
        for(unsigned int n=0;n<neighbors.size();n++)
        {
            const Particle& neighbor = particles[neighbors[n]];
            glm::vec2 viscVel = p.vel - neighbor.vel;
            viscAcc[i] += viscVel*viscKernel->laplace(glm::distance(p.pos,neighbor.pos),search_radius);
        }
    }

    //Update Final Velocity
    for(unsigned int i=0;i<particles.size();i++)
    {
        std::vector<unsigned int> neighbors = spatial_struct->findNeighbors(i,search_radius);
        Particle& p = particles[i];
        p.vel += viscosity*viscAcc[i];
    }
}

float Solver::computeDensity(const Particle& particle, const std::vector<unsigned int> &neighbors)
{
    float density = 0.0f;
    for(unsigned int p=0;p<neighbors.size();p++)
    {
        const Particle& neighbor = particles[neighbors[p]];
        density += densityKernel->evaluate(glm::distance(particle.proj_pos,neighbor.pos),search_radius);
    }
    return density;
}

float Solver::computePressureGradientSum(const Particle& particle, const std::vector<unsigned int> &neighbors)
{
    float pressureSum = 0.0f;

    glm::vec2 out_pressure = glm::vec2(0.0f, 0.0f);

    for(int i=0;i<neighbors.size();i++)
    {
        const Particle& neighbor = particles[neighbors[i]];
        glm::vec2 distVec = particle.proj_pos - neighbor.pos;
        out_pressure += gradientKernel->gradient(glm::length(distVec), search_radius, glm::normalize(distVec));
    }
    out_pressure = (1.0f/resting_density) * out_pressure;
    pressureSum = glm::dot(out_pressure,out_pressure);

    for(int i=0;i<neighbors.size();i++)
    {
        const Particle& neighbor = particles[neighbors[i]];
        glm::vec2 distVec = particle.proj_pos - neighbor.pos;
        glm::vec2 in_pressure = -(1.0f/resting_density) * gradientKernel->gradient(glm::length(distVec), search_radius, glm::normalize(distVec));
        pressureSum += glm::dot(in_pressure,in_pressure);
    }
    return pressureSum;
}
