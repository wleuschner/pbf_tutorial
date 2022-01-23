#include "linearsearch.h"

LinearSearch::LinearSearch(std::vector<Particle>& particles)
    : AbstractSpatialStruct(particles)
{

}

void LinearSearch::rebuild(float width, float height, float radius)
{
}

std::vector<unsigned int> LinearSearch::findNeighbors(unsigned int index, float r)
{
    const Particle particle = particles[index];
    std::vector<unsigned int> neighbors;
    for(unsigned int i=0;i<particles.size();i++)
    {
        if(i==index) continue;
        Particle& p = particles[i];
        float d = glm::length(p.pos-particle.proj_pos);
        if(d<r && p.pos!=particle.proj_pos)
        {
            neighbors.push_back(i);
        }
    }
    return neighbors;
}
