#include "radixsort.h"

RadixSort::RadixSort(std::vector<Particle>& particles)
  : AbstractSpatialStruct(particles)
{
}

RadixSort::~RadixSort()
{

}

void RadixSort::rebuild(float width, float height, float radius)
{
    discrete_width = ceil(width/radius);
    discrete_height = ceil(height/radius);
    unsigned int num_buckets = discrete_width * discrete_height;

    std::vector<std::atomic<unsigned int>> buckets(num_buckets);
    #pragma omp parallel for
    for(unsigned int i=0;i<num_buckets;i++)
    {
        buckets[i] = 0;
    }

    //Compute bucket sizes
    #pragma omp parallel for
    for(unsigned int i=0;i<particles.size();i++)
    {
        Particle& particle = particles[i];
        int bucked_id = (particle.pos.x/discrete_width) + (particle.pos.y/discrete_height) * discrete_width;
        particle.bucked_id = bucked_id;
        buckets[bucked_id].fetch_add(1);
    }

    //Compute histogram
    histogram.resize(num_buckets);
    unsigned int offset = 0;
    for(unsigned int i=0;i<num_buckets;i++)
    {
        offset += buckets[i];
        histogram[i] = offset;
    }

    #pragma omp parallel for
    for(unsigned int i=0;i<num_buckets;i++)
    {
        buckets[i] = 0;
    }

    if(particlesSorted.size()!=particles.size())
    {
        particlesSorted.resize(particles.size());
    }
    #pragma omp parallel for
    for(unsigned int i=0;i<particles.size();i++)
    {
        Particle& particle = particles[i];
        unsigned int bucket_id = particle.bucked_id;
        unsigned int histogram_offset = bucket_id>0 ? histogram[bucket_id-1] : 0;
        unsigned int offset = histogram_offset + buckets[bucket_id].fetch_add(1);
        particlesSorted[offset] = particle;
    }
    particles.swap(particlesSorted);
}

std::vector<unsigned int> RadixSort::findNeighbors(unsigned int index, float r)
{
    const Particle particle = particles[index];
    std::vector<unsigned int> neighbors;
    float r_squared = r*r;
    for(int y=-1;y<=1;y++)
    {
        int y_offset = particle.bucked_id/discrete_width + y;
        if(y_offset<0 || y_offset>=discrete_height) continue;
        for(int x=-1;x<=1;x++)
        {
            int x_offset = particle.bucked_id%discrete_width + x;
            if(x_offset<0 || x_offset>=discrete_width) continue;

            unsigned int array_offset = x_offset+y_offset*discrete_width;
            unsigned int range_begin = x_offset>0 ? histogram[array_offset - 1] : 0;
            unsigned int range_end = histogram[array_offset];
            for(unsigned int i=range_begin;i<range_end;i++)
            {
                if(i==index) continue;
                Particle& p = particles[i];
                glm::vec2 distVec = p.pos-particle.proj_pos;
                float d = glm::dot(distVec,distVec);
                if(d<r_squared && p.pos!=particle.proj_pos)
                {
                    neighbors.push_back(i);
                }
            }
        }
    }
    return neighbors;
}
