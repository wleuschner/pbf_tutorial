#ifndef ABSTRACTSPATIALSTRUCT_H
#define ABSTRACTSPATIALSTRUCT_H
#include <vector>
#include "particle.h"

class AbstractSpatialStruct
{
public:
    AbstractSpatialStruct(std::vector<Particle>& particles);
    virtual ~AbstractSpatialStruct();
    virtual void rebuild(float width, float height, float radius) = 0;
    virtual std::vector<unsigned int> findNeighbors(unsigned int index, float r) = 0;
protected:
    std::vector<Particle>& particles;
};

#endif // ABSTRACTSPATIALSTRUCT_H
