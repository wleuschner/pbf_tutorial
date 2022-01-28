#ifndef LINEARSEARCH_H
#define LINEARSEARCH_H
#include "abstractspatialstruct.h"


class LinearSearch : public AbstractSpatialStruct
{
public:
    LinearSearch(std::vector<Particle>& particles);
    ~LinearSearch();

    void rebuild(float width, float height, float radius) override;
    std::vector<unsigned int> findNeighbors(unsigned int index, float r) override;
};

#endif // LINEARSEARCH_H
