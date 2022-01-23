#ifndef RADIXSORT_H
#define RADIXSORT_H
#include "abstractspatialstruct.h"
#include <atomic>

class RadixSort : public AbstractSpatialStruct
{
public:
    RadixSort(std::vector<Particle>& particles);

    void rebuild(float width, float height, float radius) override;
    std::vector<unsigned int> findNeighbors(unsigned int index, float r) override;
private:
    unsigned int discrete_width;
    unsigned int discrete_height;
    std::vector<unsigned int> histogram;
};

#endif // LINEARSEARCH_H
