#ifndef ABSTRACTKERNEL_H
#define ABSTRACTKERNEL_H
#include <glm/vec2.hpp>

class AbstractKernel
{
public:
    AbstractKernel();
    virtual ~AbstractKernel();

    virtual float evaluate(float r, float h) = 0;
    virtual glm::vec2 gradient(float r, float h,const glm::vec2 n) = 0;
    virtual float laplace(float r, float h) = 0;
};

#endif // ABSTRACTKERNEL_H
