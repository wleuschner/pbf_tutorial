#ifndef SPIKYKERNEL_H
#define SPIKYKERNEL_H
#include<abstractkernel.h>

class SpikyKernel : public AbstractKernel
{
public:
    SpikyKernel();
    float evaluate(float r, float h) override;
    glm::vec2 gradient(float r, float h, const glm::vec2 n) override;
    float laplace(float r, float h) override;
};

#endif // POLY6KERNEL_H
