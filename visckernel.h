#ifndef VISCKERNEL_H
#define VISC6KERNEL_H
#include<abstractkernel.h>

class ViscKernel : public AbstractKernel
{
public:
    ViscKernel();
    float evaluate(float r, float h) override;
    glm::vec2 gradient(float r, float h, const glm::vec2 n) override;
    float laplace(float r, float h) override;
};

#endif // POLY6KERNEL_H
