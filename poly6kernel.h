#ifndef POLY6KERNEL_H
#define POLY6KERNEL_H
#include<abstractkernel.h>

class Poly6Kernel : public AbstractKernel
{
public:
    Poly6Kernel();
    float evaluate(float r, float h) override;
    glm::vec2 gradient(float r, float h, const glm::vec2 n) override;
    float laplace(float r, float h) override;
};

#endif // POLY6KERNEL_H
