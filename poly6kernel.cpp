#include "poly6kernel.h"
#include <cmath>

Poly6Kernel::Poly6Kernel()
{

}

Poly6Kernel::~Poly6Kernel()
{

}

float Poly6Kernel::evaluate(float r, float h)
{
    if(r==0.0f)
    {
        return 0.0f;
    }
    return (315/(64*M_PI*std::pow(h,9)))*std::pow(h*h-r*r,3);
}

glm::vec2 Poly6Kernel::gradient(float r, float h, const glm::vec2 n)
{

}

float Poly6Kernel::laplace(float r, float h)
{

}
