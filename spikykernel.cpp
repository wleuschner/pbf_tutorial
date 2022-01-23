#include "spikykernel.h"
#include <cmath>

SpikyKernel::SpikyKernel()
{

}

float SpikyKernel::evaluate(float r, float h)
{
    return (15/(M_PI*std::pow(h,6)))*std::pow(h-r,3);
}

glm::vec2 SpikyKernel::gradient(float r, float h,const glm::vec2 n)
{
    if(r==0.0f)
    {
        return glm::vec2(0.0f);
    }
    return ((float)(-45.0f/(M_PI*std::pow(h,6))*std::pow(h-r,2)))*n;
}

float SpikyKernel::laplace(float r, float h)
{

}
