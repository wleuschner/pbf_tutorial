#include "visckernel.h"
#include <cmath>

ViscKernel::ViscKernel()
{

}

float ViscKernel::evaluate(float r, float h)
{
    return 0.0f;
}

glm::vec2 ViscKernel::gradient(float r, float h, const glm::vec2 n)
{
    return glm::vec2(0.0f,0.0f);
}

float ViscKernel::laplace(float r, float h)
{

    if(r==0.0f)
    {
        return 0.0f;
    }
    return 45.0f/(M_PI*std::pow(h,6))*(h-r);
}
