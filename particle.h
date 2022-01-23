#ifndef PARTICLE_H
#define PARTICLE_H
#include <glm/glm.hpp>

struct Particle
{
    unsigned int bucked_id;
    glm::vec2 pos;
    glm::vec2 vel;
    glm::vec2 proj_pos;
    glm::vec2 delta;
    float lambda;
    float pressure;
};

#endif // PARTICLE_H
