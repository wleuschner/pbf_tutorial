#ifndef PARTICLE_H
#define PARTICLE_H
#include <glm/glm.hpp>

struct Particle
{
    glm::vec2 pos;
    glm::vec2 vel;
    glm::vec2 proj_pos;
    glm::vec2 delta;
    float lambda;
    float pressure;
    unsigned int bucked_id;
    char pad[2];
};

#endif // PARTICLE_H
