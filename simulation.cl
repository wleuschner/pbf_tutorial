#pragma OPENCL EXTENSION cl_khr_global_int32_base_atomics : enable
typedef struct
{
    float2 pos;
    float2 vel;
    float2 proj_pos;
    float2 delta;
    float lambda;
    float pressure;
    uint id;
} Particle;

typedef struct
{
    float2 begin;
    float2 end;
    float2 normal;
} Boundary;

float poly6_kernel(float r, float h)
{
    return (315/(64*M_PI_F*pown(h,9)))*pown(h*h-r*r,3);
}

float2 spiky_kernel_gradient(float r, float h, float2 n)
{
    return ((float)(-45.0f/(M_PI_F*pown(h,6))*pown(h-r,2)))*n;
}

float visc_kernel(float r, float h)
{
    return 45.0f/(M_PI_F*pown(h,6))*(h-r);
}


float3 intersect_line_line(const float2 a, const float2 b,const float2 c, const float2 d)
{
    float t = ((a.x-c.x)*(c.y-d.y)-(a.y-c.y)*(c.x-d.x))/((a.x-b.x)*(c.y-d.y)-(a.y-b.y)*(c.x-d.x));
    return (float3)(a+t*(b-a),t);
}

void __kernel bucket_count(global Particle* particles, global int* counters, uint2 dims)
{
    int index = get_global_id(0);
    Particle* particle = &particles[index];
    particle->id = particle->pos.x/dims.x + (particle->pos.y/dims.y)*dims.x;
    atomic_add(&counters[particle->id],1);
}

void __kernel histogram(global int* counters, global int* histogram, global int* collect)
{
    int index = get_global_id(0);
    histogram[index] = work_group_scan_inclusive_add(counters[index]);
    if(get_local_id(0)==get_local_size(0)-1)
    {
        collect[get_group_id(0)] = histogram[index];
    }
}
void __kernel histogram2(global int* counters, global int* histogram)
{
    int index = get_global_id(0);
    histogram[index] = work_group_scan_inclusive_add(counters[index]);
}
void __kernel histogram3(global int* counters, global int* histogram, global int* temp)
{
    int index = get_local_id(0) + (get_group_id(0)+1)*get_local_size(0);
    histogram[index] += counters[get_group_id(0)];
}

void __kernel sort_particles(global Particle* particles_front, global Particle* particles_back, global uint* offsets, global uint* temp)
{
    uint index = get_global_id(0);
    Particle* particle = &particles_front[index];
    uint counter_offset = particle->id>0 ? offsets[particle->id-1] : 0;
    uint offset = counter_offset + atomic_add(&temp[particle->id],1);
    particles_back[offset] = *particle;
}

void __kernel add_external_forces(global Particle* particles, float timestep, float2 gravity)
{
    Particle* particle = &particles[get_global_id(0)];
    particle->vel += timestep*gravity;
    particle->proj_pos = particle->pos + timestep * particle->vel;
}

void __kernel compute_lambda(global Particle* particles, uint2 dims, global int* histogram, float radius, float resting_density, float artificial_density)
{
        uint index = get_global_id(0);
        Particle* particle = &particles[index];
        float density = 0.0f;
        float pressureSum = 0.0f;
        float2 out_pressure = (float2)(0.0f,0.0f);
        float r_squared = radius*radius;
        for(int y=-1;y<=1;y++)
        {
            int y_offset = particle->id/dims.x + y;
            if(y_offset<0 || y_offset>=dims.y) continue;
            for(int x=-1;x<=1;x++)
            {
                int x_offset = particle->id%dims.x + x;
                if(x_offset<0 || x_offset>=dims.x) continue;

                uint array_offset = x_offset+y_offset*dims.x;
                uint range_begin = x_offset>0 ? histogram[array_offset - 1] : 0;
                uint range_end = histogram[array_offset];
                for(uint i=range_begin;i<range_end;i++)
                {
                    if(i==index) continue;
                    Particle* p = &particles[i];
                    float2 distVec = p->pos-particle->proj_pos;
                    float d = dot(distVec,distVec);
                    if(d<r_squared)
                    {
                        density += poly6_kernel(length(distVec),radius);
                        float2 in_pressure = -(1.0f/resting_density) * spiky_kernel_gradient(length(distVec), radius, normalize(distVec));
                        pressureSum += dot(in_pressure,in_pressure);
                        out_pressure += spiky_kernel_gradient(length(distVec), radius, normalize(distVec));
                    }
                }
            }
        }
        out_pressure = (1.0f/resting_density) * out_pressure;
        pressureSum += dot(out_pressure,out_pressure);
        float densityConstraint = density/resting_density;
        particle->pressure = density/resting_density;
        particle->lambda = -densityConstraint/(pressureSum + artificial_density);
}

void __kernel update_delta(global uint* histogram, uint2 dims, float search_radius, float resting_density, global Particle* particles)
{
    int index = get_global_id(0);
    Particle* particle = &particles[index];
    float2 delta = (float2)(0.0f,0.0f);
    float2 particle_proj_pos = particle->proj_pos;
    float2 particle_lambda = particle->lambda;
    float r_squared = search_radius*search_radius;
    for(int y=-1;y<=1;y++)
    {
        int y_offset = particle->id/dims.x + y;
        if(y_offset<0 || y_offset>=dims.y) continue;
        for(int x=-1;x<=1;x++)
        {
            int x_offset = particle->id%dims.x + x;
            if(x_offset<0 || x_offset>=dims.x) continue;
            uint array_offset = x_offset+y_offset*dims.x;
            uint range_begin = x_offset>0 ? histogram[array_offset - 1] : 0;
            uint range_end = histogram[array_offset];
            for(uint i=range_begin;i<range_end;i++)
            {
                if(i==index) continue;
                Particle* p = &particles[i];
                float2 distVec = -p->pos+particle_proj_pos;
                float d = dot(distVec,distVec);
                if(d<r_squared)
                {
                    float sCorr = -0.1f * pown((poly6_kernel(length(distVec),search_radius)/poly6_kernel(0.3*search_radius,search_radius)),4);
                    delta += (1.0f/resting_density) * (particle_lambda + p->lambda + sCorr) * spiky_kernel_gradient(length(distVec),search_radius,normalize(distVec));
                }
            }
        }
    }
    particle->delta = delta;
}

void __kernel check_boundary_collision(global Particle* particles, global Boundary* boundaries, uint num_boundaries)
{
    int index = get_global_id(0);
    Particle* particle = &particles[index];
    for(int i=0;i<num_boundaries;i++)
    {
        Boundary* bound = &boundaries[i];
        float d1 = dot(bound->normal,bound->begin-particle->pos);
        float d2 = dot(bound->normal,bound->begin-((particle->proj_pos + particle->delta)));
        if((sign(d1)!=sign(d2)))
        {
            float si = sign(d2);
            float3 pos = intersect_line_line(bound->begin,
                                            bound->end,
                                            particle->pos,
                                            particle->proj_pos + particle->delta);
            if(pos.z<0.0 || pos.z>1.0) continue;
            particle->delta = -particle->proj_pos + ((float2)(pos.x,pos.y) + 0.1f * si*bound->normal);
        }
    }
}
void __kernel update_projected_positions(global Particle* particles)
{
    int index = get_global_id(0);
    Particle* particle = &particles[index];
    particle->proj_pos = particle->proj_pos + particle->delta;
}

void __kernel update_particles(global Particle* particles, float timestep, uint2 dims, global uint* histogram, float visc, float radius)
{
    int index = get_global_id(0);
    Particle* particle = &particles[index];
    particle->vel = (1.0f/timestep) * (particle->proj_pos-particle->pos);
    barrier(CLK_GLOBAL_MEM_FENCE);
    float2 viscAcc = (float2)(0.0f,0.0f);
    float r_squared = radius * radius;
    for(int y=-1;y<=1;y++)
    {
        int y_offset = particle->id/dims.x + y;
        if(y_offset<0 || y_offset>=dims.y) continue;
        for(int x=-1;x<=1;x++)
        {
            int x_offset = particle->id%dims.x + x;
            if(x_offset<0 || x_offset>=dims.x) continue;

            uint array_offset = x_offset+y_offset*dims.x;
            uint range_begin = x_offset>0 ? histogram[array_offset - 1] : 0;
            uint range_end = histogram[array_offset];
            for(uint i=range_begin;i<range_end;i++)
            {
                if(i==index) continue;
                Particle* p = &particles[i];
                float2 distVec = p->pos-particle->proj_pos;
                float d = dot(distVec,distVec);
                if(d<r_squared)
                {
                    float2 viscVel = p->vel - particle->vel;
                    viscAcc += viscVel*visc_kernel(distance(particle->pos,p->pos),radius);
                }
            }
        }
    }
    barrier(CLK_LOCAL_MEM_FENCE);
    particle->vel += visc*viscAcc;
    particle->pos = particle->proj_pos;
}
