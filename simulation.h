#ifndef SIMULATION_H
#define SIMULATION_H

static const char* simulation_source =
"#pragma OPENCL EXTENSION cl_khr_global_int32_base_atomics : enable\n"
"typedef struct\n"
"{\n"
"    float2 pos;\n"
"    float2 vel;\n"
"    float2 proj_pos;\n"
"    float2 delta;\n"
"    float lambda;\n"
"    float pressure;\n"
"    uint id;\n"
"} Particle;\n"
"\n"
"typedef struct\n"
"{\n"
"    float2 begin;\n"
"    float2 end;\n"
"    float2 normal;\n"
"} Boundary;\n"
"\n"
"float poly6_kernel(float r, float h)\n"
"{\n"
"    return (315/(64*M_PI_F*pown(h,9)))*pown(h*h-r*r,3);\n"
"}\n"
"\n"
"float2 spiky_kernel_gradient(float r, float h, float2 n)\n"
"{\n"
"    return ((float)(-45.0f/(M_PI_F*pown(h,6))*pown(h-r,2)))*n;\n"
"}\n"
"\n"
"float visc_kernel(float r, float h)\n"
"{\n"
"    return 45.0f/(M_PI_F*pown(h,6))*(h-r);\n"
"}\n"
"\n"
"\n"
"float3 intersect_line_line(const float2 a, const float2 b,const float2 c, const float2 d)\n"
"{\n"
"    float t = ((a.x-c.x)*(c.y-d.y)-(a.y-c.y)*(c.x-d.x))/((a.x-b.x)*(c.y-d.y)-(a.y-b.y)*(c.x-d.x));\n"
"    return (float3)(a+t*(b-a),t);\n"
"}\n"
"\n"
"void __kernel bucket_count(global Particle* particles, global int* counters, uint2 dims)\n"
"{\n"
"    int index = get_global_id(0);\n"
"    Particle* particle = &particles[index];\n"
"    particle->id = particle->pos.x/dims.x + (particle->pos.y/dims.y)*dims.x;\n"
"    atomic_add(&counters[particle->id],1);\n"
"}\n"
"\n"
"void __kernel histogram(global int* counters, global int* histogram, global int* collect)\n"
"{\n"
"    int index = get_global_id(0);"
"    histogram[index] = work_group_scan_inclusive_add(counters[index]);\n"
"    if(get_local_id(0)==get_local_size(0)-1)\n"
"    {\n"
"        collect[get_group_id(0)] = histogram[index];"
"    }\n"
"}\n"
"void __kernel histogram2(global int* counters, global int* histogram)\n"
"{\n"
"    int index = get_global_id(0);"
"    histogram[index] = work_group_scan_inclusive_add(counters[index]);\n"
"}\n"
"void __kernel histogram3(global int* counters, global int* histogram, global int* temp)\n"
"{\n"
"    int index = get_local_id(0) + (get_group_id(0)+1)*get_local_size(0);"
"    histogram[index] += counters[get_group_id(0)];\n"
"}\n"
"\n"
"void __kernel sort_particles(global Particle* particles_front, global Particle* particles_back, global uint* offsets, global uint* temp)\n"
"{\n"
"    uint index = get_global_id(0);\n"
"    Particle* particle = &particles_front[index];\n"
"    uint counter_offset = particle->id>0 ? offsets[particle->id-1] : 0;\n"
"    uint offset = counter_offset + atomic_add(&temp[particle->id],1);\n"
"    particles_back[offset] = *particle;"
"}\n"
"\n"
"void __kernel add_external_forces(global Particle* particles, float timestep, float2 gravity)\n"
"{\n"
"    Particle* particle = &particles[get_global_id(0)];\n"
"    particle->vel += timestep*gravity;\n"
"    particle->proj_pos = particle->pos + timestep * particle->vel;\n"
"}\n"
"\n"
"void __kernel compute_lambda(global Particle* particles, uint2 dims, global int* histogram, float radius, float resting_density, float artificial_density)\n"
"{\n"
"        uint index = get_global_id(0);\n"
"        Particle* particle = &particles[index];\n"
"        float density = 0.0f;\n"
"        float pressureSum = 0.0f;\n"
"        float2 out_pressure = (float2)(0.0f,0.0f);\n"
"        float r_squared = radius*radius;\n"
"        for(int y=-1;y<=1;y++)\n"
"        {\n"
"            int y_offset = particle->id/dims.x + y;\n"
"            if(y_offset<0 || y_offset>=dims.y) continue;\n"
"            for(int x=-1;x<=1;x++)\n"
"            {\n"
"                int x_offset = particle->id%dims.x + x;\n"
"                if(x_offset<0 || x_offset>=dims.x) continue;\n"
"\n"
"                uint array_offset = x_offset+y_offset*dims.x;\n"
"                uint range_begin = x_offset>0 ? histogram[array_offset - 1] : 0;\n"
"                uint range_end = histogram[array_offset];\n"
"                for(uint i=range_begin;i<range_end;i++)\n"
"                {\n"
"                    if(i==index) continue;\n"
"                    Particle* p = &particles[i];\n"
"                    float2 distVec = p->pos-particle->proj_pos;\n"
"                    float d = dot(distVec,distVec);\n"
"                    if(d<r_squared)\n"
"                    {\n"
"                        density += poly6_kernel(length(distVec),radius);\n"
"                        float2 in_pressure = -(1.0f/resting_density) * spiky_kernel_gradient(length(distVec), radius, normalize(distVec));\n"
"                        pressureSum += dot(in_pressure,in_pressure);\n"
"                        out_pressure += spiky_kernel_gradient(length(distVec), radius, normalize(distVec));\n"
"                    }\n"
"                }\n"
"            }\n"
"        }\n"
"        out_pressure = (1.0f/resting_density) * out_pressure;\n"
"        pressureSum += dot(out_pressure,out_pressure);\n"
"        float densityConstraint = density/resting_density;\n"
"        particle->pressure = density/resting_density;\n"
"        particle->lambda = -densityConstraint/(pressureSum + artificial_density);\n"
"}\n"
"\n"
"void __kernel update_delta(global uint* histogram, uint2 dims, float search_radius, float resting_density, global Particle* particles)\n"
"{\n"
"    int index = get_global_id(0);\n"
"    Particle* particle = &particles[index];\n"
"    float2 delta = (float2)(0.0f,0.0f);\n"
"    float2 particle_proj_pos = particle->proj_pos;\n"
"    float2 particle_lambda = particle->lambda;\n"
"    float r_squared = search_radius*search_radius;\n"
"    for(int y=-1;y<=1;y++)\n"
"    {\n"
"        int y_offset = particle->id/dims.x + y;\n"
"        if(y_offset<0 || y_offset>=dims.y) continue;\n"
"        for(int x=-1;x<=1;x++)\n"
"        {\n"
"            int x_offset = particle->id%dims.x + x;\n"
"            if(x_offset<0 || x_offset>=dims.x) continue;\n"
"            uint array_offset = x_offset+y_offset*dims.x;\n"
"            uint range_begin = x_offset>0 ? histogram[array_offset - 1] : 0;\n"
"            uint range_end = histogram[array_offset];\n"
"            for(uint i=range_begin;i<range_end;i++)\n"
"            {\n"
"                if(i==index) continue;\n"
"                Particle* p = &particles[i];\n"
"                float2 distVec = -p->pos+particle_proj_pos;\n"
"                float d = dot(distVec,distVec);\n"
"                if(d<r_squared)\n"
"                {\n"
"                    float sCorr = -0.1f * pown((poly6_kernel(length(distVec),search_radius)/poly6_kernel(0.3*search_radius,search_radius)),4);\n"
"                    delta += (1.0f/resting_density) * (particle_lambda + p->lambda + sCorr) * spiky_kernel_gradient(length(distVec),search_radius,normalize(distVec));\n"
"                }\n"
"            }\n"
"        }\n"
"    }\n"
"    particle->delta = delta;\n"
"}\n"
"\n"
"void __kernel check_boundary_collision(global Particle* particles, global Boundary* boundaries, uint num_boundaries)\n"
"{\n"
"    int index = get_global_id(0);\n"
"    Particle* particle = &particles[index];\n"
"    for(int i=0;i<num_boundaries;i++)\n"
"    {\n"
"        Boundary* bound = &boundaries[i];\n"
"        float d1 = dot(bound->normal,bound->begin-particle->pos);\n"
"        float d2 = dot(bound->normal,bound->begin-((particle->proj_pos + particle->delta)));\n"
"        if((sign(d1)!=sign(d2)))\n"
"        {\n"
"            float si = sign(d2);\n"
"            float3 pos = intersect_line_line(bound->begin,\n"
"                                            bound->end,\n"
"                                            particle->pos,\n"
"                                            particle->proj_pos + particle->delta);\n"
"            if(pos.z<0.0 || pos.z>1.0) continue;\n"
"            particle->delta = -particle->proj_pos + ((float2)(pos.x,pos.y) + 0.1f * si*bound->normal);\n"
"        }\n"
"    }\n"
"}\n"
"void __kernel update_projected_positions(global Particle* particles)\n"
"{\n"
"    int index = get_global_id(0);\n"
"    Particle* particle = &particles[index];\n"
"    particle->proj_pos = particle->proj_pos + particle->delta;\n"
"}\n"
"\n"
"void __kernel update_particles(global Particle* particles, float timestep, uint2 dims, global uint* histogram, float visc, float radius)\n"
"{\n"
"    int index = get_global_id(0);\n"
"    Particle* particle = &particles[index];"
"    particle->vel = (1.0f/timestep) * (particle->proj_pos-particle->pos);\n"
"    barrier(CLK_GLOBAL_MEM_FENCE);\n"
"    float2 viscAcc = (float2)(0.0f,0.0f);"
"    float r_squared = radius * radius;\n"
"    for(int y=-1;y<=1;y++)\n"
"    {\n"
"        int y_offset = particle->id/dims.x + y;\n"
"        if(y_offset<0 || y_offset>=dims.y) continue;\n"
"        for(int x=-1;x<=1;x++)\n"
"        {\n"
"            int x_offset = particle->id%dims.x + x;\n"
"            if(x_offset<0 || x_offset>=dims.x) continue;\n"
"\n"
"            uint array_offset = x_offset+y_offset*dims.x;\n"
"            uint range_begin = x_offset>0 ? histogram[array_offset - 1] : 0;\n"
"            uint range_end = histogram[array_offset];\n"
"            for(uint i=range_begin;i<range_end;i++)\n"
"            {\n"
"                if(i==index) continue;\n"
"                Particle* p = &particles[i];\n"
"                float2 distVec = p->pos-particle->proj_pos;\n"
"                float d = dot(distVec,distVec);\n"
"                if(d<r_squared)\n"
"                {\n"
"                    float2 viscVel = p->vel - particle->vel;\n"
"                    viscAcc += viscVel*visc_kernel(distance(particle->pos,p->pos),radius);\n"
"                }\n"
"            }\n"
"        }\n"
"    }\n"
"    barrier(CLK_LOCAL_MEM_FENCE);\n"
"    particle->vel += visc*viscAcc;\n"
"    particle->pos = particle->proj_pos;"
"}\n"
"";

#endif
