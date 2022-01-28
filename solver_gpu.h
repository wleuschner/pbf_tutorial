#ifndef SOLVER_GPU_H
#define SOLVER_GPU_H
#include "abstractsolver.h"
#include "particle.h"
#include "boundary.h"
#define CL_TARGET_OPENCL_VERSION 210
#include <CL/cl.h>

class SolverGPU : public AbstractSolver
{
public:
    SolverGPU(unsigned int width, unsigned int height);
    ~SolverGPU();

    void solve() override;

private:
    void rebuildParticleBufferGPU();
    void rebuildBoundaryBufferGPU();
    void rebuildBucketBufferGPU();

    unsigned int last_num_buckets;
    unsigned int last_num_particles;
    unsigned int last_num_boundaries;

    cl_context context;
    cl_command_queue queue;
    cl_device_id device_id;
    cl_program program;

    cl_kernel bucket_count;
    cl_kernel histogram;
    cl_kernel histogram2;
    cl_kernel histogram3;
    cl_kernel sort_particles;
    cl_kernel add_external_forces;
    cl_kernel compute_lambda;
    cl_kernel update_delta;
    cl_kernel check_boundary_collision;
    cl_kernel update_projected_positions;
    cl_kernel update_particles;


    cl_mem particles_buffer_front;
    cl_mem particles_buffer_back;
    cl_mem counter_buffer;
    cl_mem histogram_buffer;
    cl_mem scratch_buffer;
    cl_mem scratch_buffer2;

    cl_mem boundaries_buffer;
};

#endif // SOLVER_GPU_H
