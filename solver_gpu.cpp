#include "solver_gpu.h"
#include "simulation.h"
#include <string.h>
#include <stdio.h>
#include <glm/gtc/type_ptr.hpp>

void create_context_callback(const char* message, const void* data, size_t data_size, void* userdata)
{
    printf("%s\n",message);
    exit(-1);
}

void create_kernel_callback(cl_program program, void* userdata)
{
    int build_status;
    clGetProgramBuildInfo(program, (cl_device_id)userdata, CL_PROGRAM_BUILD_STATUS, 0, &build_status, 0);
    if(build_status!=CL_BUILD_SUCCESS)
    {
        size_t length;
        char buffer[2048];
        clGetProgramBuildInfo(program, (cl_device_id)userdata, CL_PROGRAM_BUILD_LOG, sizeof(buffer), buffer, &length);
        printf("--- Build log ---\n %s\n",buffer);
        exit(build_status);
    }
}

SolverGPU::SolverGPU(unsigned int width, unsigned int height) : AbstractSolver(width, height)
{
    last_num_buckets = 0;
    last_num_particles = 0;
    last_num_boundaries = 0;

    // Get platform and device information
    cl_platform_id platform_id = NULL;
    cl_uint ret_num_devices;
    cl_uint ret_num_platforms;
    cl_int ret = clGetPlatformIDs(1, &platform_id, &ret_num_platforms);
    ret = clGetDeviceIDs( platform_id, CL_DEVICE_TYPE_GPU, 1,
            &device_id, &ret_num_devices);

    cl_context_properties props[] = {
        CL_CONTEXT_PLATFORM,(cl_context_properties) platform_id,
        0
    };

    context = clCreateContext(props,1,&device_id,create_context_callback,0,&ret);

    queue = clCreateCommandQueueWithProperties(context,device_id,0,&ret);

    size_t source_len = strlen(simulation_source);
    program = clCreateProgramWithSource(context,1,&simulation_source,&source_len,nullptr);

    clBuildProgram(program,1,&device_id,"-cl-std=CL2.0",create_kernel_callback,device_id);
    bucket_count = clCreateKernel(program,"bucket_count",nullptr);
    histogram = clCreateKernel(program,"histogram",nullptr);
    histogram2 = clCreateKernel(program,"histogram2",nullptr);
    histogram3 = clCreateKernel(program,"histogram3",nullptr);
    sort_particles = clCreateKernel(program,"sort_particles",nullptr);
    compute_lambda = clCreateKernel(program,"compute_lambda",nullptr);
    update_delta = clCreateKernel(program,"update_delta",nullptr);
    check_boundary_collision = clCreateKernel(program,"check_boundary_collision",nullptr);
    update_projected_positions = clCreateKernel(program,"update_projected_positions",nullptr);

    add_external_forces = clCreateKernel(program,"add_external_forces",nullptr);
    update_particles = clCreateKernel(program,"update_particles",nullptr);

}

SolverGPU::~SolverGPU()
{
    clReleaseMemObject(particles_buffer_front);
    clReleaseMemObject(particles_buffer_back);
    clReleaseMemObject(counter_buffer);
    clReleaseMemObject(histogram_buffer);
    clReleaseMemObject(scratch_buffer);
    clReleaseMemObject(scratch_buffer2);

    clReleaseKernel(bucket_count);
    clReleaseKernel(histogram);
    clReleaseKernel(histogram2);
    clReleaseKernel(histogram3);
    clReleaseKernel(sort_particles);
    clReleaseKernel(add_external_forces);
    clReleaseKernel(compute_lambda);
    clReleaseKernel(update_delta);
    clReleaseKernel(update_projected_positions);
    clReleaseKernel(update_particles);

    clReleaseProgram(program);
    clReleaseDevice(device_id);
    clReleaseContext(context);
}

void SolverGPU::rebuildParticleBufferGPU()
{
    if(last_num_particles!=0)
    {
        clReleaseMemObject(particles_buffer_front);
        clReleaseMemObject(particles_buffer_back);
    }
    particles_buffer_front = clCreateBuffer(context,CL_MEM_READ_WRITE|CL_MEM_USE_HOST_PTR,sizeof(Particle)*particles.size(),particles.data(),nullptr);
    particles_buffer_back = clCreateBuffer(context,CL_MEM_READ_WRITE|CL_MEM_USE_HOST_PTR,sizeof(Particle)*particles.size(),particles.data(),nullptr);
}

void SolverGPU::rebuildBoundaryBufferGPU()
{
    if(last_num_boundaries!=0)
    {
        clReleaseMemObject(boundaries_buffer);
    }
    boundaries_buffer = clCreateBuffer(context,CL_MEM_READ_ONLY|CL_MEM_USE_HOST_PTR,boundaries.size() * sizeof(Boundary), boundaries.data(), nullptr);
}

void SolverGPU::rebuildBucketBufferGPU()
{
    size_t num_buckets = std::ceil(domain_width/search_radius)*std::ceil(domain_height/search_radius);
    if(last_num_buckets!=0)
    {
        clReleaseMemObject(counter_buffer);
        clReleaseMemObject(histogram_buffer);
        clReleaseMemObject(scratch_buffer);
        clReleaseMemObject(scratch_buffer2);

    }
    counter_buffer = clCreateBuffer(context,CL_MEM_READ_WRITE,num_buckets,nullptr,nullptr);
    histogram_buffer = clCreateBuffer(context,CL_MEM_READ_WRITE,num_buckets,nullptr,nullptr);
    scratch_buffer = clCreateBuffer(context,CL_MEM_READ_WRITE,num_buckets,nullptr,nullptr);
    scratch_buffer2 = clCreateBuffer(context,CL_MEM_READ_WRITE,num_buckets,nullptr,nullptr);
}

void SolverGPU::solve()
{
    size_t num_particles = particles.size();
    size_t num_borders = boundaries.size();
    glm::uvec2 dims(std::ceil(domain_width/search_radius),std::ceil(domain_height/search_radius));
    size_t num_buckets = dims.x * dims.y;

    if(num_buckets!=last_num_buckets)
    {
        rebuildBucketBufferGPU();
    }
    if(num_particles!=last_num_particles)
    {
        rebuildParticleBufferGPU();
    }
    if(num_borders!=last_num_boundaries)
    {
        rebuildBoundaryBufferGPU();
    }

    unsigned int pattern = 0;

    clEnqueueFillBuffer(queue, counter_buffer, &pattern, sizeof(pattern),0, num_buckets * sizeof(unsigned int),0,nullptr,nullptr);
    clEnqueueFillBuffer(queue, scratch_buffer, &pattern, sizeof(pattern),0, num_buckets * sizeof(unsigned int),0,nullptr,nullptr);
    clEnqueueFillBuffer(queue, scratch_buffer2, &pattern, sizeof(pattern),0, std::ceil(num_buckets/128) * sizeof(unsigned int),0,nullptr,nullptr);

    clSetKernelArg(bucket_count,0,sizeof(cl_mem),&particles_buffer_front);
    clSetKernelArg(bucket_count,1,sizeof(cl_mem),&counter_buffer);
    clSetKernelArg(bucket_count,2,sizeof(glm::uvec2),glm::value_ptr(dims));
    clEnqueueNDRangeKernel(queue, bucket_count,1,nullptr,&num_particles, nullptr, 0, nullptr, nullptr);

    size_t local_size = 128;
    for(unsigned i=0;i<1;i++)
    {
        size_t global_size = num_buckets;
        clSetKernelArg(histogram,0,sizeof(cl_mem),&counter_buffer);
        clSetKernelArg(histogram,1,sizeof(cl_mem),&histogram_buffer);
        clSetKernelArg(histogram,2,sizeof(cl_mem),&scratch_buffer);
        clEnqueueNDRangeKernel(queue, histogram,1,nullptr,&global_size, &local_size, 0, nullptr, nullptr);

        global_size = std::ceil(num_buckets/128.0);
        clSetKernelArg(histogram2,0,sizeof(cl_mem),&scratch_buffer);
        clSetKernelArg(histogram2,1,sizeof(cl_mem),&scratch_buffer2);
        clEnqueueNDRangeKernel(queue, histogram2,1,nullptr,&global_size, &local_size, 0, nullptr, nullptr);

        global_size=num_buckets-128;
        clSetKernelArg(histogram3,0,sizeof(cl_mem),&scratch_buffer2);
        clSetKernelArg(histogram3,1,sizeof(cl_mem),&histogram_buffer);
        clSetKernelArg(histogram3,2,sizeof(cl_mem),&counter_buffer);

        clEnqueueNDRangeKernel(queue, histogram3,1,nullptr,&global_size, &local_size, 0, nullptr, nullptr);
    }

    clEnqueueFillBuffer(queue, scratch_buffer, &pattern, sizeof(pattern),0, num_buckets * sizeof(unsigned int),0,nullptr,nullptr);

    clSetKernelArg(sort_particles,0,sizeof(cl_mem),&particles_buffer_front);
    clSetKernelArg(sort_particles,1,sizeof(cl_mem),&particles_buffer_back);
    clSetKernelArg(sort_particles,2,sizeof(cl_mem),&histogram_buffer);
    clSetKernelArg(sort_particles,3,sizeof(cl_mem),&scratch_buffer);
    clEnqueueNDRangeKernel(queue, sort_particles,1,nullptr,&num_particles, nullptr, 0, nullptr, nullptr);

    clSetKernelArg(add_external_forces,0,sizeof(cl_mem),&particles_buffer_back);
    clSetKernelArg(add_external_forces,1,sizeof(float),&timestep);
    clSetKernelArg(add_external_forces,2,sizeof(glm::vec2),glm::value_ptr(gravity));
    clEnqueueNDRangeKernel(queue, add_external_forces,1,nullptr,&num_particles, nullptr, 0, nullptr, nullptr);

    for(unsigned int i=0;i<iterations;i++)
    {
        clSetKernelArg(compute_lambda,0,sizeof(cl_mem),&particles_buffer_back);
        clSetKernelArg(compute_lambda,1,sizeof(glm::uvec2),glm::value_ptr(dims));
        clSetKernelArg(compute_lambda,2,sizeof(cl_mem),&histogram_buffer);
        clSetKernelArg(compute_lambda,3,sizeof(float),&search_radius);
        clSetKernelArg(compute_lambda,4,sizeof(float),&resting_density);
        clSetKernelArg(compute_lambda,5,sizeof(float),&artificial_density);
        clEnqueueNDRangeKernel(queue, compute_lambda,1,nullptr,&num_particles, nullptr, 0, nullptr, nullptr);

        clSetKernelArg(update_delta,0,sizeof(cl_mem),&histogram_buffer);
        clSetKernelArg(update_delta,1,sizeof(glm::uvec2),glm::value_ptr(dims));
        clSetKernelArg(update_delta,2,sizeof(float),&search_radius);
        clSetKernelArg(update_delta,3,sizeof(float),&resting_density);
        clSetKernelArg(update_delta,4,sizeof(cl_mem),&particles_buffer_back);
        clEnqueueNDRangeKernel(queue, update_delta,1,nullptr,&num_particles, nullptr, 0, nullptr, nullptr);

        clSetKernelArg(check_boundary_collision,0,sizeof(cl_mem),&particles_buffer_back);
        clSetKernelArg(check_boundary_collision,1,sizeof(cl_mem),&boundaries_buffer);
        clSetKernelArg(check_boundary_collision,2,sizeof(unsigned int),&num_borders);
        clEnqueueNDRangeKernel(queue, check_boundary_collision,1,nullptr,&num_particles, nullptr, 0, nullptr, nullptr);

        clSetKernelArg(update_projected_positions,0,sizeof(cl_mem),&particles_buffer_back);
        clEnqueueNDRangeKernel(queue, update_projected_positions,1,nullptr,&num_particles, nullptr, 0, nullptr, nullptr);
    }
    clSetKernelArg(update_particles,0,sizeof(cl_mem),&particles_buffer_back);
    clSetKernelArg(update_particles,1,sizeof(float),&timestep);
    clEnqueueNDRangeKernel(queue, update_particles,1,nullptr,&num_particles, nullptr, 0, nullptr, nullptr);
    clEnqueueReadBuffer(queue, particles_buffer_back, true, 0, sizeof(Particle)*num_particles,particles.data(),0,nullptr,nullptr);
    clFinish(queue);

    cl_mem temp = particles_buffer_front;
    particles_buffer_front = particles_buffer_back;
    particles_buffer_back = temp;

    last_num_particles = num_particles;
    last_num_boundaries = num_borders;
    last_num_buckets = num_buckets;
}
