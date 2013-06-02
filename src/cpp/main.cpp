#define VIENNACL_WITH_OPENCL

#include "viennacl/ocl/backend.hpp"

#include "SimulationRunner.hpp"
#include "GUI.hpp"

int main(int argc, char **argv)
{
    Glib::RefPtr<Gtk::Application> app = Gtk::Application::create(
        argc, argv, "org.bnkwsk.simulator");

    const char *complexPointwiseMultiplyAddCode = 
"__kernel void ComplexPointwiseMultiplyAdd(__global const float *a, __global const float *b, __global float *c, const unsigned int parts)"
"{"
"    const int id = get_global_id(0); /* block position 0 -> 2B - 1 */"
"    if(id % 2 == 0)"
"        for(unsigned int p = 0; p < parts; ++p)"
"            c[id] += a[id] * b[p * parts + id] - a[id + 1] * b[p * parts + id + 1];"
"    else"
"        for(unsigned int p = 0; p < parts; ++p)"
"            c[id] += a[id - 1] * b[p * parts + id] + a[id] * b[p * parts + id - 1];"
"}";
    viennacl::ocl::program & program = viennacl::ocl::current_context().add_program(complexPointwiseMultiplyAddCode, "ComplexPointwiseMultiplyAdd");
    program.add_kernel("ComplexPointwiseMultiplyAdd");

    SimulationRunner runner;
    GUI gui(app, runner);
    runner.join();

    return 0;
}
