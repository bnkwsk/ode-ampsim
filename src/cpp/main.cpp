#define VIENNACL_WITH_OPENCL

#include <iostream>

#include "SimulationRunner.hpp"
#include "GUI.hpp"


int main(int argc, char **argv)
{
    Glib::RefPtr<Gtk::Application> app = Gtk::Application::create(
        argc, argv, "org.bnkwsk.simulator");
    SimulationRunner runner;
    GUI gui(app, runner);
    runner.join();
    return 0;
}