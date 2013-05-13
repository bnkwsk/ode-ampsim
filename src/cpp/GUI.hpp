#pragma once

#include <gtkmm.h>

#include "SimulationRunner.hpp"

class GUI
{
    Glib::RefPtr<Gtk::Application> app;
    Glib::RefPtr<Gtk::Builder> builder;
    Gtk::Window* mainWindow;

    Gtk::FileChooserButton* inputFileChooserButton;
    Gtk::FileChooserButton* impulseResponseFileChooserButton;
    Gtk::FileChooserDialog* outputFileChooserDialog;

    std::string inputFilePath;
    std::string impulseResponseFilePath;
    std::string outputFilePath;

    sigc::connection progressTimer;

    SimulationRunner &runner;

    void outputFileChooserDialogOKButtonClicked()
    {
        Gtk::Label* outputFileDisplayLabel;
        builder->get_widget("outputFileDisplayLabel", outputFileDisplayLabel);
        outputFilePath = outputFileChooserDialog->get_filename();
        outputFileChooserDialog->hide();
        outputFileDisplayLabel->set_text(outputFilePath);
    }

    void inputFileChanged()
    {
        inputFilePath = inputFileChooserButton->get_filename();
    }

    void impulseResponseFileChanged()
    {
        impulseResponseFilePath = impulseResponseFileChooserButton->get_filename();
    }

    bool updateProgressBar()
    {
        Gtk::ProgressBar *progressBar;
        builder->get_widget("progressBar", progressBar);
        if(runner.isRunning())
        {
            progressBar->set_fraction(runner.getProgress());
            return true;
        }
        else {
            progressBar->set_fraction(0.0);
            return false;
        }

    }

    void runButtonClicked()
    {
        Gtk::MessageDialog dialog(*mainWindow, "File settings error",
            false, Gtk::MESSAGE_ERROR);
        if(inputFilePath.length() == 0)
        {
            dialog.set_secondary_text("There is no input file selected");
            dialog.run();
            return;
        }
        if(impulseResponseFilePath.length() == 0)
        {
            dialog.set_secondary_text("There is no impulse response file selected");
            dialog.run();
            return;
        }
        if(outputFilePath.length() == 0)
        {
            dialog.set_secondary_text("There is no output file selected");
            dialog.run();
            return;
        }
        runner.run(inputFilePath.c_str(), impulseResponseFilePath.c_str(), outputFilePath.c_str());
        progressTimer = Glib::signal_timeout().connect(sigc::mem_fun(this, &GUI::updateProgressBar), 500);
    }

    void stopSimulationRequested()
    {
        runner.interrupt();
        if (!progressTimer.empty())
            progressTimer.disconnect();
        runner.join();
        updateProgressBar();
    }

    void quitCallback()
    {
        stopSimulationRequested();
        app->quit();
    }

    public:
        GUI(Glib::RefPtr<Gtk::Application> _app, SimulationRunner &_runner) : app(_app), runner(_runner)
        {
            builder = Gtk::Builder::create_from_file("gui.glade");
            builder->get_widget("mainWindow", mainWindow);

            Gtk::Button* closeButton;
            builder->get_widget("closeButton", closeButton);
            closeButton->signal_clicked().connect(
                sigc::mem_fun(this, &GUI::stopSimulationRequested));
            mainWindow->signal_hide().connect(
                sigc::mem_fun(this, &GUI::quitCallback));

            builder->get_widget("inputFile", inputFileChooserButton);
            inputFileChooserButton->signal_selection_changed().connect(
                sigc::mem_fun(this, &GUI::inputFileChanged));
            builder->get_widget("impulseFile", impulseResponseFileChooserButton);
            impulseResponseFileChooserButton->signal_selection_changed().connect(
                sigc::mem_fun(this, &GUI::impulseResponseFileChanged));

            Gtk::Button* outputFileButton;
            builder->get_widget("outputFileButton", outputFileButton);
            builder->get_widget(
                "outputFileChooserDialog", outputFileChooserDialog);
            outputFileButton->signal_clicked().connect(
                sigc::mem_fun(
                    outputFileChooserDialog, &Gtk::FileChooserDialog::show));

            Gtk::Button* outputFileChooserDialogCloseButton;
            builder->get_widget(
                "outputFileChooserDialogCloseButton",
                outputFileChooserDialogCloseButton);
            outputFileChooserDialogCloseButton->signal_clicked().connect(
                sigc::mem_fun(
                    outputFileChooserDialog, &Gtk::FileChooserDialog::hide));

            Gtk::Button* outputFileChooserDialogOKButton;
            builder->get_widget(
                "outputFileChooserDialogOKButton",
                outputFileChooserDialogOKButton);
            outputFileChooserDialogOKButton->signal_clicked().connect(
                sigc::mem_fun(
                    this, &GUI::outputFileChooserDialogOKButtonClicked));
            outputFileChooserDialog->signal_file_activated().connect(
                sigc::mem_fun(
                    this, &GUI::outputFileChooserDialogOKButtonClicked));

            Gtk::Button* runButton;
            builder->get_widget("runButton", runButton);
            runButton->signal_clicked().connect(
                sigc::mem_fun(this, &GUI::runButtonClicked));

            app->run(*mainWindow);
        }
};