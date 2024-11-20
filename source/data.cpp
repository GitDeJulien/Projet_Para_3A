#ifndef _DATA_CPP

#include "data.h"
#include <fstream>
#include <iostream>
//#include <sstream>
#include <cmath>
#include <stdlib.h>
#include <string>



// Constructeur
Data::Data(std::string file_name)
{

    //Data file reading
    std::ifstream file(file_name);

    if (!file.is_open()) {
        std::cerr << "Error: Could not open file " << file_name << std::endl;
        return;
    }

    std::string line;
    while (std::getline(file, line)) {
        // Ignore empty lines and labels (lines starting with '#')
        if (line.empty() || line[0] == '#') continue;

        // Parse key-value pairs
        size_t equal_pos = line.find('=');
        if (equal_pos != std::string::npos) {
            std::string key = line.substr(0, equal_pos - 1);
            std::string string_value = line.substr(equal_pos + 1);
            //double double_value = std::stod(line.substr(equal_pos + 1));
            double double_value = atof(line.substr(equal_pos + 1).c_str());
            int int_value = atoi(line.substr(equal_pos + 1).c_str());


            // Assign values based on the key
            if (key == "Lx") this->_Lx = double_value;
            else if (key == "Ly") this->_Ly = double_value;
            else if (key == "Nx") this->_Nx = double_value;
            else if (key == "Ny") this->_Ny = double_value;
            else if (key == "t0") this->_t0 = double_value;
            else if (key == "niter") this->_niter = double_value;
            else if (key == "dt") this->_dt = double_value;
            else if (key == "cfl") this->_cfl = double_value;
            else if (key == "D") this->_diffusionCoeff = double_value;
            else if (key == "key_time_scheme") this->_key_TimeScheme = int_value;
            else if (key == "key_space_scheme") this->_key_SpaceScheme = int_value;
            else if (key == "key_RightBoundCond") this->_key_RightBoundCond = int_value;
            else if (key == "key_LeftBoundCond") this->_key_LeftBoundCond = int_value;
            else if (key == "key_UpBoundCond") this->_key_UpBoundCond = int_value;
            else if (key == "key_DownBoundCond") this->_key_DownBoundCond = int_value;
            else if (key == "key_SourceTerme") this->_key_SourceTerme = int_value;
            else if (key == "key_InitialCondition") this->_key_InitialCondition = int_value;
            else if (key == "outputPath") this->_outputPath = line.substr(equal_pos + 2).c_str();
            else {
                std::cerr << "Error: This parameters is not valid : " << key << std::endl;
                exit (EXIT_FAILURE);
            }

        }
    }

    file.close();

    this->_hx = this->_Lx / (this->_Nx+1);
    this->_hy = this->_Ly / (this->_Ny+1);

    if (this->_key_TimeScheme==1){ //Euler explicite => CFL
        this->_dt = this->_cfl*pow(this->_hy,2)*pow(this->_hx,2)/(2*this->_diffusionCoeff*(pow(this->_hx,2)+pow(this->_hy,2)));
    }

    if (this->_dt<pow(10,-5)){
        std::cout << "Time step lower bownd have been exceded" << std::endl;
        std::cout << "Please, reduce diffusion coeficient 'D' " << std::endl;
        exit(EXIT_FAILURE);
    }

    // Number of time iteration
    this->_tfinal = _dt*_niter; 
    int nb_iterations = (int)ceil((this->_tfinal-this->_t0)/this->_dt);
    // To set: _tfinal = _t0 + nb_iterations*_dt
    this->_dt = (this->_tfinal-this->_t0) / nb_iterations;

}

void Data::display_parameters() const {
    std::cout << "---------------------------------------------" << std::endl;
    std::cout << "-----------------PARAMETERS------------------" << std::endl;
    std::cout << "---------------------------------------------" << std::endl;
    std::cout << "Space Parameters:" << std::endl;
    std::cout << "Lx = " << _Lx << ", Ly = " << _Ly << std::endl;
    std::cout << "Nx = " << _Nx << ", Ny = " << _Ny << std::endl;
    std::cout << "hx = " << _hx << ", hy = " << _hy << std::endl;
    std::cout << "---------------------------------------------" << std::endl;


    std::cout << "Time Parameters:" << std::endl;
    std::cout << "t0 = " << _t0 << ", nb_iter = " << _niter << ", dt = " << _dt << ", tfinal = " << _tfinal << std::endl;
    std::cout << "---------------------------------------------" << std::endl;



    std::cout << "Boundary condition:" << std::endl;
    if (_key_LeftBoundCond == 1) std::cout << "Left = 0" << std::endl;
    else if (_key_LeftBoundCond == 2) std::cout << "Left = sin(x)+cos(y)" << std::endl;
    else if (_key_LeftBoundCond == 3) std::cout << "Left = 1" << std::endl;
    else {
        std::cerr << "Error: The Left boundary condition key" << _key_LeftBoundCond << "is not referenced" << std::endl;
        exit (EXIT_FAILURE);
    }

    if (_key_RightBoundCond == 1) std::cout << "Right = 0" << std::endl;
    else if (_key_RightBoundCond == 2) std::cout << "Right = sin(x)+cos(y)" << std::endl;
    else if (_key_RightBoundCond == 3) std::cout << "Right = 1" << std::endl;
    else {
        std::cerr << "Error: The Right boundary condition key" << _key_RightBoundCond << "is not referenced" << std::endl;
        exit (EXIT_FAILURE);
    }

    if (_key_UpBoundCond == 1) std::cout << "Up = 0" << std::endl;
    else if (_key_UpBoundCond == 2) std::cout << "Up = sin(x)+cos(y)" << std::endl;
    else if (_key_UpBoundCond == 3) std::cout << "Up = 1" << std::endl;
    else {
        std::cerr << "Error: The Up boundary condition key" << _key_UpBoundCond << "is not referenced" << std::endl;
        exit (EXIT_FAILURE);
    }

    if (_key_DownBoundCond == 1) std::cout << "Down = 0" << std::endl;
    else if (_key_DownBoundCond == 2) std::cout << "Down = sin(x)+cos(y)" << std::endl;
    else if (_key_DownBoundCond == 3) std::cout << "Down = 1" << std::endl;
    else {
        std::cerr << "Error: The Down boundary condition key" << _key_DownBoundCond << "is not referenced" << std::endl;
        exit (EXIT_FAILURE);
    }
    std::cout << "---------------------------------------------" << std::endl;
    

    if (_key_TimeScheme == 1) std::cout << "Time scheme: " << "Explicite Euler" << std::endl;
    else if (_key_TimeScheme == 2) std::cout << "Time scheme: " << "Implicte Euler " << std::endl;
    else if (_key_TimeScheme == 3) std::cout << "Time scheme: " << "Crank-Nicholson" << std::endl;
    else {
        std::cerr << "Error: The time scheme key" << _key_TimeScheme << "is not referenced" << std::endl;
        exit (EXIT_FAILURE);
    }
    std::cout << "---------------------------------------------" << std::endl;

    if (_key_SpaceScheme == 1) std::cout << "Space scheme: " << "Laplacian Centered" << std::endl;
    else {
        std::cerr << "Error: The time scheme key" << _key_SpaceScheme << "is not referenced" << std::endl;
        exit (EXIT_FAILURE);
    }
    std::cout << "---------------------------------------------" << std::endl;

    if (_key_SourceTerme == 1) std::cout << "Source terme: " << "f = 2*(x-x*x+y-y*y)" << std::endl;
    else if (_key_SourceTerme == 2) std::cout << "Source terme: " << "f = sin(x) + cos(y)" << std::endl;
    else if (_key_SourceTerme == 3) std::cout << "Source terme: " << "f = exp(-(x-Lx/2)^2)*exp(-(y-Ly/2)^2)*cos(pi/2*t)" << std::endl;
    else {
        std::cerr << "Error: The source terme key" << _key_SourceTerme << "is not referenced" << std::endl;
        exit (EXIT_FAILURE);
    }
    std::cout << "---------------------------------------------" << std::endl;

    if (_key_InitialCondition == 1) std::cout << "Initial condition: " << "u0 = ExactSol(x,y,0)" << std::endl;
    else if (_key_InitialCondition == 2) std::cout << "Initial condition: " << "u0 = sin(x) + cos(y)" << std::endl;
    else if (_key_InitialCondition == 3) std::cout << "Initial condition: " << "u0 = ..." << std::endl;
    else {
        std::cerr << "Error: The initial condition key" << _key_InitialCondition << "is not referenced" << std::endl;
        exit (EXIT_FAILURE);
    }
    std::cout << "---------------------------------------------" << std::endl;


    std::cout << "Diffusion coefficient:" << std::endl;
    std::cout << "D = " << _diffusionCoeff << std::endl;
    std::cout << "---------------------------------------------" << std::endl;
    std::cout << "output path = " << _outputPath << std::endl;
    std::cout <<  std::endl;

}


#define _DATA_CPP
#endif