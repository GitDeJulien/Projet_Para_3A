#ifndef _DATA_H
#define _DATA_H

#include <string>
#include <iostream>

class Data {

    private:
        
        //Space
        double _Lx;
        double _Ly;
        int _Nx;
        int _Ny;
        double _hx;
        double _hy;
        int _N_pts;

        //Time
        double _t0;
        int _niter;
        double _dt;
        double _cfl;
        double _tfinal;

        //Time scheme
        int _key_TimeScheme;

        //Space scheme
        int _key_SpaceScheme;

        //Boundary conditions key
        int _key_RightBoundCond;
        int _key_LeftBoundCond;
        int _key_DownBoundCond;
        int _key_UpBoundCond;

        //Boundary condition for Schwarz method
        int _key_Schwarz_Bounds;

        //Source terme key
        int _key_SourceTerme;

        //Initial condition key
        int _key_InitialCondition;

        //Diffusion coefficient D
        double _diffusionCoeff;

        //Output path
        std::string _outputPath;

        


    public:

        // Constructeur
        Data(std::string file_name);

        //Getter
        const double & Get_diffusion_coeff() const {return _diffusionCoeff;};
        const double & Get_Lx() const {return _Lx;};
        const double & Get_Ly() const {return _Ly;};
        const int & Get_Nx() const {return _Nx;};
        const int & Get_Ny() const {return _Ny;};
        const int & Get_N_pts() const {return _N_pts;};
        const double & Get_hx() const {return _hx;};
        const double & Get_hy() const {return _hy;};
        const double & Get_t0() const {return _t0;};
        const int & Get_niter() const {return _niter;};
        const double & Get_dt() const {return _dt;};
        const double & Get_cfl() const {return _cfl;};
        const double & Get_tfinal() const {return _tfinal;};
        const int & Get_TimeScheme() const {return _key_TimeScheme;};
        const int & Get_SpaceScheme() const {return _key_SpaceScheme;};
        const int & Get_key_RightBoundCond() const {return _key_RightBoundCond;};
        const int & Get_key_LeftBoundCond() const {return _key_LeftBoundCond;};
        const int & Get_key_UpBoundCond() const {return _key_UpBoundCond;};
        const int & Get_key_DownBoundCond() const {return _key_DownBoundCond;};
        const int & Get_key_Schwarz_Bounds() const {return _key_Schwarz_Bounds;};
        const int & Get_key_SourceTerme() const {return _key_SourceTerme;};
        const int & Get_key_InitialCondition() const {return _key_InitialCondition;};
        const std::string & Get_outputPath() const {return _outputPath;};


        //Display parameters
        void display_parameters() const;

};

#endif