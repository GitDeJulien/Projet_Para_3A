module data_mod

    use precision
    use toml_parser
    implicit none

    public :: config_data

    type, public :: DataType

    ![SPACE] (2D)
    real(pr) :: Lx
    real(pr) :: Ly
    integer  :: Nx
    integer  :: Ny
    integer  :: N_pts
    real(pr) :: hx
    real(pr) :: hy
    
    ![TIME] (implicite scheme)
    real(pr) :: t0
    integer  :: niter
    real(pr) :: dt
    real(pr) :: tfinal
    real(pr) :: cfl
    
    ![Diffusion coefficient]
    real(pr) :: D

    ![BC SCHWARZ METHODE] (1:Dirichlet, 2:Robin)
    integer  :: BC_Schwarz
    integer  :: overlap
    
    ![CAS TEST]
    integer  :: cas

    ![PARALLEL]
    integer  :: rank
    integer  :: n_proc
    integer  :: jbeg
    integer  :: jend
    integer  :: l_top
    integer  :: l_bot



    end type DataType

contains

    subroutine config_data(data, filename)
        character(len=*), intent(in) :: filename
        type(DataType), intent(inout) :: data

        !Local
        integer :: nb_iteration

        call parse_toml(filename, "Lx", data%Lx)
        call parse_toml(filename, "Ly", data%Ly)
        call parse_toml(filename, "Nx", data%Nx)
        call parse_toml(filename, "Ny", data%Ny)

        call parse_toml(filename, "t0", data%t0)
        call parse_toml(filename, "niter", data%niter)
        call parse_toml(filename, "dt", data%dt)
        call parse_toml(filename, "cfl", data%cfl)

        call parse_toml(filename, "D", data%D)
        call parse_toml(filename, "BC_Schwarz", data%BC_Schwarz)
        call parse_toml(filename, "overlap", data%overlap)
        call parse_toml(filename, "cas", data%cas)

        data%hx = data%Lx / (data%Nx+1)
        data%hy = data%Ly / (data%Ny+1)
        data%N_pts = data%Nx*data%Ny

        data%tfinal = data%dt*data%niter
        nb_iteration = CEILING((data%tfinal-data%t0)/data%dt, kind=selected_int_kind(8)) 
        data%dt = (data%tfinal-data%t0) / nb_iteration


    end subroutine config_data
    
end module data_mod