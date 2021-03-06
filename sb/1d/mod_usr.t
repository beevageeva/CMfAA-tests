module mod_usr
  use mod_hd
  implicit none

contains

  subroutine usr_init()
    unit_length        = 3.0857d18
    unit_temperature   = 1d7
    unit_numberdensity = 1d0
    usr_init_one_grid => initonegrid_usr
    usr_aux_output    => specialvar_output
    usr_add_aux_names => specialvarnames_output 

    call set_coordinate_system("spherical")
    !call set_coordinate_system("Cartesian_1D")
    call hd_activate()

  end subroutine usr_init

  subroutine initonegrid_usr(ixI^L,ixO^L,w,x)
  ! initialize one grid
    use mod_physics
    use mod_constants
    integer, intent(in) :: ixI^L, ixO^L
    double precision, intent(in) :: x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)

    double precision :: rbs,xc1
    rbs=1D18/unit_length

    !print*, "UNITTIME", unit_time

    w(ixO^S,rho_)=1.d0*mp_cgs/unit_density
    w(ixO^S,p_)=1d2*mp_cgs/unit_pressure
    xc1=(xprobmin1+xprobmax1)*0.5d0
    where((x(ixO^S,1)-xc1)**2<rbs**2)
      w(ixO^S,p_)=((1d51/(4d0/3d0 * dpi * 1d18**3))*(hd_gamma-1.0)) /unit_pressure
    endwhere
    w(ixO^S,mom(:))=0.d0

    call phys_to_conserved(ixI^L,ixO^L,w,x)

  end subroutine initonegrid_usr

  subroutine specialvar_output(ixI^L,ixO^L,w,x,normconv)
  ! this subroutine can be used in convert, to add auxiliary variables to the
  ! converted output file, for further analysis using tecplot, paraview, ....
  ! these auxiliary values need to be stored in the nw+1:nw+nwauxio slots
  !
  ! the array normconv can be filled in the (nw+1:nw+nwauxio) range with
  ! corresponding normalization values (default value 1)
    use mod_physics
    integer, intent(in)                :: ixI^L,ixO^L
    double precision, intent(in)       :: x(ixI^S,1:ndim)
    double precision                   :: w(ixI^S,nw+nwauxio)
    double precision                   :: normconv(0:nw+nwauxio)

    double precision                   :: tmp(ixI^S) 

    call phys_get_pthermal(w,x,ixI^L,ixO^L,tmp)
    ! output the temperature p/rho
    w(ixO^S,nw+1)=tmp(ixO^S)/w(ixO^S,rho_)

  end subroutine specialvar_output

  subroutine specialvarnames_output(varnames)
  ! newly added variables need to be concatenated with the w_names/primnames string
    character(len=*) :: varnames
    varnames='Te'

  end subroutine specialvarnames_output

end module mod_usr
