! blast wave
module mod_usr
  use mod_hd
  implicit none
  double precision :: Edot1,Rhodot1,Edot2,Rhodot2,Rstar,Vstar

  double precision :: rhoISM,TISM


contains

  subroutine usr_init()

   use mod_global_parameters
 
    unit_length        = 3.0857D18
    unit_temperature   = 1.0d7**2.0d0/ (kb_cgs/mp_cgs)
    unit_numberdensity = 10.0**(-25)/mp_cgs

   
    usr_set_parameters  => initglobaldata_usr
    usr_init_one_grid   => initonegrid_usr
   ! usr_refine_grid     => specialrefine_grid
    usr_source          => special_source
    usr_aux_output      => specialvar_output
    usr_add_aux_names   => specialvarnames_output
   ! usr_var_for_errest  => myvar_for_errest



    call set_coordinate_system("Cartesian")

    call hd_activate()
 
  end subroutine usr_init
 
 subroutine initglobaldata_usr()
    use mod_global_parameters
    
    double precision :: physR, physV

    hd_gamma=5.0d0/3.0d0
    Rstar=8 !because the unit_length = 1pc
    Vstar = (4/3D0)*dpi*Rstar**3 
    physR = Rstar * unit_length
    physV = (4/3D0)*dpi*physR**3 

    Edot1 = ((1.5e51/physV)/unit_pressure)/((4e6 * const_years)/unit_time)
    Edot2 = ((1.7e51/physV)/unit_pressure)/((1e6 * const_years)/unit_time)
    Rhodot1  = (((15.0*const_msun)/physV)/unit_density)/((4e6 * const_years)/unit_time)
    Rhodot2  = (((40.0*const_msun)/physV)/unit_density)/((1e6 * const_years)/unit_time)

    rhoISM= (10.0 * mp_cgs)/unit_density 
    TISM  = 121.0/unit_temperature


  if(mype==0) then
      print *, 'unit_density = ', unit_density
      print *, 'unit_pressure = ', unit_pressure
      print *, 'unit_velocity = ', unit_velocity
      print *, 'unit_time = ', unit_time
  end if

  end subroutine initglobaldata_usr



  subroutine initonegrid_usr(ixI^L,ixO^L,w,x)
  ! initialize one grid
    use mod_physics
    use mod_global_parameters

    integer, intent(in) :: ixI^L, ixO^L
    double precision, intent(in) :: x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)

    !double precision :: rbs,xc1,xc2
    logical, save:: first=.true.

    if (first) then
       if (mype==0) then
          print *,'2D HD blast wave in Cartesian coordinate'
       end if
       first=.false.
    end if
    w(ixO^S,rho_)=rhoISM
    w(ixO^S,p_)=w(ixO^S,rho_)*TISM
  !  rbs=0.2d0
  !  xc1=(xprobmin1+xprobmax1)*0.5d0
  !  xc2=(xprobmin2+xprobmax2)*0.5d0
  !  where((x(ixO^S,1)-xc1)**2+(x(ixO^S,2)-xc2)**2<rbs**2)
  !    w(ixO^S,p_)=100.d0
  !  endwhere
    w(ixO^S,mom(:))=0.d0

    call phys_to_conserved(ixI^L,ixO^L,w,x)

  end subroutine initonegrid_usr


 subroutine special_source(qdt,ixI^L,ixO^L,iw^LIM,qtC,wCT,qt,w,x)
    use mod_global_parameters
    integer, intent(in) :: ixI^L, ixO^L, iw^LIM
    double precision, intent(in) :: qdt, qtC, qt
    double precision, intent(in) :: x(ixI^S,1:ndim), wCT(ixI^S,1:nw)
    double precision, intent(inout) :: w(ixI^S,1:nw)

    double precision :: rad(ixI^S), cosTh(ixI^S), sinTh(ixI^S)

    double precision :: dE,dRho
    ! use of special source as an internal boundary....

    if(qt < (5e6 * const_years)/unit_time) then 
      rad(ixO^S) = dsqrt(x(ixO^S,1)**2 + x(ixO^S,2)**2)
      cosTh(ixO^S) = x(ixO^S,2)/rad(ixO^S)
      sinTh(ixO^S) = x(ixO^S,1)/rad(ixO^S)
      if(qt < (4e6 * const_years)/unit_time) then 
        dE = Edot1
        dRho = Rhodot1
      else  
        dE = Edot2
        dRho = Rhodot2
      endif

      where (rad(ixO^S)< Rstar)
        w(ixO^S,rho_)  = w(ixO^S,rho_)+ qdt*dRho
        w(ixO^S,e_)  = w(ixO^S,e_)+ qdt*dE     
      end where
  endif
  end subroutine special_source


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
    w(ixO^S,nw+1)=tmp(ixO^S)/w(ixO^S,rho_)/Tscale

  end subroutine specialvar_output

  subroutine specialvarnames_output(varnames)
  ! newly added variables need to be concatenated with the w_names/primnames string
    character(len=*) :: varnames
    varnames='Te'

  end subroutine specialvarnames_output

end module mod_usr
