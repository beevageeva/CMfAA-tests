! blast wave
module mod_usr
  use mod_hd
  implicit none
logical :: hmmm=.true.
contains

  subroutine usr_init()

!    unit_length        = 3.0857D18 !(parsecs)
!    unit_time          = 31557600D6 !(Myrs)
    print *, 'unit_time = ', unit_time
     print *, 'unit_length = ', unit_length
     print *, 'unit_density = ', unit_density
     print *, 'mp_cgs = ', mp_cgs
 
    usr_init_one_grid => initonegrid_usr
    usr_source          => special_source

    usr_aux_output    => specialvar_output
    usr_add_aux_names => specialvarnames_output 

    call set_coordinate_system("Cartesian")

    call hd_activate()

    end subroutine usr_init

  subroutine initonegrid_usr(ixI^L,ixO^L,w,x)
  ! initialize one grid
    use mod_physics
    integer, intent(in) :: ixI^L, ixO^L
    double precision, intent(in) :: x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)

    double precision :: rbs,xc1,xc2
  

    
       if (mype==0) then
          print *,'2D HD blast wave in Cartesian coordinate'

    end if
    w(ixO^S,rho_)=1.d0*10*mp_cgs
    w(ixO^S,p_)=w(ixO^S,rho_)*121.d0
    rbs=0.2d0
    xc1=(xprobmin1+xprobmax1)*0.5d0
    xc2=(xprobmin2+xprobmax2)*0.5d0
!    where((x(ixO^S,1)-xc1)**2+(x(ixO^S,2)-xc2)**2<rbs**2)
!      w(ixO^S,p_)=100.d0
!    endwhere
    w(ixO^S,mom(:))=0.d0


      call phys_to_conserved(ixI^L,ixO^L,w,x)
       

  end subroutine initonegrid_usr

  subroutine special_source(qdt,ixI^L,ixO^L,iw^LIM,qtC,wCT,qt,w,x)
    use mod_global_parameters
    integer, intent(in) :: ixI^L, ixO^L, iw^LIM
    double precision, intent(in) :: qdt, qtC, qt
    double precision, intent(in) :: x(ixI^S,1:ndim), wCT(ixI^S,1:nw)
    double precision, intent(inout) :: w(ixI^S,1:nw)
    double precision :: rad(ixI^S), xc1,xc2
  logical, save:: first=.true.
    ! use of special source as an internal boundary....
     xc1=(xprobmin1+xprobmax1)*0.5d0
     xc2=(xprobmin2+xprobmax2)*0.5d0

    rad(ixO^S) = dsqrt((x(ixO^S,1)-xc1)**2 + (x(ixO^S,2)-xc2)**2)
    where ( rad(ixO^S)< 1 )
      w(ixO^S,rho_)  =  w(ixO^S,rho_)+ qdt*2.5

      w(ixO^S,e_)= w(ixO^S,e_) + qdt
    end where

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
    w(ixO^S,nw+1)=tmp(ixO^S)/w(ixO^S,rho_)

  end subroutine specialvar_output

  subroutine specialvarnames_output(varnames)
  ! newly added variables need to be concatenated with the w_names/primnames string
    character(len=*) :: varnames
    varnames='Te'

  end subroutine specialvarnames_output

end module mod_usr
