module patmo_photoRates
contains

  !**************
  subroutine computePhotoRates(tau)
    use patmo_commons
    use patmo_parameters
    implicit none
    real*8,intent(in)::tau(photoBinsNumber,cellsNumber)

    !COS -> CO + S
    krate(:,38) = integrateXsec(1, tau(:,:))

    !O3 -> O2 + O
    krate(:,39) = integrateXsec(2, tau(:,:))

    !O2 -> O + O
    krate(:,40) = integrateXsec(3, tau(:,:))

    !CS2 -> CS + S
    krate(:,41) = integrateXsec(4, tau(:,:))

    !SO3 -> SO2 + O
    krate(:,42) = integrateXsec(5, tau(:,:))

    !SO2 -> SO + O
    krate(:,43) = integrateXsec(6, tau(:,:))

    !H2S -> SH + H
    krate(:,44) = integrateXsec(7, tau(:,:))

    !SO -> S + O
    krate(:,45) = integrateXsec(8, tau(:,:))

  end subroutine computePhotoRates

  !*************
  function integrateXsec(index,tau)
    use patmo_parameters
    use patmo_commons
    use patmo_constants
    implicit none
    integer,intent(in)::index
    real*8,intent(in)::tau(photoBinsNumber,cellsNumber)
    real*8::integrateXsec(cellsNumber),dE !declare dE JV
    integer::j,i !Added integer i JV

        dE=0.5!0.04956071 !nm value of dE also added JV
    !loop on cells (stride photobins)
    do j=1,cellsNumber
      integrateXsec(j) = 0.5*(sum(xsecAll(:,index)*photoFlux(:)*exp(-tau(:,j)/cos(57.3))*dE)) 
         !Photodissociation rate constant calculation changed in V4 
         !Added diurnal angle and avg planetary solar zenith angle JV
    end do

  end function integrateXsec

end module patmo_photoRates
