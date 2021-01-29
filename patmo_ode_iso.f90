module patmo_ode
contains
  subroutine fex(neq,tt,nin,dy)
    use patmo_commons
    use patmo_constants
    use patmo_parameters
    use patmo_utils
    use patmo_budget, only: budget
    implicit none
    integer,intent(in)::neq
    real*8,intent(in)::tt,nin(neqAll)
    real*8,intent(out)::dy(neqAll)
    real*8::d_hp(cellsNumber,speciesNumber)
    real*8::d_hm(cellsNumber,speciesNumber)
    real*8::k_hp(cellsNumber)
    real*8::k_hm(cellsNumber)
    real*8::dzz_hp(cellsNumber),dzz_hm(cellsNumber)
    real*8::kzz_hp(cellsNumber),kzz_hm(cellsNumber)
    real*8::prem(cellsNumber)
    real*8::n(cellsNumber,speciesNumber)
    real*8::dn(cellsNumber,speciesNumber)
    real*8::Tgas(cellsNumber)
    real*8::n_p(cellsNumber,speciesNumber)
    real*8::n_m(cellsNumber,speciesNumber)
    real*8::m(speciesNumber),ngas(cellsNumber)
    real*8::ngas_hp(cellsNumber),ngas_hm(cellsNumber)
    real*8::ngas_p(cellsNumber),ngas_m(cellsNumber)
    real*8::Tgas_hp(cellsNumber),Tgas_hm(cellsNumber)
    real*8::Tgas_p(cellsNumber),Tgas_m(cellsNumber)
    real*8::ngas_hpp(cellsNumber)
    real*8::ngas_hmm(cellsNumber)
    real*8::ngas_hpz(cellsNumber)
    real*8::ngas_hmz(cellsNumber)
    real*8::therm_hp(cellsNumber)
    real*8::therm_hm(cellsNumber)
    real*8::dzzh_hp(cellsNumber)
    real*8::dzzh_hm(cellsNumber)
    real*8::iTgas_hp(cellsNumber)
    real*8::iTgas_hm(cellsNumber)
    real*8::OCStotal(cellsNumber),SO2total(cellsNumber)
    real*8::H2Stotal(cellsNumber),CH3SCH3total(cellsNumber)
    integer::i,j

    !get mass of individual species
    m(:) = getSpeciesMass()

    !roll chemistry
    do i=1,speciesNumber
      n(:,i) = nin((i-1)*cellsNumber+1:(i*cellsNumber))
    end do

    !local copy of Tgas
    Tgas(:) = nin((positionTgas-1)*cellsNumber+1:(positionTgas*cellsNumber))
    ngas(:) = nTotAll(:)

    !forward grid points
    do j=1,cellsNumber-1
      dzz_hp(j) = .5d0*(diffusionDzz(j)+diffusionDzz(j+1))
      kzz_hp(j) = .5d0*(eddyKzz(j)+eddyKzz(j+1))
      Tgas_hp(j) = .5d0*(Tgas(j)+Tgas(j+1))
      Tgas_p(j) = Tgas(j+1)
      ngas_p(j) = ngas(j+1)
      ngas_hp(j) = .5d0*(ngas(j)+ngas(j+1))
      n_p(j,:) = n(j+1,:)
    end do

    !forward grid points: boundary conditions
    dzz_hp(cellsNumber) = 0d0
    kzz_hp(cellsNumber) = 0d0
    Tgas_hp(cellsNumber) = Tgas_hp(cellsNumber-1)
    Tgas_p(cellsNumber) = Tgas_p(cellsNumber-1)
    ngas_p(cellsNumber) = ngas_p(cellsNumber-1)
    ngas_hp(cellsNumber) = ngas_hp(cellsNumber-1)
    n_p(cellsNumber,:) = n_p(cellsNumber-1,:)

    !bakcward grid points
    do j=2,cellsNumber
      dzz_hm(j) = .5d0*(diffusionDzz(j)+diffusionDzz(j-1))
      kzz_hm(j) = .5d0*(eddyKzz(j)+eddyKzz(j-1))
      Tgas_hm(j) = .5d0*(Tgas(j)+Tgas(j-1))
      Tgas_m(j) = Tgas(j-1)
      ngas_m(j) = ngas(j-1)
      ngas_hm(j) = .5d0*(ngas(j)+ngas(j-1))
      n_m(j,:) = n(j-1,:)
    end do

    !backward grid points: boundary conditions
    dzz_hm(1) = 0d0
    kzz_hm(1) = 0d0
    Tgas_hm(1) = Tgas_hm(2)
    Tgas_m(1) = Tgas_m(2)
    ngas_m(1) = ngas_m(2)
    ngas_hm(1) = ngas_hm(2)
    n_m(1,:) = n_m(2,:)

    !eqn.24 of Rimmer+Helling (2015), http://arxiv.org/abs/1510.07052
    therm_hp(:) = thermalDiffusionFactor/Tgas_hp(:)*(Tgas_p(:)-Tgas(:))
    therm_hm(:) = thermalDiffusionFactor/Tgas_hm(:)*(Tgas_m(:)-Tgas(:))
    dzzh_hp(:) = 0.5d0*dzz_hp(:)*idh2(:)
    dzzh_hm(:) = 0.5d0*dzz_hm(:)*idh2(:)
    iTgas_hp(:) = 1d0/Tgas_hp(:)
    iTgas_hm(:) = 1d0/Tgas_hm(:)
    do i=1,speciesNumber
      prem(:) = (meanMolecularMass-m(i))*gravity/kboltzmann*gridSpace(:)
      d_hp(:,i) =  dzzh_hp(:) &
          * (prem(:)*iTgas_hp(:) &
          - therm_hp(:))
      d_hm(:,i) = dzzh_hm(:) &
          * (prem(:)*iTgas_hm(:) &
          - therm_hm(:))
    end do

    k_hp(:) = (kzz_hp(:)+dzz_hp(:))*idh2(:)
    k_hm(:) = (kzz_hm(:)+dzz_hm(:))*idh2(:)

    !chemistry module
    dn(:,:) = 0d0

!    if (end_of_run) then
!        budget(:,patmo_idx_33SO2,1) = &
!        + krate(:,58)*n(:,patmo_idx_33SO)*n(:,patmo_idx_O3) &
!        + krate(:,62)*n(:,patmo_idx_33SO)*n(:,patmo_idx_O2) &
!        + krate(:,66)*n(:,patmo_idx_33SO)*n(:,patmo_idx_OH) &
!        + krate(:,90)*n(:,patmo_idx_33HSO)*n(:,patmo_idx_O2) &
!        + krate(:,98)*n(:,patmo_idx_33HSO2)*n(:,patmo_idx_O2) &
!        + krate(:,118)*n(:,patmo_idx_33H2SO4) &
!        + krate(:,128)*n(:,patmo_idx_33CH3SCH3)*n(:,patmo_idx_O) &
!        + krate(:,132)*n(:,patmo_idx_33CH3SCH3)*n(:,patmo_idx_OH) &
!        + krate(:,136)*n(:,patmo_idx_33CH3SCH3)*n(:,patmo_idx_OH) &
!        + krate(:,140)*n(:,patmo_idx_33SO2_1) &
!        + krate(:,144)*n(:,patmo_idx_33SO2_1) &
!        + krate(:,160)*n(:,patmo_idx_33SO2_3) &
!        + krate(:,164)*n(:,patmo_idx_33SO2_3) &
!        + krate(:,188)*n(:,patmo_idx_33SO3) &
!        + krate(:,290)*n(:,patmo_idx_OH)*n(:,patmo_idx_33SO3) &
!        + krate(:,294)*n(:,patmo_idx_33SO3)*n(:,patmo_idx_O2) &
!        + krate(:,314)*n(:,patmo_idx_33SO3) &
!        + krate(:,318)*n(:,patmo_idx_33HSO3) &
!        + krate(:,332)*n(:,patmo_idx_33SO4)
!
!        budget(:,patmo_idx_33SO2,2) = &
!        - krate(:,82)*n(:,patmo_idx_33SO2)*n(:,patmo_idx_HO2) &
!        - krate(:,86)*n(:,patmo_idx_33SO2)*n(:,patmo_idx_O3) &
!        - krate(:,106)*n(:,patmo_idx_33SO2)*n(:,patmo_idx_O) &
!        - krate(:,110)*n(:,patmo_idx_33SO2)*n(:,patmo_idx_OH) &
!        - krate(:,124)*n(:,patmo_idx_33SO2) &
!        - krate(:,198)*n(:,patmo_idx_33SO2) &
!        - krate(:,199)*n(:,patmo_idx_33SO2) &
!        - krate(:,200)*n(:,patmo_idx_33SO2) &
!        - krate(:,266)*n(:,patmo_idx_33SO2)*n(:,patmo_idx_O2) &
!        - krate(:,270)*n(:,patmo_idx_33SO2)*n(:,patmo_idx_O) &
!        - krate(:,274)*n(:,patmo_idx_33SO2)*n(:,patmo_idx_H) &
!        - krate(:,298)*n(:,patmo_idx_33SO2)*n(:,patmo_idx_OH) &
!        - krate(:,306)*n(:,patmo_idx_HO2)*n(:,patmo_idx_33SO2) &
!        - krate(:,326)*n(:,patmo_idx_33SO2)*n(:,patmo_idx_OH)*n(:,patmo_idx_OH) &
!        - krate(:,336)*n(:,patmo_idx_33SO2) &
!        - krate(:,340)*n(:,patmo_idx_33SO2) &
!        - krate(:,344)*n(:,patmo_idx_33SO2)*n(:,patmo_idx_CH4O3S) &
!        - krate(:,348)*n(:,patmo_idx_33SO2)*n(:,patmo_idx_N2) &
!        - krate(:,352)*n(:,patmo_idx_33SO2) &
!        - krate(:,368)*n(:,patmo_idx_33SO2) &
!        - krate(:,372)*n(:,patmo_idx_33SO2)*n(:,patmo_idx_N2)
!    endif
!
    dn(:,patmo_idx_33SO2) = &
        + krate(:,58)*n(:,patmo_idx_33SO)*n(:,patmo_idx_O3) &
        + krate(:,62)*n(:,patmo_idx_33SO)*n(:,patmo_idx_O2) &
        + krate(:,66)*n(:,patmo_idx_33SO)*n(:,patmo_idx_OH) &
        - krate(:,82)*n(:,patmo_idx_33SO2)*n(:,patmo_idx_HO2) &
        - krate(:,86)*n(:,patmo_idx_33SO2)*n(:,patmo_idx_O3) &
        + krate(:,90)*n(:,patmo_idx_33HSO)*n(:,patmo_idx_O2) &
        + krate(:,98)*n(:,patmo_idx_33HSO2)*n(:,patmo_idx_O2) &
        - krate(:,106)*n(:,patmo_idx_33SO2)*n(:,patmo_idx_O) &
        - krate(:,110)*n(:,patmo_idx_33SO2)*n(:,patmo_idx_OH) &
        + krate(:,118)*n(:,patmo_idx_33H2SO4) &
        - krate(:,124)*n(:,patmo_idx_33SO2) &
        + 0.99*(krate(:,128)*n(:,patmo_idx_33CH3SCH3)*n(:,patmo_idx_O)) &
        + 0.99*(krate(:,132)*n(:,patmo_idx_33CH3SCH3)*n(:,patmo_idx_OH)) &
        + 0.99*(krate(:,136)*n(:,patmo_idx_33CH3SCH3)*n(:,patmo_idx_OH)) &
        + krate(:,140)*n(:,patmo_idx_33SO2_1) &
        + krate(:,144)*n(:,patmo_idx_33SO2_1) &
        + krate(:,160)*n(:,patmo_idx_33SO2_3) &
        + krate(:,164)*n(:,patmo_idx_33SO2_3) &
        + krate(:,188)*n(:,patmo_idx_33SO3) &
        - krate(:,198)*n(:,patmo_idx_33SO2) &
        - krate(:,199)*n(:,patmo_idx_33SO2) &
        - krate(:,200)*n(:,patmo_idx_33SO2) &
        - krate(:,266)*n(:,patmo_idx_33SO2)*n(:,patmo_idx_O2) &
        - krate(:,270)*n(:,patmo_idx_33SO2)*n(:,patmo_idx_O) &
        - krate(:,274)*n(:,patmo_idx_33SO2)*n(:,patmo_idx_H) &
        + krate(:,290)*n(:,patmo_idx_OH)*n(:,patmo_idx_33SO3) &
        + krate(:,294)*n(:,patmo_idx_33SO3)*n(:,patmo_idx_O2) &
        - krate(:,298)*n(:,patmo_idx_33SO2)*n(:,patmo_idx_OH) &
        - krate(:,306)*n(:,patmo_idx_HO2)*n(:,patmo_idx_33SO2) &
        + krate(:,314)*n(:,patmo_idx_33SO3) &
        + krate(:,318)*n(:,patmo_idx_33HSO3) &
        - krate(:,326)*n(:,patmo_idx_33SO2)*n(:,patmo_idx_OH)*n(:,patmo_idx_OH) &
        + krate(:,332)*n(:,patmo_idx_33SO4) &
        - krate(:,336)*n(:,patmo_idx_33SO2) &
        - krate(:,340)*n(:,patmo_idx_33SO2) &
        - krate(:,344)*n(:,patmo_idx_33SO2)*n(:,patmo_idx_CH4O3S) &
        - krate(:,348)*n(:,patmo_idx_33SO2)*n(:,patmo_idx_N2) &
        - krate(:,352)*n(:,patmo_idx_33SO2) &
        - krate(:,368)*n(:,patmo_idx_33SO2) &
        - krate(:,372)*n(:,patmo_idx_33SO2)*n(:,patmo_idx_N2)

    dn(:,patmo_idx_34HSO2) = &
        - krate(:,99)*n(:,patmo_idx_34HSO2)*n(:,patmo_idx_O2) &
        + krate(:,307)*n(:,patmo_idx_HO2)*n(:,patmo_idx_34SO2)

    dn(:,patmo_idx_32HSO3) = &
        - krate(:,101)*n(:,patmo_idx_32HSO3)*n(:,patmo_idx_O2) &
        + krate(:,109)*n(:,patmo_idx_32SO2)*n(:,patmo_idx_OH) &
        + krate(:,309)*n(:,patmo_idx_HO2)*n(:,patmo_idx_32SO3) &
        - krate(:,317)*n(:,patmo_idx_32HSO3)

    dn(:,patmo_idx_HO2) = 0d0
        ! - krate(:,41)*n(:,patmo_idx_32H2S)*n(:,patmo_idx_HO2) &
        ! - krate(:,42)*n(:,patmo_idx_33H2S)*n(:,patmo_idx_HO2) &
        ! - krate(:,43)*n(:,patmo_idx_34H2S)*n(:,patmo_idx_HO2) &
        ! - krate(:,44)*n(:,patmo_idx_36H2S)*n(:,patmo_idx_HO2) &
        ! - krate(:,81)*n(:,patmo_idx_32SO2)*n(:,patmo_idx_HO2) &
        ! - krate(:,82)*n(:,patmo_idx_33SO2)*n(:,patmo_idx_HO2) &
        ! - krate(:,83)*n(:,patmo_idx_34SO2)*n(:,patmo_idx_HO2) &
        ! - krate(:,84)*n(:,patmo_idx_36SO2)*n(:,patmo_idx_HO2) &
        ! + krate(:,97)*n(:,patmo_idx_32HSO2)*n(:,patmo_idx_O2) &
        ! + krate(:,98)*n(:,patmo_idx_33HSO2)*n(:,patmo_idx_O2) &
        ! + krate(:,99)*n(:,patmo_idx_34HSO2)*n(:,patmo_idx_O2) &
        ! + krate(:,100)*n(:,patmo_idx_36HSO2)*n(:,patmo_idx_O2) &
        ! + krate(:,101)*n(:,patmo_idx_32HSO3)*n(:,patmo_idx_O2) &
        ! + krate(:,102)*n(:,patmo_idx_33HSO3)*n(:,patmo_idx_O2) &
        ! + krate(:,103)*n(:,patmo_idx_34HSO3)*n(:,patmo_idx_O2) &
        ! + krate(:,104)*n(:,patmo_idx_36HSO3)*n(:,patmo_idx_O2) &
        ! + krate(:,249)*n(:,patmo_idx_H2O)*n(:,patmo_idx_32HSO) &
        ! + krate(:,250)*n(:,patmo_idx_H2O)*n(:,patmo_idx_33HSO) &
        ! + krate(:,251)*n(:,patmo_idx_H2O)*n(:,patmo_idx_34HSO) &
        ! + krate(:,252)*n(:,patmo_idx_H2O)*n(:,patmo_idx_36HSO) &
        ! + krate(:,289)*n(:,patmo_idx_OH)*n(:,patmo_idx_32SO3) &
        ! + krate(:,290)*n(:,patmo_idx_OH)*n(:,patmo_idx_33SO3) &
        ! + krate(:,291)*n(:,patmo_idx_OH)*n(:,patmo_idx_34SO3) &
        ! + krate(:,292)*n(:,patmo_idx_OH)*n(:,patmo_idx_36SO3) &
        ! - krate(:,305)*n(:,patmo_idx_HO2)*n(:,patmo_idx_32SO2) &
        ! - krate(:,306)*n(:,patmo_idx_HO2)*n(:,patmo_idx_33SO2) &
        ! - krate(:,307)*n(:,patmo_idx_HO2)*n(:,patmo_idx_34SO2) &
        ! - krate(:,308)*n(:,patmo_idx_HO2)*n(:,patmo_idx_36SO2) &
        ! - krate(:,309)*n(:,patmo_idx_HO2)*n(:,patmo_idx_32SO3) &
        ! - krate(:,310)*n(:,patmo_idx_HO2)*n(:,patmo_idx_33SO3) &
        ! - krate(:,311)*n(:,patmo_idx_HO2)*n(:,patmo_idx_34SO3) &
        ! - krate(:,312)*n(:,patmo_idx_HO2)*n(:,patmo_idx_36SO3)

    dn(:,patmo_idx_36HSO3) = &
        - krate(:,104)*n(:,patmo_idx_36HSO3)*n(:,patmo_idx_O2) &
        + krate(:,112)*n(:,patmo_idx_36SO2)*n(:,patmo_idx_OH) &
        + krate(:,312)*n(:,patmo_idx_HO2)*n(:,patmo_idx_36SO3) &
        - krate(:,320)*n(:,patmo_idx_36HSO3)

    dn(:,patmo_idx_N) = 0d0
        ! + krate(:,122)*n(:,patmo_idx_N2) &
        ! + krate(:,122)*n(:,patmo_idx_N2) &
        ! - krate(:,330)*n(:,patmo_idx_N)*n(:,patmo_idx_N) &
        ! - krate(:,330)*n(:,patmo_idx_N)*n(:,patmo_idx_N)

!    if (end_of_run) then
!        budget(:,patmo_idx_33SO4,1) =  + krate(:,124)*n(:,patmo_idx_33SO2)
!
!        budget(:,patmo_idx_33SO4,2) =  - krate(:,332)*n(:,patmo_idx_33SO4)
!    endif
!
    dn(:,patmo_idx_33SO4) = &
        + krate(:,124)*n(:,patmo_idx_33SO2) &
        - krate(:,332)*n(:,patmo_idx_33SO4)

    dn(:,patmo_idx_33S) = &
        + krate(:,26)*n(:,patmo_idx_33CS)*n(:,patmo_idx_O) &
        - krate(:,70)*n(:,patmo_idx_33S)*n(:,patmo_idx_O2) &
        - krate(:,74)*n(:,patmo_idx_33S)*n(:,patmo_idx_O3) &
        - krate(:,78)*n(:,patmo_idx_33S)*n(:,patmo_idx_OH) &
        + krate(:,183)*n(:,patmo_idx_33CS2) &
        + krate(:,197)*n(:,patmo_idx_33COS) &
        + krate(:,203)*n(:,patmo_idx_33SO) &
        - krate(:,234)*n(:,patmo_idx_CO)*n(:,patmo_idx_33S) &
        + krate(:,278)*n(:,patmo_idx_33SO)*n(:,patmo_idx_O) &
        + krate(:,282)*n(:,patmo_idx_O2)*n(:,patmo_idx_33SO) &
        + krate(:,286)*n(:,patmo_idx_H)*n(:,patmo_idx_33SO)

    if (end_of_run) then


        budget(:,patmo_idx_34COS,1) = &
        + (0.83*(krate(:,11)*n(:,patmo_idx_34CS2)*n(:,patmo_idx_OH))) &
        + krate(:,19)*n(:,patmo_idx_34CS)*n(:,patmo_idx_O2) &
        + krate(:,23)*n(:,patmo_idx_34CS)*n(:,patmo_idx_O3) &
        + krate(:,211)*n(:,patmo_idx_CO2)*n(:,patmo_idx_34SH) &
        + krate(:,215)*n(:,patmo_idx_CO)*n(:,patmo_idx_34SO) &
        + 0.007*(krate(:,133)*n(:,patmo_idx_34CH3SCH3)*n(:,patmo_idx_OH)) &
        + 0.007*(krate(:,137)*n(:,patmo_idx_34CH3SCH3)*n(:,patmo_idx_OH)) &
        + 0.007*(krate(:,129)*n(:,patmo_idx_34CH3SCH3)*n(:,patmo_idx_O))

        budget(:,patmo_idx_34COS,2) = &
        - krate(:,3)*n(:,patmo_idx_34COS)*n(:,patmo_idx_OH) &
        - krate(:,7)*n(:,patmo_idx_34COS)*n(:,patmo_idx_O) &
        - krate(:,219)*n(:,patmo_idx_34SH)*n(:,patmo_idx_34COS) &
        - krate(:,186)*n(:,patmo_idx_34COS) & 
        - krate(:,227)*n(:,patmo_idx_34COS)*n(:,patmo_idx_O) &
        - krate(:,231)*n(:,patmo_idx_34COS)*n(:,patmo_idx_O2) 

        budget(:,patmo_idx_34COS,9) = &
        + krate(:,186)*n(:,patmo_idx_34COS) 
    endif
   
    dn(:,patmo_idx_34COS) = &
        - krate(:,3)*n(:,patmo_idx_34COS)*n(:,patmo_idx_OH) &
        - krate(:,7)*n(:,patmo_idx_34COS)*n(:,patmo_idx_O) &
        + (0.83*(krate(:,11)*n(:,patmo_idx_34CS2)*n(:,patmo_idx_OH))) &
        + krate(:,19)*n(:,patmo_idx_34CS)*n(:,patmo_idx_O2) &
        + krate(:,23)*n(:,patmo_idx_34CS)*n(:,patmo_idx_O3) &
        - krate(:,186)*n(:,patmo_idx_34COS) &
        + krate(:,211)*n(:,patmo_idx_CO2)*n(:,patmo_idx_34SH) &
        + krate(:,215)*n(:,patmo_idx_CO)*n(:,patmo_idx_34SO) &
        - krate(:,219)*n(:,patmo_idx_34SH)*n(:,patmo_idx_34COS) &
        - krate(:,227)*n(:,patmo_idx_34COS)*n(:,patmo_idx_O) &
        - krate(:,231)*n(:,patmo_idx_34COS)*n(:,patmo_idx_O2) &
        + 0.007*(krate(:,133)*n(:,patmo_idx_34CH3SCH3)*n(:,patmo_idx_OH)) &
        + 0.007*(krate(:,137)*n(:,patmo_idx_34CH3SCH3)*n(:,patmo_idx_OH)) &
	+ 0.007*(krate(:,129)*n(:,patmo_idx_34CH3SCH3)*n(:,patmo_idx_O))

    dn(:,patmo_idx_34SO) = &
        + krate(:,7)*n(:,patmo_idx_34COS)*n(:,patmo_idx_O) &
        + krate(:,15)*n(:,patmo_idx_34CS2)*n(:,patmo_idx_O) &
        + krate(:,47)*n(:,patmo_idx_34SH)*n(:,patmo_idx_O) &
        + krate(:,51)*n(:,patmo_idx_34SH)*n(:,patmo_idx_O2) &
        - krate(:,59)*n(:,patmo_idx_34SO)*n(:,patmo_idx_O3) &
        - krate(:,63)*n(:,patmo_idx_34SO)*n(:,patmo_idx_O2) &
        - krate(:,67)*n(:,patmo_idx_34SO)*n(:,patmo_idx_OH) &
        + krate(:,71)*n(:,patmo_idx_34S)*n(:,patmo_idx_O2) &
        + krate(:,75)*n(:,patmo_idx_34S)*n(:,patmo_idx_O3) &
        + krate(:,79)*n(:,patmo_idx_34S)*n(:,patmo_idx_OH) &
        + krate(:,181)*n(:,patmo_idx_34SO2) &
        - krate(:,195)*n(:,patmo_idx_34SO) &
        - krate(:,215)*n(:,patmo_idx_CO)*n(:,patmo_idx_34SO) &
        - krate(:,223)*n(:,patmo_idx_34CS)*n(:,patmo_idx_34SO) &
        - krate(:,255)*n(:,patmo_idx_H)*n(:,patmo_idx_34SO) &
        - krate(:,259)*n(:,patmo_idx_OH)*n(:,patmo_idx_34SO) &
        + krate(:,267)*n(:,patmo_idx_34SO2)*n(:,patmo_idx_O2) &
        + krate(:,271)*n(:,patmo_idx_34SO2)*n(:,patmo_idx_O) &
        + krate(:,275)*n(:,patmo_idx_34SO2)*n(:,patmo_idx_H) &
        - krate(:,279)*n(:,patmo_idx_34SO)*n(:,patmo_idx_O) &
        - krate(:,283)*n(:,patmo_idx_O2)*n(:,patmo_idx_34SO) &
        - krate(:,287)*n(:,patmo_idx_H)*n(:,patmo_idx_34SO)

    dn(:,patmo_idx_36SO) = &
        + krate(:,8)*n(:,patmo_idx_36COS)*n(:,patmo_idx_O) &
        + krate(:,16)*n(:,patmo_idx_36CS2)*n(:,patmo_idx_O) &
        + krate(:,48)*n(:,patmo_idx_36SH)*n(:,patmo_idx_O) &
        + krate(:,52)*n(:,patmo_idx_36SH)*n(:,patmo_idx_O2) &
        - krate(:,60)*n(:,patmo_idx_36SO)*n(:,patmo_idx_O3) &
        - krate(:,64)*n(:,patmo_idx_36SO)*n(:,patmo_idx_O2) &
        - krate(:,68)*n(:,patmo_idx_36SO)*n(:,patmo_idx_OH) &
        + krate(:,72)*n(:,patmo_idx_36S)*n(:,patmo_idx_O2) &
        + krate(:,76)*n(:,patmo_idx_36S)*n(:,patmo_idx_O3) &
        + krate(:,80)*n(:,patmo_idx_36S)*n(:,patmo_idx_OH) &
        + krate(:,175)*n(:,patmo_idx_36SO2) &
        - krate(:,204)*n(:,patmo_idx_36SO) &
        - krate(:,216)*n(:,patmo_idx_CO)*n(:,patmo_idx_36SO) &
        - krate(:,224)*n(:,patmo_idx_36CS)*n(:,patmo_idx_36SO) &
        - krate(:,256)*n(:,patmo_idx_H)*n(:,patmo_idx_36SO) &
        - krate(:,260)*n(:,patmo_idx_OH)*n(:,patmo_idx_36SO) &
        + krate(:,268)*n(:,patmo_idx_36SO2)*n(:,patmo_idx_O2) &
        + krate(:,272)*n(:,patmo_idx_36SO2)*n(:,patmo_idx_O) &
        + krate(:,276)*n(:,patmo_idx_36SO2)*n(:,patmo_idx_H) &
        - krate(:,280)*n(:,patmo_idx_36SO)*n(:,patmo_idx_O) &
        - krate(:,284)*n(:,patmo_idx_O2)*n(:,patmo_idx_36SO) &
        - krate(:,288)*n(:,patmo_idx_H)*n(:,patmo_idx_36SO)

    if (end_of_run) then
            budget(:,patmo_idx_34CS2,1) = &
                      + krate(:,219)*n(:,patmo_idx_34SH)*n(:,patmo_idx_34COS) &
                   + krate(:,223)*n(:,patmo_idx_34CS)*n(:,patmo_idx_34SO)

            budget(:,patmo_idx_34CS2,2) = &
                - krate(:,11)*n(:,patmo_idx_34CS2)*n(:,patmo_idx_OH) &
                - krate(:,15)*n(:,patmo_idx_34CS2)*n(:,patmo_idx_O) &
                - krate(:,193)*n(:,patmo_idx_34CS2) 
     
     endif



    dn(:,patmo_idx_34CS2) = &
        - krate(:,11)*n(:,patmo_idx_34CS2)*n(:,patmo_idx_OH) &
        - krate(:,15)*n(:,patmo_idx_34CS2)*n(:,patmo_idx_O) &
        - krate(:,193)*n(:,patmo_idx_34CS2) &
        + krate(:,219)*n(:,patmo_idx_34SH)*n(:,patmo_idx_34COS) &
        + krate(:,223)*n(:,patmo_idx_34CS)*n(:,patmo_idx_34SO)

    dn(:,patmo_idx_CO2) = 0d0
        ! + krate(:,1)*n(:,patmo_idx_32COS)*n(:,patmo_idx_OH) &
        ! + krate(:,2)*n(:,patmo_idx_33COS)*n(:,patmo_idx_OH) &
        ! + krate(:,3)*n(:,patmo_idx_34COS)*n(:,patmo_idx_OH) &
        ! + krate(:,4)*n(:,patmo_idx_36COS)*n(:,patmo_idx_OH) &
        ! - krate(:,209)*n(:,patmo_idx_CO2)*n(:,patmo_idx_32SH) &
        ! - krate(:,210)*n(:,patmo_idx_CO2)*n(:,patmo_idx_33SH) &
        ! - krate(:,211)*n(:,patmo_idx_CO2)*n(:,patmo_idx_34SH) &
        ! - krate(:,212)*n(:,patmo_idx_CO2)*n(:,patmo_idx_36SH)

    dn(:,patmo_idx_33CS) = &
        + krate(:,14)*n(:,patmo_idx_33CS2)*n(:,patmo_idx_O) &
        - krate(:,18)*n(:,patmo_idx_33CS)*n(:,patmo_idx_O2) &
        - krate(:,22)*n(:,patmo_idx_33CS)*n(:,patmo_idx_O3) &
        - krate(:,26)*n(:,patmo_idx_33CS)*n(:,patmo_idx_O) &
        + krate(:,183)*n(:,patmo_idx_33CS2) &
        - krate(:,222)*n(:,patmo_idx_33CS)*n(:,patmo_idx_33SO) &
        + krate(:,226)*n(:,patmo_idx_33COS)*n(:,patmo_idx_O) &
        + krate(:,230)*n(:,patmo_idx_33COS)*n(:,patmo_idx_O2) &
        + krate(:,234)*n(:,patmo_idx_CO)*n(:,patmo_idx_33S)

    dn(:,patmo_idx_33SO) = &
        + krate(:,6)*n(:,patmo_idx_33COS)*n(:,patmo_idx_O) &
        + krate(:,14)*n(:,patmo_idx_33CS2)*n(:,patmo_idx_O) &
        + krate(:,46)*n(:,patmo_idx_33SH)*n(:,patmo_idx_O) &
        + krate(:,50)*n(:,patmo_idx_33SH)*n(:,patmo_idx_O2) &
        - krate(:,58)*n(:,patmo_idx_33SO)*n(:,patmo_idx_O3) &
        - krate(:,62)*n(:,patmo_idx_33SO)*n(:,patmo_idx_O2) &
        - krate(:,66)*n(:,patmo_idx_33SO)*n(:,patmo_idx_OH) &
        + krate(:,70)*n(:,patmo_idx_33S)*n(:,patmo_idx_O2) &
        + krate(:,74)*n(:,patmo_idx_33S)*n(:,patmo_idx_O3) &
        + krate(:,78)*n(:,patmo_idx_33S)*n(:,patmo_idx_OH) &
        + krate(:,199)*n(:,patmo_idx_33SO2) &
        - krate(:,203)*n(:,patmo_idx_33SO) &
        - krate(:,214)*n(:,patmo_idx_CO)*n(:,patmo_idx_33SO) &
        - krate(:,222)*n(:,patmo_idx_33CS)*n(:,patmo_idx_33SO) &
        - krate(:,254)*n(:,patmo_idx_H)*n(:,patmo_idx_33SO) &
        - krate(:,258)*n(:,patmo_idx_OH)*n(:,patmo_idx_33SO) &
        + krate(:,266)*n(:,patmo_idx_33SO2)*n(:,patmo_idx_O2) &
        + krate(:,270)*n(:,patmo_idx_33SO2)*n(:,patmo_idx_O) &
        + krate(:,274)*n(:,patmo_idx_33SO2)*n(:,patmo_idx_H) &
        - krate(:,278)*n(:,patmo_idx_33SO)*n(:,patmo_idx_O) &
        - krate(:,282)*n(:,patmo_idx_O2)*n(:,patmo_idx_33SO) &
        - krate(:,286)*n(:,patmo_idx_H)*n(:,patmo_idx_33SO)

    dn(:,patmo_idx_32SO2_3) = &
        + krate(:,151)*n(:,patmo_idx_32SO2_1) &
        + krate(:,155)*n(:,patmo_idx_32SO2_1) &
        - krate(:,159)*n(:,patmo_idx_32SO2_3) &
        - krate(:,163)*n(:,patmo_idx_32SO2_3) &
        - krate(:,167)*n(:,patmo_idx_32SO2_3)*n(:,patmo_idx_O2) &
        - krate(:,171)*n(:,patmo_idx_32SO2_3)*n(:,patmo_idx_32SO2) &
        + krate(:,190)*n(:,patmo_idx_32SO2) &
        - krate(:,359)*n(:,patmo_idx_32SO2_3)*n(:,patmo_idx_N2) &
        - krate(:,363)*n(:,patmo_idx_32SO2_3) &
        + krate(:,367)*n(:,patmo_idx_32SO2) &
        + krate(:,371)*n(:,patmo_idx_32SO2)*n(:,patmo_idx_N2) &
        + krate(:,375)*n(:,patmo_idx_32SO3)*n(:,patmo_idx_O) &
        + krate(:,379)*n(:,patmo_idx_32SO3)*n(:,patmo_idx_32SO)

    dn(:,patmo_idx_32S) = &
        + krate(:,25)*n(:,patmo_idx_32CS)*n(:,patmo_idx_O) &
        - krate(:,69)*n(:,patmo_idx_32S)*n(:,patmo_idx_O2) &
        - krate(:,73)*n(:,patmo_idx_32S)*n(:,patmo_idx_O3) &
        - krate(:,77)*n(:,patmo_idx_32S)*n(:,patmo_idx_OH) &
        + krate(:,179)*n(:,patmo_idx_32COS) &
        + krate(:,196)*n(:,patmo_idx_32SO) &
        + krate(:,202)*n(:,patmo_idx_32CS2) &
        - krate(:,233)*n(:,patmo_idx_CO)*n(:,patmo_idx_32S) &
        + krate(:,277)*n(:,patmo_idx_32SO)*n(:,patmo_idx_O) &
        + krate(:,281)*n(:,patmo_idx_O2)*n(:,patmo_idx_32SO) &
        + krate(:,285)*n(:,patmo_idx_H)*n(:,patmo_idx_32SO)
    
!    if (end_of_run) then
!        budget(:,patmo_idx_36SO2,1) = &
!        + krate(:,60)*n(:,patmo_idx_36SO)*n(:,patmo_idx_O3) &
!        + krate(:,64)*n(:,patmo_idx_36SO)*n(:,patmo_idx_O2) &
!        + krate(:,68)*n(:,patmo_idx_36SO)*n(:,patmo_idx_OH) &
!        + krate(:,92)*n(:,patmo_idx_36HSO)*n(:,patmo_idx_O2) &
!        + krate(:,100)*n(:,patmo_idx_36HSO2)*n(:,patmo_idx_O2) &
!        + krate(:,120)*n(:,patmo_idx_36H2SO4) &
!        + krate(:,130)*n(:,patmo_idx_36CH3SCH3)*n(:,patmo_idx_O) &
!        + krate(:,134)*n(:,patmo_idx_36CH3SCH3)*n(:,patmo_idx_OH) &
!        + krate(:,138)*n(:,patmo_idx_36CH3SCH3)*n(:,patmo_idx_OH) &
!        + krate(:,142)*n(:,patmo_idx_36SO2_1) &
!        + krate(:,146)*n(:,patmo_idx_36SO2_1) &
!        + krate(:,162)*n(:,patmo_idx_36SO2_3) &
!        + krate(:,166)*n(:,patmo_idx_36SO2_3) &
!        + krate(:,178)*n(:,patmo_idx_36SO3) &
!        + krate(:,292)*n(:,patmo_idx_OH)*n(:,patmo_idx_36SO3) &
!        + krate(:,296)*n(:,patmo_idx_36SO3)*n(:,patmo_idx_O2) &
!        + krate(:,316)*n(:,patmo_idx_36SO3) &
!        + krate(:,320)*n(:,patmo_idx_36HSO3) &
!        + krate(:,334)*n(:,patmo_idx_36SO4) 
!        
!        budget(:,patmo_idx_36SO2,2) = &
!        - krate(:,84)*n(:,patmo_idx_36SO2)*n(:,patmo_idx_HO2) &
!        - krate(:,88)*n(:,patmo_idx_36SO2)*n(:,patmo_idx_O3) &
!        - krate(:,108)*n(:,patmo_idx_36SO2)*n(:,patmo_idx_O) &
!        - krate(:,112)*n(:,patmo_idx_36SO2)*n(:,patmo_idx_OH) &
!        - krate(:,126)*n(:,patmo_idx_36SO2) &
!        - krate(:,175)*n(:,patmo_idx_36SO2) &
!        - krate(:,176)*n(:,patmo_idx_36SO2) &
!        - krate(:,177)*n(:,patmo_idx_36SO2) &
!        - krate(:,268)*n(:,patmo_idx_36SO2)*n(:,patmo_idx_O2) &
!        - krate(:,272)*n(:,patmo_idx_36SO2)*n(:,patmo_idx_O) &
!        - krate(:,276)*n(:,patmo_idx_36SO2)*n(:,patmo_idx_H) &
!        - krate(:,300)*n(:,patmo_idx_36SO2)*n(:,patmo_idx_OH) &
!        - krate(:,308)*n(:,patmo_idx_HO2)*n(:,patmo_idx_36SO2) &
!        - krate(:,328)*n(:,patmo_idx_36SO2)*n(:,patmo_idx_OH)*n(:,patmo_idx_OH) &
!        - krate(:,338)*n(:,patmo_idx_36SO2) &
!        - krate(:,342)*n(:,patmo_idx_36SO2) &
!        - krate(:,346)*n(:,patmo_idx_36SO2)*n(:,patmo_idx_CH4O3S) &
!        - krate(:,350)*n(:,patmo_idx_36SO2)*n(:,patmo_idx_N2) &
!        - krate(:,354)*n(:,patmo_idx_36SO2) &
!        - krate(:,370)*n(:,patmo_idx_36SO2) &
!        - krate(:,374)*n(:,patmo_idx_36SO2)*n(:,patmo_idx_N2)
!    endif
!
    dn(:,patmo_idx_36SO2) = &
        + krate(:,60)*n(:,patmo_idx_36SO)*n(:,patmo_idx_O3) &
        + krate(:,64)*n(:,patmo_idx_36SO)*n(:,patmo_idx_O2) &
        + krate(:,68)*n(:,patmo_idx_36SO)*n(:,patmo_idx_OH) &
        - krate(:,84)*n(:,patmo_idx_36SO2)*n(:,patmo_idx_HO2) &
        - krate(:,88)*n(:,patmo_idx_36SO2)*n(:,patmo_idx_O3) &
        + krate(:,92)*n(:,patmo_idx_36HSO)*n(:,patmo_idx_O2) &
        + krate(:,100)*n(:,patmo_idx_36HSO2)*n(:,patmo_idx_O2) &
        - krate(:,108)*n(:,patmo_idx_36SO2)*n(:,patmo_idx_O) &
        - krate(:,112)*n(:,patmo_idx_36SO2)*n(:,patmo_idx_OH) &
        + krate(:,120)*n(:,patmo_idx_36H2SO4) &
        - krate(:,126)*n(:,patmo_idx_36SO2) &
        + 0.99*(krate(:,130)*n(:,patmo_idx_36CH3SCH3)*n(:,patmo_idx_O)) &
        + 0.99*(krate(:,134)*n(:,patmo_idx_36CH3SCH3)*n(:,patmo_idx_OH)) &
        + 0.99*(krate(:,138)*n(:,patmo_idx_36CH3SCH3)*n(:,patmo_idx_OH)) &
        + krate(:,142)*n(:,patmo_idx_36SO2_1) &
        + krate(:,146)*n(:,patmo_idx_36SO2_1) &
        + krate(:,162)*n(:,patmo_idx_36SO2_3) &
        + krate(:,166)*n(:,patmo_idx_36SO2_3) &
        - krate(:,175)*n(:,patmo_idx_36SO2) &
        - krate(:,176)*n(:,patmo_idx_36SO2) &
        - krate(:,177)*n(:,patmo_idx_36SO2) &
        + krate(:,178)*n(:,patmo_idx_36SO3) &
        - krate(:,268)*n(:,patmo_idx_36SO2)*n(:,patmo_idx_O2) &
        - krate(:,272)*n(:,patmo_idx_36SO2)*n(:,patmo_idx_O) &
        - krate(:,276)*n(:,patmo_idx_36SO2)*n(:,patmo_idx_H) &
        + krate(:,292)*n(:,patmo_idx_OH)*n(:,patmo_idx_36SO3) &
        + krate(:,296)*n(:,patmo_idx_36SO3)*n(:,patmo_idx_O2) &
        - krate(:,300)*n(:,patmo_idx_36SO2)*n(:,patmo_idx_OH) &
        - krate(:,308)*n(:,patmo_idx_HO2)*n(:,patmo_idx_36SO2) &
        + krate(:,316)*n(:,patmo_idx_36SO3) &
        + krate(:,320)*n(:,patmo_idx_36HSO3) &
        - krate(:,328)*n(:,patmo_idx_36SO2)*n(:,patmo_idx_OH)*n(:,patmo_idx_OH) &
        + krate(:,334)*n(:,patmo_idx_36SO4) &
        - krate(:,338)*n(:,patmo_idx_36SO2) &
        - krate(:,342)*n(:,patmo_idx_36SO2) &
        - krate(:,346)*n(:,patmo_idx_36SO2)*n(:,patmo_idx_CH4O3S) &
        - krate(:,350)*n(:,patmo_idx_36SO2)*n(:,patmo_idx_N2) &
        - krate(:,354)*n(:,patmo_idx_36SO2) &
        - krate(:,370)*n(:,patmo_idx_36SO2) &
        - krate(:,374)*n(:,patmo_idx_36SO2)*n(:,patmo_idx_N2)

    dn(:,patmo_idx_32HSO) = &
        + krate(:,41)*n(:,patmo_idx_32H2S)*n(:,patmo_idx_HO2) &
        + krate(:,53)*n(:,patmo_idx_32SH)*n(:,patmo_idx_O3) &
        - krate(:,89)*n(:,patmo_idx_32HSO)*n(:,patmo_idx_O2) &
        - krate(:,93)*n(:,patmo_idx_32HSO)*n(:,patmo_idx_O3) &
        - krate(:,249)*n(:,patmo_idx_H2O)*n(:,patmo_idx_32HSO) &
        - krate(:,261)*n(:,patmo_idx_32HSO)*n(:,patmo_idx_O2) &
        + krate(:,297)*n(:,patmo_idx_32SO2)*n(:,patmo_idx_OH) &
        + krate(:,301)*n(:,patmo_idx_O2)*n(:,patmo_idx_O2)*n(:,patmo_idx_32SH)

    dn(:,patmo_idx_34S) = &
        + krate(:,27)*n(:,patmo_idx_34CS)*n(:,patmo_idx_O) &
        - krate(:,71)*n(:,patmo_idx_34S)*n(:,patmo_idx_O2) &
        - krate(:,75)*n(:,patmo_idx_34S)*n(:,patmo_idx_O3) &
        - krate(:,79)*n(:,patmo_idx_34S)*n(:,patmo_idx_OH) &
        + krate(:,186)*n(:,patmo_idx_34COS) &
        + krate(:,193)*n(:,patmo_idx_34CS2) &
        + krate(:,195)*n(:,patmo_idx_34SO) &
        - krate(:,235)*n(:,patmo_idx_CO)*n(:,patmo_idx_34S) &
        + krate(:,279)*n(:,patmo_idx_34SO)*n(:,patmo_idx_O) &
        + krate(:,283)*n(:,patmo_idx_O2)*n(:,patmo_idx_34SO) &
        + krate(:,287)*n(:,patmo_idx_H)*n(:,patmo_idx_34SO)

!H2SO4 chemistry budget

    if (end_of_run) then 
        budget(:,patmo_idx_34H2SO4,1) = &
        + krate(:,115)*n(:,patmo_idx_34SO3)*n(:,patmo_idx_H2O) &
        + krate(:,327)*n(:,patmo_idx_34SO2)*n(:,patmo_idx_OH)*n(:,patmo_idx_OH)
        
        budget(:,patmo_idx_34H2SO4,8) = &
        + krate(:,327)*n(:,patmo_idx_34SO2)*n(:,patmo_idx_OH)*n(:,patmo_idx_OH)

        budget(:,patmo_idx_34H2SO4,9) = &
        + krate(:,115)*n(:,patmo_idx_34SO3)*n(:,patmo_idx_H2O)

        budget(:,patmo_idx_34H2SO4,2) = &
        + krate(:,119)*n(:,patmo_idx_34H2SO4) &
        + krate(:,323)*n(:,patmo_idx_34H2SO4) 

        budget(:,patmo_idx_34H2SO4,10) = &
        + krate(:,119)*n(:,patmo_idx_34H2SO4) 

    endif

    dn(:,patmo_idx_34H2SO4) = &
        + krate(:,115)*n(:,patmo_idx_34SO3)*n(:,patmo_idx_H2O) &
        - krate(:,119)*n(:,patmo_idx_34H2SO4) &
        - krate(:,323)*n(:,patmo_idx_34H2SO4) &
        + krate(:,327)*n(:,patmo_idx_34SO2)*n(:,patmo_idx_OH)*n(:,patmo_idx_OH)

    dn(:,patmo_idx_H2O) = 0d0
        ! + krate(:,29)*n(:,patmo_idx_32H2S)*n(:,patmo_idx_OH) &
        ! + krate(:,30)*n(:,patmo_idx_33H2S)*n(:,patmo_idx_OH) &
        ! + krate(:,31)*n(:,patmo_idx_34H2S)*n(:,patmo_idx_OH) &
        ! + krate(:,32)*n(:,patmo_idx_36H2S)*n(:,patmo_idx_OH) &
        ! + krate(:,41)*n(:,patmo_idx_32H2S)*n(:,patmo_idx_HO2) &
        ! + krate(:,42)*n(:,patmo_idx_33H2S)*n(:,patmo_idx_HO2) &
        ! + krate(:,43)*n(:,patmo_idx_34H2S)*n(:,patmo_idx_HO2) &
        ! + krate(:,44)*n(:,patmo_idx_36H2S)*n(:,patmo_idx_HO2) &
        ! - krate(:,113)*n(:,patmo_idx_32SO3)*n(:,patmo_idx_H2O) &
        ! - krate(:,114)*n(:,patmo_idx_33SO3)*n(:,patmo_idx_H2O) &
        ! - krate(:,115)*n(:,patmo_idx_34SO3)*n(:,patmo_idx_H2O) &
        ! - krate(:,116)*n(:,patmo_idx_36SO3)*n(:,patmo_idx_H2O) &
        ! - krate(:,237)*n(:,patmo_idx_H2O)*n(:,patmo_idx_32SH) &
        ! - krate(:,238)*n(:,patmo_idx_H2O)*n(:,patmo_idx_33SH) &
        ! - krate(:,239)*n(:,patmo_idx_H2O)*n(:,patmo_idx_34SH) &
        ! - krate(:,240)*n(:,patmo_idx_H2O)*n(:,patmo_idx_36SH) &
        ! - krate(:,249)*n(:,patmo_idx_H2O)*n(:,patmo_idx_32HSO) &
        ! - krate(:,250)*n(:,patmo_idx_H2O)*n(:,patmo_idx_33HSO) &
        ! - krate(:,251)*n(:,patmo_idx_H2O)*n(:,patmo_idx_34HSO) &
        ! - krate(:,252)*n(:,patmo_idx_H2O)*n(:,patmo_idx_36HSO) &
        ! + krate(:,321)*n(:,patmo_idx_32H2SO4) &
        ! + krate(:,322)*n(:,patmo_idx_33H2SO4) &
        ! + krate(:,323)*n(:,patmo_idx_34H2SO4) &
        ! + krate(:,324)*n(:,patmo_idx_36H2SO4)

    dn(:,patmo_idx_36H2S) = &
        - krate(:,32)*n(:,patmo_idx_36H2S)*n(:,patmo_idx_OH) &
        - krate(:,36)*n(:,patmo_idx_36H2S)*n(:,patmo_idx_O) &
        - krate(:,40)*n(:,patmo_idx_36H2S)*n(:,patmo_idx_H) &
        - krate(:,44)*n(:,patmo_idx_36H2S)*n(:,patmo_idx_HO2) &
        - krate(:,187)*n(:,patmo_idx_36H2S) &
        + krate(:,240)*n(:,patmo_idx_H2O)*n(:,patmo_idx_36SH) &
        + krate(:,244)*n(:,patmo_idx_OH)*n(:,patmo_idx_36SH) &
        + krate(:,248)*n(:,patmo_idx_H2)*n(:,patmo_idx_36SH) &
        + krate(:,252)*n(:,patmo_idx_H2O)*n(:,patmo_idx_36HSO)

    dn(:,patmo_idx_34CS) = &
        + krate(:,15)*n(:,patmo_idx_34CS2)*n(:,patmo_idx_O) &
        - krate(:,19)*n(:,patmo_idx_34CS)*n(:,patmo_idx_O2) &
        - krate(:,23)*n(:,patmo_idx_34CS)*n(:,patmo_idx_O3) &
        - krate(:,27)*n(:,patmo_idx_34CS)*n(:,patmo_idx_O) &
        + krate(:,193)*n(:,patmo_idx_34CS2) &
        - krate(:,223)*n(:,patmo_idx_34CS)*n(:,patmo_idx_34SO) &
        + krate(:,227)*n(:,patmo_idx_34COS)*n(:,patmo_idx_O) &
        + krate(:,231)*n(:,patmo_idx_34COS)*n(:,patmo_idx_O2) &
        + krate(:,235)*n(:,patmo_idx_CO)*n(:,patmo_idx_34S)

    dn(:,patmo_idx_33HSO2) = &
        - krate(:,98)*n(:,patmo_idx_33HSO2)*n(:,patmo_idx_O2) &
        + krate(:,306)*n(:,patmo_idx_HO2)*n(:,patmo_idx_33SO2)

    dn(:,patmo_idx_34SO2_3) = &
        + krate(:,153)*n(:,patmo_idx_34SO2_1) &
        + krate(:,157)*n(:,patmo_idx_34SO2_1) &
        - krate(:,161)*n(:,patmo_idx_34SO2_3) &
        - krate(:,165)*n(:,patmo_idx_34SO2_3) &
        - krate(:,169)*n(:,patmo_idx_34SO2_3)*n(:,patmo_idx_O2) &
        - krate(:,173)*n(:,patmo_idx_34SO2_3)*n(:,patmo_idx_32SO2) &
        + krate(:,180)*n(:,patmo_idx_34SO2) &
        - krate(:,361)*n(:,patmo_idx_34SO2_3)*n(:,patmo_idx_N2) &
        - krate(:,365)*n(:,patmo_idx_34SO2_3) &
        + krate(:,369)*n(:,patmo_idx_34SO2) &
        + krate(:,373)*n(:,patmo_idx_34SO2)*n(:,patmo_idx_N2) &
        + krate(:,377)*n(:,patmo_idx_34SO3)*n(:,patmo_idx_O) &
        + krate(:,381)*n(:,patmo_idx_34SO3)*n(:,patmo_idx_32SO)

    dn(:,patmo_idx_33SO2_3) = &
        + krate(:,152)*n(:,patmo_idx_33SO2_1) &
        + krate(:,156)*n(:,patmo_idx_33SO2_1) &
        - krate(:,160)*n(:,patmo_idx_33SO2_3) &
        - krate(:,164)*n(:,patmo_idx_33SO2_3) &
        - krate(:,168)*n(:,patmo_idx_33SO2_3)*n(:,patmo_idx_O2) &
        - krate(:,172)*n(:,patmo_idx_33SO2_3)*n(:,patmo_idx_32SO2) &
        + krate(:,198)*n(:,patmo_idx_33SO2) &
        - krate(:,360)*n(:,patmo_idx_33SO2_3)*n(:,patmo_idx_N2) &
        - krate(:,364)*n(:,patmo_idx_33SO2_3) &
        + krate(:,368)*n(:,patmo_idx_33SO2) &
        + krate(:,372)*n(:,patmo_idx_33SO2)*n(:,patmo_idx_N2) &
        + krate(:,376)*n(:,patmo_idx_33SO3)*n(:,patmo_idx_O) &
        + krate(:,380)*n(:,patmo_idx_33SO3)*n(:,patmo_idx_32SO)

    dn(:,patmo_idx_33CH3SCH3) = &
        - krate(:,128)*n(:,patmo_idx_33CH3SCH3)*n(:,patmo_idx_O) &
        - krate(:,132)*n(:,patmo_idx_33CH3SCH3)*n(:,patmo_idx_OH) &
        - krate(:,136)*n(:,patmo_idx_33CH3SCH3)*n(:,patmo_idx_OH) &
        + krate(:,336)*n(:,patmo_idx_33SO2) &
        + krate(:,340)*n(:,patmo_idx_33SO2) &
        + krate(:,344)*n(:,patmo_idx_33SO2)*n(:,patmo_idx_CH4O3S) 
       ! - (7.940d-6)*n(:,patmo_idx_33CH3SCH3)*n(:,patmo_idx_OH) 

    dn(:,patmo_idx_32SO2_1) = &
        - krate(:,139)*n(:,patmo_idx_32SO2_1) &
        - krate(:,143)*n(:,patmo_idx_32SO2_1) &
        - krate(:,147)*n(:,patmo_idx_32SO2_1)*n(:,patmo_idx_32SO2) &
        - krate(:,151)*n(:,patmo_idx_32SO2_1) &
        - krate(:,155)*n(:,patmo_idx_32SO2_1) &
        + krate(:,192)*n(:,patmo_idx_32SO2) &
        + krate(:,347)*n(:,patmo_idx_32SO2)*n(:,patmo_idx_N2) &
        + krate(:,351)*n(:,patmo_idx_32SO2) &
        + krate(:,355)*n(:,patmo_idx_32SO3)*n(:,patmo_idx_32SO) &
        + krate(:,359)*n(:,patmo_idx_32SO2_3)*n(:,patmo_idx_N2) &
        + krate(:,363)*n(:,patmo_idx_32SO2_3)

    dn(:,patmo_idx_CO) = 0d0
        ! + krate(:,5)*n(:,patmo_idx_32COS)*n(:,patmo_idx_O) &
        ! + krate(:,6)*n(:,patmo_idx_33COS)*n(:,patmo_idx_O) &
        ! + krate(:,7)*n(:,patmo_idx_34COS)*n(:,patmo_idx_O) &
        ! + krate(:,8)*n(:,patmo_idx_36COS)*n(:,patmo_idx_O) &
        ! + krate(:,25)*n(:,patmo_idx_32CS)*n(:,patmo_idx_O) &
        ! + krate(:,26)*n(:,patmo_idx_33CS)*n(:,patmo_idx_O) &
        ! + krate(:,27)*n(:,patmo_idx_34CS)*n(:,patmo_idx_O) &
        ! + krate(:,28)*n(:,patmo_idx_36CS)*n(:,patmo_idx_O) &
        ! + krate(:,179)*n(:,patmo_idx_32COS) &
        ! + krate(:,186)*n(:,patmo_idx_34COS) &
        ! + krate(:,197)*n(:,patmo_idx_33COS) &
        ! + krate(:,207)*n(:,patmo_idx_36COS) &
        ! - krate(:,213)*n(:,patmo_idx_CO)*n(:,patmo_idx_32SO) &
        ! - krate(:,214)*n(:,patmo_idx_CO)*n(:,patmo_idx_33SO) &
        ! - krate(:,215)*n(:,patmo_idx_CO)*n(:,patmo_idx_34SO) &
        ! - krate(:,216)*n(:,patmo_idx_CO)*n(:,patmo_idx_36SO) &
        ! - krate(:,233)*n(:,patmo_idx_CO)*n(:,patmo_idx_32S) &
        ! - krate(:,234)*n(:,patmo_idx_CO)*n(:,patmo_idx_33S) &
        ! - krate(:,235)*n(:,patmo_idx_CO)*n(:,patmo_idx_34S) &
        ! - krate(:,236)*n(:,patmo_idx_CO)*n(:,patmo_idx_36S)

    if(end_of_run)then
            budget(:,patmo_idx_32H2S,1) = &
                     + krate(:,237)*n(:,patmo_idx_H2O)*n(:,patmo_idx_32SH) &
                     + krate(:,241)*n(:,patmo_idx_OH)*n(:,patmo_idx_32SH) &
                     + krate(:,245)*n(:,patmo_idx_H2)*n(:,patmo_idx_32SH) &
                     + krate(:,249)*n(:,patmo_idx_H2O)*n(:,patmo_idx_32HSO)

            budget(:,patmo_idx_32H2S,2) = &
                     - krate(:,29)*n(:,patmo_idx_32H2S)*n(:,patmo_idx_OH) &
                     - krate(:,33)*n(:,patmo_idx_32H2S)*n(:,patmo_idx_O) &
                     - krate(:,37)*n(:,patmo_idx_32H2S)*n(:,patmo_idx_H) &
                     - krate(:,41)*n(:,patmo_idx_32H2S)*n(:,patmo_idx_HO2) &
     endif



    dn(:,patmo_idx_32H2S) = &
        - krate(:,29)*n(:,patmo_idx_32H2S)*n(:,patmo_idx_OH) &
        - krate(:,33)*n(:,patmo_idx_32H2S)*n(:,patmo_idx_O) &
        - krate(:,37)*n(:,patmo_idx_32H2S)*n(:,patmo_idx_H) &
        - krate(:,41)*n(:,patmo_idx_32H2S)*n(:,patmo_idx_HO2) &
        - krate(:,206)*n(:,patmo_idx_32H2S) &
        + krate(:,237)*n(:,patmo_idx_H2O)*n(:,patmo_idx_32SH) &
        + krate(:,241)*n(:,patmo_idx_OH)*n(:,patmo_idx_32SH) &
        + krate(:,245)*n(:,patmo_idx_H2)*n(:,patmo_idx_32SH) &
        + krate(:,249)*n(:,patmo_idx_H2O)*n(:,patmo_idx_32HSO)

    dn(:,patmo_idx_33H2SO4) = &
        + krate(:,114)*n(:,patmo_idx_33SO3)*n(:,patmo_idx_H2O) &
        - krate(:,118)*n(:,patmo_idx_33H2SO4) &
        - krate(:,322)*n(:,patmo_idx_33H2SO4) &
        + krate(:,326)*n(:,patmo_idx_33SO2)*n(:,patmo_idx_OH)*n(:,patmo_idx_OH)

!    if (end_of_run) then
!        budget(:,patmo_idx_36SO4,1) = + krate(:,126)*n(:,patmo_idx_36SO2)
!        
!        budget(:,patmo_idx_36SO4,2) = - krate(:,334)*n(:,patmo_idx_36SO4)
!    endif
!
    dn(:,patmo_idx_36SO4) = &
        + krate(:,126)*n(:,patmo_idx_36SO2) &
        - krate(:,334)*n(:,patmo_idx_36SO4)

    dn(:,patmo_idx_O2) = 0d0
        ! - krate(:,17)*n(:,patmo_idx_32CS)*n(:,patmo_idx_O2) &
        ! - krate(:,18)*n(:,patmo_idx_33CS)*n(:,patmo_idx_O2) &
        ! - krate(:,19)*n(:,patmo_idx_34CS)*n(:,patmo_idx_O2) &
        ! - krate(:,20)*n(:,patmo_idx_36CS)*n(:,patmo_idx_O2) &
        ! + krate(:,21)*n(:,patmo_idx_32CS)*n(:,patmo_idx_O3) &
        ! + krate(:,22)*n(:,patmo_idx_33CS)*n(:,patmo_idx_O3) &
        ! + krate(:,23)*n(:,patmo_idx_34CS)*n(:,patmo_idx_O3) &
        ! + krate(:,24)*n(:,patmo_idx_36CS)*n(:,patmo_idx_O3) &
        ! - krate(:,49)*n(:,patmo_idx_32SH)*n(:,patmo_idx_O2) &
        ! - krate(:,50)*n(:,patmo_idx_33SH)*n(:,patmo_idx_O2) &
        ! - krate(:,51)*n(:,patmo_idx_34SH)*n(:,patmo_idx_O2) &
        ! - krate(:,52)*n(:,patmo_idx_36SH)*n(:,patmo_idx_O2) &
        ! + krate(:,53)*n(:,patmo_idx_32SH)*n(:,patmo_idx_O3) &
        ! + krate(:,54)*n(:,patmo_idx_33SH)*n(:,patmo_idx_O3) &
        ! + krate(:,55)*n(:,patmo_idx_34SH)*n(:,patmo_idx_O3) &
        ! + krate(:,56)*n(:,patmo_idx_36SH)*n(:,patmo_idx_O3) &
        ! + krate(:,57)*n(:,patmo_idx_32SO)*n(:,patmo_idx_O3) &
        ! + krate(:,58)*n(:,patmo_idx_33SO)*n(:,patmo_idx_O3) &
        ! + krate(:,59)*n(:,patmo_idx_34SO)*n(:,patmo_idx_O3) &
        ! + krate(:,60)*n(:,patmo_idx_36SO)*n(:,patmo_idx_O3) &
        ! - krate(:,61)*n(:,patmo_idx_32SO)*n(:,patmo_idx_O2) &
        ! - krate(:,62)*n(:,patmo_idx_33SO)*n(:,patmo_idx_O2) &
        ! - krate(:,63)*n(:,patmo_idx_34SO)*n(:,patmo_idx_O2) &
        ! - krate(:,64)*n(:,patmo_idx_36SO)*n(:,patmo_idx_O2) &
        ! - krate(:,69)*n(:,patmo_idx_32S)*n(:,patmo_idx_O2) &
        ! - krate(:,70)*n(:,patmo_idx_33S)*n(:,patmo_idx_O2) &
        ! - krate(:,71)*n(:,patmo_idx_34S)*n(:,patmo_idx_O2) &
        ! - krate(:,72)*n(:,patmo_idx_36S)*n(:,patmo_idx_O2) &
        ! + krate(:,73)*n(:,patmo_idx_32S)*n(:,patmo_idx_O3) &
        ! + krate(:,74)*n(:,patmo_idx_33S)*n(:,patmo_idx_O3) &
        ! + krate(:,75)*n(:,patmo_idx_34S)*n(:,patmo_idx_O3) &
        ! + krate(:,76)*n(:,patmo_idx_36S)*n(:,patmo_idx_O3) &
        ! + krate(:,85)*n(:,patmo_idx_32SO2)*n(:,patmo_idx_O3) &
        ! + krate(:,86)*n(:,patmo_idx_33SO2)*n(:,patmo_idx_O3) &
        ! + krate(:,87)*n(:,patmo_idx_34SO2)*n(:,patmo_idx_O3) &
        ! + krate(:,88)*n(:,patmo_idx_36SO2)*n(:,patmo_idx_O3) &
        ! - krate(:,89)*n(:,patmo_idx_32HSO)*n(:,patmo_idx_O2) &
        ! - krate(:,90)*n(:,patmo_idx_33HSO)*n(:,patmo_idx_O2) &
        ! - krate(:,91)*n(:,patmo_idx_34HSO)*n(:,patmo_idx_O2) &
        ! - krate(:,92)*n(:,patmo_idx_36HSO)*n(:,patmo_idx_O2) &
        ! + krate(:,93)*n(:,patmo_idx_32HSO)*n(:,patmo_idx_O3) &
        ! + krate(:,93)*n(:,patmo_idx_32HSO)*n(:,patmo_idx_O3) &
        ! + krate(:,94)*n(:,patmo_idx_33HSO)*n(:,patmo_idx_O3) &
        ! + krate(:,94)*n(:,patmo_idx_33HSO)*n(:,patmo_idx_O3) &
        ! + krate(:,95)*n(:,patmo_idx_34HSO)*n(:,patmo_idx_O3) &
        ! + krate(:,95)*n(:,patmo_idx_34HSO)*n(:,patmo_idx_O3) &
        ! + krate(:,96)*n(:,patmo_idx_36HSO)*n(:,patmo_idx_O3) &
        ! + krate(:,96)*n(:,patmo_idx_36HSO)*n(:,patmo_idx_O3) &
        ! - krate(:,97)*n(:,patmo_idx_32HSO2)*n(:,patmo_idx_O2) &
        ! - krate(:,98)*n(:,patmo_idx_33HSO2)*n(:,patmo_idx_O2) &
        ! - krate(:,99)*n(:,patmo_idx_34HSO2)*n(:,patmo_idx_O2) &
        ! - krate(:,100)*n(:,patmo_idx_36HSO2)*n(:,patmo_idx_O2) &
        ! - krate(:,101)*n(:,patmo_idx_32HSO3)*n(:,patmo_idx_O2) &
        ! - krate(:,102)*n(:,patmo_idx_33HSO3)*n(:,patmo_idx_O2) &
        ! - krate(:,103)*n(:,patmo_idx_34HSO3)*n(:,patmo_idx_O2) &
        ! - krate(:,104)*n(:,patmo_idx_36HSO3)*n(:,patmo_idx_O2) &
        ! - krate(:,121)*n(:,patmo_idx_O)*n(:,patmo_idx_O2) &
        ! - krate(:,167)*n(:,patmo_idx_32SO2_3)*n(:,patmo_idx_O2) &
        ! - krate(:,168)*n(:,patmo_idx_33SO2_3)*n(:,patmo_idx_O2) &
        ! - krate(:,169)*n(:,patmo_idx_34SO2_3)*n(:,patmo_idx_O2) &
        ! - krate(:,170)*n(:,patmo_idx_36SO2_3)*n(:,patmo_idx_O2) &
        ! + krate(:,184)*n(:,patmo_idx_O3) &
        ! - krate(:,185)*n(:,patmo_idx_O2) &
        ! + krate(:,225)*n(:,patmo_idx_32COS)*n(:,patmo_idx_O) &
        ! + krate(:,226)*n(:,patmo_idx_33COS)*n(:,patmo_idx_O) &
        ! + krate(:,227)*n(:,patmo_idx_34COS)*n(:,patmo_idx_O) &
        ! + krate(:,228)*n(:,patmo_idx_36COS)*n(:,patmo_idx_O) &
        ! - krate(:,229)*n(:,patmo_idx_32COS)*n(:,patmo_idx_O2) &
        ! - krate(:,230)*n(:,patmo_idx_33COS)*n(:,patmo_idx_O2) &
        ! - krate(:,231)*n(:,patmo_idx_34COS)*n(:,patmo_idx_O2) &
        ! - krate(:,232)*n(:,patmo_idx_36COS)*n(:,patmo_idx_O2) &
        ! + krate(:,257)*n(:,patmo_idx_OH)*n(:,patmo_idx_32SO) &
        ! + krate(:,258)*n(:,patmo_idx_OH)*n(:,patmo_idx_33SO) &
        ! + krate(:,259)*n(:,patmo_idx_OH)*n(:,patmo_idx_34SO) &
        ! + krate(:,260)*n(:,patmo_idx_OH)*n(:,patmo_idx_36SO) &
        ! - krate(:,261)*n(:,patmo_idx_32HSO)*n(:,patmo_idx_O2) &
        ! - krate(:,262)*n(:,patmo_idx_33HSO)*n(:,patmo_idx_O2) &
        ! - krate(:,263)*n(:,patmo_idx_34HSO)*n(:,patmo_idx_O2) &
        ! - krate(:,264)*n(:,patmo_idx_36HSO)*n(:,patmo_idx_O2) &
        ! - krate(:,265)*n(:,patmo_idx_32SO2)*n(:,patmo_idx_O2) &
        ! - krate(:,266)*n(:,patmo_idx_33SO2)*n(:,patmo_idx_O2) &
        ! - krate(:,267)*n(:,patmo_idx_34SO2)*n(:,patmo_idx_O2) &
        ! - krate(:,268)*n(:,patmo_idx_36SO2)*n(:,patmo_idx_O2) &
        ! + krate(:,269)*n(:,patmo_idx_32SO2)*n(:,patmo_idx_O) &
        ! + krate(:,270)*n(:,patmo_idx_33SO2)*n(:,patmo_idx_O) &
        ! + krate(:,271)*n(:,patmo_idx_34SO2)*n(:,patmo_idx_O) &
        ! + krate(:,272)*n(:,patmo_idx_36SO2)*n(:,patmo_idx_O) &
        ! + krate(:,277)*n(:,patmo_idx_32SO)*n(:,patmo_idx_O) &
        ! + krate(:,278)*n(:,patmo_idx_33SO)*n(:,patmo_idx_O) &
        ! + krate(:,279)*n(:,patmo_idx_34SO)*n(:,patmo_idx_O) &
        ! + krate(:,280)*n(:,patmo_idx_36SO)*n(:,patmo_idx_O) &
        ! - krate(:,281)*n(:,patmo_idx_O2)*n(:,patmo_idx_32SO) &
        ! - krate(:,282)*n(:,patmo_idx_O2)*n(:,patmo_idx_33SO) &
        ! - krate(:,283)*n(:,patmo_idx_O2)*n(:,patmo_idx_34SO) &
        ! - krate(:,284)*n(:,patmo_idx_O2)*n(:,patmo_idx_36SO) &
        ! - krate(:,293)*n(:,patmo_idx_32SO3)*n(:,patmo_idx_O2) &
        ! - krate(:,294)*n(:,patmo_idx_33SO3)*n(:,patmo_idx_O2) &
        ! - krate(:,295)*n(:,patmo_idx_34SO3)*n(:,patmo_idx_O2) &
        ! - krate(:,296)*n(:,patmo_idx_36SO3)*n(:,patmo_idx_O2) &
        ! + krate(:,297)*n(:,patmo_idx_32SO2)*n(:,patmo_idx_OH) &
        ! + krate(:,298)*n(:,patmo_idx_33SO2)*n(:,patmo_idx_OH) &
        ! + krate(:,299)*n(:,patmo_idx_34SO2)*n(:,patmo_idx_OH) &
        ! + krate(:,300)*n(:,patmo_idx_36SO2)*n(:,patmo_idx_OH) &
        ! - krate(:,301)*n(:,patmo_idx_O2)*n(:,patmo_idx_O2)*n(:,patmo_idx_32SH) &
        ! - krate(:,301)*n(:,patmo_idx_O2)*n(:,patmo_idx_O2)*n(:,patmo_idx_32SH) &
        ! - krate(:,302)*n(:,patmo_idx_O2)*n(:,patmo_idx_O2)*n(:,patmo_idx_33SH) &
        ! - krate(:,302)*n(:,patmo_idx_O2)*n(:,patmo_idx_O2)*n(:,patmo_idx_33SH) &
        ! - krate(:,303)*n(:,patmo_idx_O2)*n(:,patmo_idx_O2)*n(:,patmo_idx_34SH) &
        ! - krate(:,303)*n(:,patmo_idx_O2)*n(:,patmo_idx_O2)*n(:,patmo_idx_34SH) &
        ! - krate(:,304)*n(:,patmo_idx_O2)*n(:,patmo_idx_O2)*n(:,patmo_idx_36SH) &
        ! - krate(:,304)*n(:,patmo_idx_O2)*n(:,patmo_idx_O2)*n(:,patmo_idx_36SH) &
        ! + krate(:,305)*n(:,patmo_idx_HO2)*n(:,patmo_idx_32SO2) &
        ! + krate(:,306)*n(:,patmo_idx_HO2)*n(:,patmo_idx_33SO2) &
        ! + krate(:,307)*n(:,patmo_idx_HO2)*n(:,patmo_idx_34SO2) &
        ! + krate(:,308)*n(:,patmo_idx_HO2)*n(:,patmo_idx_36SO2) &
        ! + krate(:,309)*n(:,patmo_idx_HO2)*n(:,patmo_idx_32SO3) &
        ! + krate(:,310)*n(:,patmo_idx_HO2)*n(:,patmo_idx_33SO3) &
        ! + krate(:,311)*n(:,patmo_idx_HO2)*n(:,patmo_idx_34SO3) &
        ! + krate(:,312)*n(:,patmo_idx_HO2)*n(:,patmo_idx_36SO3) &
        ! + krate(:,329)*n(:,patmo_idx_O3) &
        ! + krate(:,375)*n(:,patmo_idx_32SO3)*n(:,patmo_idx_O) &
        ! + krate(:,376)*n(:,patmo_idx_33SO3)*n(:,patmo_idx_O) &
        ! + krate(:,377)*n(:,patmo_idx_34SO3)*n(:,patmo_idx_O) &
        ! + krate(:,378)*n(:,patmo_idx_36SO3)*n(:,patmo_idx_O)

    dn(:,patmo_idx_36S) = &
        + krate(:,28)*n(:,patmo_idx_36CS)*n(:,patmo_idx_O) &
        - krate(:,72)*n(:,patmo_idx_36S)*n(:,patmo_idx_O2) &
        - krate(:,76)*n(:,patmo_idx_36S)*n(:,patmo_idx_O3) &
        - krate(:,80)*n(:,patmo_idx_36S)*n(:,patmo_idx_OH) &
        + krate(:,204)*n(:,patmo_idx_36SO) &
        + krate(:,205)*n(:,patmo_idx_36CS2) &
        + krate(:,207)*n(:,patmo_idx_36COS) &
        - krate(:,236)*n(:,patmo_idx_CO)*n(:,patmo_idx_36S) &
        + krate(:,280)*n(:,patmo_idx_36SO)*n(:,patmo_idx_O) &
        + krate(:,284)*n(:,patmo_idx_O2)*n(:,patmo_idx_36SO) &
        + krate(:,288)*n(:,patmo_idx_H)*n(:,patmo_idx_36SO)

    dn(:,patmo_idx_36H2SO4) = &
        + krate(:,116)*n(:,patmo_idx_36SO3)*n(:,patmo_idx_H2O) &
        - krate(:,120)*n(:,patmo_idx_36H2SO4) &
        - krate(:,324)*n(:,patmo_idx_36H2SO4) &
        + krate(:,328)*n(:,patmo_idx_36SO2)*n(:,patmo_idx_OH)*n(:,patmo_idx_OH)

    dn(:,patmo_idx_32SO) = &
        + krate(:,5)*n(:,patmo_idx_32COS)*n(:,patmo_idx_O) &
        + krate(:,13)*n(:,patmo_idx_32CS2)*n(:,patmo_idx_O) &
        + krate(:,45)*n(:,patmo_idx_32SH)*n(:,patmo_idx_O) &
        + krate(:,49)*n(:,patmo_idx_32SH)*n(:,patmo_idx_O2) &
        - krate(:,57)*n(:,patmo_idx_32SO)*n(:,patmo_idx_O3) &
        - krate(:,61)*n(:,patmo_idx_32SO)*n(:,patmo_idx_O2) &
        - krate(:,65)*n(:,patmo_idx_32SO)*n(:,patmo_idx_OH) &
        + krate(:,69)*n(:,patmo_idx_32S)*n(:,patmo_idx_O2) &
        + krate(:,73)*n(:,patmo_idx_32S)*n(:,patmo_idx_O3) &
        + krate(:,77)*n(:,patmo_idx_32S)*n(:,patmo_idx_OH) &
        + krate(:,147)*n(:,patmo_idx_32SO2_1)*n(:,patmo_idx_32SO2) &
        + krate(:,148)*n(:,patmo_idx_33SO2_1)*n(:,patmo_idx_32SO2) &
        + krate(:,149)*n(:,patmo_idx_34SO2_1)*n(:,patmo_idx_32SO2) &
        + krate(:,150)*n(:,patmo_idx_36SO2_1)*n(:,patmo_idx_32SO2) &
        + krate(:,171)*n(:,patmo_idx_32SO2_3)*n(:,patmo_idx_32SO2) &
        + krate(:,172)*n(:,patmo_idx_33SO2_3)*n(:,patmo_idx_32SO2) &
        + krate(:,173)*n(:,patmo_idx_34SO2_3)*n(:,patmo_idx_32SO2) &
        + krate(:,174)*n(:,patmo_idx_36SO2_3)*n(:,patmo_idx_32SO2) &
        + krate(:,191)*n(:,patmo_idx_32SO2) &
        - krate(:,196)*n(:,patmo_idx_32SO) &
        - krate(:,213)*n(:,patmo_idx_CO)*n(:,patmo_idx_32SO) &
        - krate(:,221)*n(:,patmo_idx_32CS)*n(:,patmo_idx_32SO) &
        - krate(:,253)*n(:,patmo_idx_H)*n(:,patmo_idx_32SO) &
        - krate(:,257)*n(:,patmo_idx_OH)*n(:,patmo_idx_32SO) &
        + krate(:,265)*n(:,patmo_idx_32SO2)*n(:,patmo_idx_O2) &
        + krate(:,269)*n(:,patmo_idx_32SO2)*n(:,patmo_idx_O) &
        + krate(:,273)*n(:,patmo_idx_32SO2)*n(:,patmo_idx_H) &
        - krate(:,277)*n(:,patmo_idx_32SO)*n(:,patmo_idx_O) &
        - krate(:,281)*n(:,patmo_idx_O2)*n(:,patmo_idx_32SO) &
        - krate(:,285)*n(:,patmo_idx_H)*n(:,patmo_idx_32SO) &
        - krate(:,355)*n(:,patmo_idx_32SO3)*n(:,patmo_idx_32SO) &
        - krate(:,356)*n(:,patmo_idx_33SO3)*n(:,patmo_idx_32SO) &
        - krate(:,357)*n(:,patmo_idx_34SO3)*n(:,patmo_idx_32SO) &
        - krate(:,358)*n(:,patmo_idx_36SO3)*n(:,patmo_idx_32SO) &
        - krate(:,379)*n(:,patmo_idx_32SO3)*n(:,patmo_idx_32SO) &
        - krate(:,380)*n(:,patmo_idx_33SO3)*n(:,patmo_idx_32SO) &
        - krate(:,381)*n(:,patmo_idx_34SO3)*n(:,patmo_idx_32SO) &
        - krate(:,382)*n(:,patmo_idx_36SO3)*n(:,patmo_idx_32SO)

    dn(:,patmo_idx_N2) = 0d0
        ! - krate(:,122)*n(:,patmo_idx_N2) &
        ! + krate(:,139)*n(:,patmo_idx_32SO2_1) &
        ! + krate(:,140)*n(:,patmo_idx_33SO2_1) &
        ! + krate(:,141)*n(:,patmo_idx_34SO2_1) &
        ! + krate(:,142)*n(:,patmo_idx_36SO2_1) &
        ! + krate(:,151)*n(:,patmo_idx_32SO2_1) &
        ! + krate(:,152)*n(:,patmo_idx_33SO2_1) &
        ! + krate(:,153)*n(:,patmo_idx_34SO2_1) &
        ! + krate(:,154)*n(:,patmo_idx_36SO2_1) &
        ! + krate(:,163)*n(:,patmo_idx_32SO2_3) &
        ! + krate(:,164)*n(:,patmo_idx_33SO2_3) &
        ! + krate(:,165)*n(:,patmo_idx_34SO2_3) &
        ! + krate(:,166)*n(:,patmo_idx_36SO2_3) &
        ! + krate(:,330)*n(:,patmo_idx_N)*n(:,patmo_idx_N) &
        ! - krate(:,347)*n(:,patmo_idx_32SO2)*n(:,patmo_idx_N2) &
        ! - krate(:,348)*n(:,patmo_idx_33SO2)*n(:,patmo_idx_N2) &
        ! - krate(:,349)*n(:,patmo_idx_34SO2)*n(:,patmo_idx_N2) &
        ! - krate(:,350)*n(:,patmo_idx_36SO2)*n(:,patmo_idx_N2) &
        ! - krate(:,359)*n(:,patmo_idx_32SO2_3)*n(:,patmo_idx_N2) &
        ! - krate(:,360)*n(:,patmo_idx_33SO2_3)*n(:,patmo_idx_N2) &
        ! - krate(:,361)*n(:,patmo_idx_34SO2_3)*n(:,patmo_idx_N2) &
        ! - krate(:,362)*n(:,patmo_idx_36SO2_3)*n(:,patmo_idx_N2) &
        ! - krate(:,371)*n(:,patmo_idx_32SO2)*n(:,patmo_idx_N2) &
        ! - krate(:,372)*n(:,patmo_idx_33SO2)*n(:,patmo_idx_N2) &
        ! - krate(:,373)*n(:,patmo_idx_34SO2)*n(:,patmo_idx_N2) &
        ! - krate(:,374)*n(:,patmo_idx_36SO2)*n(:,patmo_idx_N2)

    dn(:,patmo_idx_34SH) = &
        + krate(:,3)*n(:,patmo_idx_34COS)*n(:,patmo_idx_OH) &
        + krate(:,11)*n(:,patmo_idx_34CS2)*n(:,patmo_idx_OH) &
        + krate(:,31)*n(:,patmo_idx_34H2S)*n(:,patmo_idx_OH) &
        + krate(:,35)*n(:,patmo_idx_34H2S)*n(:,patmo_idx_O) &
        + krate(:,39)*n(:,patmo_idx_34H2S)*n(:,patmo_idx_H) &
        - krate(:,47)*n(:,patmo_idx_34SH)*n(:,patmo_idx_O) &
        - krate(:,51)*n(:,patmo_idx_34SH)*n(:,patmo_idx_O2) &
        - krate(:,55)*n(:,patmo_idx_34SH)*n(:,patmo_idx_O3) &
        + krate(:,95)*n(:,patmo_idx_34HSO)*n(:,patmo_idx_O3) &
        + krate(:,201)*n(:,patmo_idx_34H2S) &
        - krate(:,211)*n(:,patmo_idx_CO2)*n(:,patmo_idx_34SH) &
        - krate(:,219)*n(:,patmo_idx_34SH)*n(:,patmo_idx_34COS) &
        - krate(:,239)*n(:,patmo_idx_H2O)*n(:,patmo_idx_34SH) &
        - krate(:,243)*n(:,patmo_idx_OH)*n(:,patmo_idx_34SH) &
        - krate(:,247)*n(:,patmo_idx_H2)*n(:,patmo_idx_34SH) &
        + krate(:,255)*n(:,patmo_idx_H)*n(:,patmo_idx_34SO) &
        + krate(:,259)*n(:,patmo_idx_OH)*n(:,patmo_idx_34SO) &
        + krate(:,263)*n(:,patmo_idx_34HSO)*n(:,patmo_idx_O2) &
        - krate(:,303)*n(:,patmo_idx_O2)*n(:,patmo_idx_O2)*n(:,patmo_idx_34SH)

    dn(:,patmo_idx_32HSO2) = &
        - krate(:,97)*n(:,patmo_idx_32HSO2)*n(:,patmo_idx_O2) &
        + krate(:,305)*n(:,patmo_idx_HO2)*n(:,patmo_idx_32SO2)

    dn(:,patmo_idx_36CH3SCH3) = &
        - krate(:,130)*n(:,patmo_idx_36CH3SCH3)*n(:,patmo_idx_O) &
        - krate(:,134)*n(:,patmo_idx_36CH3SCH3)*n(:,patmo_idx_OH) &
        - krate(:,138)*n(:,patmo_idx_36CH3SCH3)*n(:,patmo_idx_OH) &
        + krate(:,338)*n(:,patmo_idx_36SO2) &
        + krate(:,342)*n(:,patmo_idx_36SO2) &
        + krate(:,346)*n(:,patmo_idx_36SO2)*n(:,patmo_idx_CH4O3S) 

    dn(:,patmo_idx_33SO3) = &
        + krate(:,82)*n(:,patmo_idx_33SO2)*n(:,patmo_idx_HO2) &
        + krate(:,86)*n(:,patmo_idx_33SO2)*n(:,patmo_idx_O3) &
        + krate(:,102)*n(:,patmo_idx_33HSO3)*n(:,patmo_idx_O2) &
        + krate(:,106)*n(:,patmo_idx_33SO2)*n(:,patmo_idx_O) &
        - krate(:,114)*n(:,patmo_idx_33SO3)*n(:,patmo_idx_H2O) &
        + krate(:,148)*n(:,patmo_idx_33SO2_1)*n(:,patmo_idx_32SO2) &
        + krate(:,168)*n(:,patmo_idx_33SO2_3)*n(:,patmo_idx_O2) &
        + krate(:,172)*n(:,patmo_idx_33SO2_3)*n(:,patmo_idx_32SO2) &
        - krate(:,188)*n(:,patmo_idx_33SO3) &
        - krate(:,290)*n(:,patmo_idx_OH)*n(:,patmo_idx_33SO3) &
        - krate(:,294)*n(:,patmo_idx_33SO3)*n(:,patmo_idx_O2) &
        - krate(:,310)*n(:,patmo_idx_HO2)*n(:,patmo_idx_33SO3) &
        - krate(:,314)*n(:,patmo_idx_33SO3) &
        + krate(:,322)*n(:,patmo_idx_33H2SO4) &
        - krate(:,356)*n(:,patmo_idx_33SO3)*n(:,patmo_idx_32SO) &
        - krate(:,376)*n(:,patmo_idx_33SO3)*n(:,patmo_idx_O) &
        - krate(:,380)*n(:,patmo_idx_33SO3)*n(:,patmo_idx_32SO)

    dn(:,patmo_idx_36HSO2) = &
        - krate(:,100)*n(:,patmo_idx_36HSO2)*n(:,patmo_idx_O2) &
        + krate(:,308)*n(:,patmo_idx_HO2)*n(:,patmo_idx_36SO2)

    dn(:,patmo_idx_34HSO3) = &
        - krate(:,103)*n(:,patmo_idx_34HSO3)*n(:,patmo_idx_O2) &
        + krate(:,111)*n(:,patmo_idx_34SO2)*n(:,patmo_idx_OH) &
        + krate(:,311)*n(:,patmo_idx_HO2)*n(:,patmo_idx_34SO3) &
        - krate(:,319)*n(:,patmo_idx_34HSO3)

    dn(:,patmo_idx_36SO2_3) = &
        + krate(:,154)*n(:,patmo_idx_36SO2_1) &
        + krate(:,158)*n(:,patmo_idx_36SO2_1) &
        - krate(:,162)*n(:,patmo_idx_36SO2_3) &
        - krate(:,166)*n(:,patmo_idx_36SO2_3) &
        - krate(:,170)*n(:,patmo_idx_36SO2_3)*n(:,patmo_idx_O2) &
        - krate(:,174)*n(:,patmo_idx_36SO2_3)*n(:,patmo_idx_32SO2) &
        + krate(:,177)*n(:,patmo_idx_36SO2) &
        - krate(:,362)*n(:,patmo_idx_36SO2_3)*n(:,patmo_idx_N2) &
        - krate(:,366)*n(:,patmo_idx_36SO2_3) &
        + krate(:,370)*n(:,patmo_idx_36SO2) &
        + krate(:,374)*n(:,patmo_idx_36SO2)*n(:,patmo_idx_N2) &
        + krate(:,378)*n(:,patmo_idx_36SO3)*n(:,patmo_idx_O) &
        + krate(:,382)*n(:,patmo_idx_36SO3)*n(:,patmo_idx_32SO)

     if(end_of_run) then
             budget(:,patmo_idx_32CS2,1) = &
                + krate(:,217)*n(:,patmo_idx_32SH)*n(:,patmo_idx_32COS) &
                + krate(:,221)*n(:,patmo_idx_32CS)*n(:,patmo_idx_32SO)
             
             budget(:,patmo_idx_32CS2,2) =  & 
                - krate(:,9)*n(:,patmo_idx_32CS2)*n(:,patmo_idx_OH) &
                - krate(:,13)*n(:,patmo_idx_32CS2)*n(:,patmo_idx_O) &
                - krate(:,202)*n(:,patmo_idx_32CS2) 
     
     endif

    dn(:,patmo_idx_32CS2) = &
        - krate(:,9)*n(:,patmo_idx_32CS2)*n(:,patmo_idx_OH) &
        - krate(:,13)*n(:,patmo_idx_32CS2)*n(:,patmo_idx_O) &
        - krate(:,202)*n(:,patmo_idx_32CS2) &
        + krate(:,217)*n(:,patmo_idx_32SH)*n(:,patmo_idx_32COS) &
        + krate(:,221)*n(:,patmo_idx_32CS)*n(:,patmo_idx_32SO)

    dn(:,patmo_idx_OH) = 0d0
        ! - krate(:,1)*n(:,patmo_idx_32COS)*n(:,patmo_idx_OH) &
        ! - krate(:,2)*n(:,patmo_idx_33COS)*n(:,patmo_idx_OH) &
        ! - krate(:,3)*n(:,patmo_idx_34COS)*n(:,patmo_idx_OH) &
        ! - krate(:,4)*n(:,patmo_idx_36COS)*n(:,patmo_idx_OH) &
        ! - krate(:,9)*n(:,patmo_idx_32CS2)*n(:,patmo_idx_OH) &
        ! - krate(:,10)*n(:,patmo_idx_33CS2)*n(:,patmo_idx_OH) &
        ! - krate(:,11)*n(:,patmo_idx_34CS2)*n(:,patmo_idx_OH) &
        ! - krate(:,12)*n(:,patmo_idx_36CS2)*n(:,patmo_idx_OH) &
        ! - krate(:,29)*n(:,patmo_idx_32H2S)*n(:,patmo_idx_OH) &
        ! - krate(:,30)*n(:,patmo_idx_33H2S)*n(:,patmo_idx_OH) &
        ! - krate(:,31)*n(:,patmo_idx_34H2S)*n(:,patmo_idx_OH) &
        ! - krate(:,32)*n(:,patmo_idx_36H2S)*n(:,patmo_idx_OH) &
        ! + krate(:,33)*n(:,patmo_idx_32H2S)*n(:,patmo_idx_O) &
        ! + krate(:,34)*n(:,patmo_idx_33H2S)*n(:,patmo_idx_O) &
        ! + krate(:,35)*n(:,patmo_idx_34H2S)*n(:,patmo_idx_O) &
        ! + krate(:,36)*n(:,patmo_idx_36H2S)*n(:,patmo_idx_O) &
        ! + krate(:,49)*n(:,patmo_idx_32SH)*n(:,patmo_idx_O2) &
        ! + krate(:,50)*n(:,patmo_idx_33SH)*n(:,patmo_idx_O2) &
        ! + krate(:,51)*n(:,patmo_idx_34SH)*n(:,patmo_idx_O2) &
        ! + krate(:,52)*n(:,patmo_idx_36SH)*n(:,patmo_idx_O2) &
        ! - krate(:,65)*n(:,patmo_idx_32SO)*n(:,patmo_idx_OH) &
        ! - krate(:,66)*n(:,patmo_idx_33SO)*n(:,patmo_idx_OH) &
        ! - krate(:,67)*n(:,patmo_idx_34SO)*n(:,patmo_idx_OH) &
        ! - krate(:,68)*n(:,patmo_idx_36SO)*n(:,patmo_idx_OH) &
        ! - krate(:,77)*n(:,patmo_idx_32S)*n(:,patmo_idx_OH) &
        ! - krate(:,78)*n(:,patmo_idx_33S)*n(:,patmo_idx_OH) &
        ! - krate(:,79)*n(:,patmo_idx_34S)*n(:,patmo_idx_OH) &
        ! - krate(:,80)*n(:,patmo_idx_36S)*n(:,patmo_idx_OH) &
        ! + krate(:,81)*n(:,patmo_idx_32SO2)*n(:,patmo_idx_HO2) &
        ! + krate(:,82)*n(:,patmo_idx_33SO2)*n(:,patmo_idx_HO2) &
        ! + krate(:,83)*n(:,patmo_idx_34SO2)*n(:,patmo_idx_HO2) &
        ! + krate(:,84)*n(:,patmo_idx_36SO2)*n(:,patmo_idx_HO2) &
        ! + krate(:,89)*n(:,patmo_idx_32HSO)*n(:,patmo_idx_O2) &
        ! + krate(:,90)*n(:,patmo_idx_33HSO)*n(:,patmo_idx_O2) &
        ! + krate(:,91)*n(:,patmo_idx_34HSO)*n(:,patmo_idx_O2) &
        ! + krate(:,92)*n(:,patmo_idx_36HSO)*n(:,patmo_idx_O2) &
        ! - krate(:,109)*n(:,patmo_idx_32SO2)*n(:,patmo_idx_OH) &
        ! - krate(:,110)*n(:,patmo_idx_33SO2)*n(:,patmo_idx_OH) &
        ! - krate(:,111)*n(:,patmo_idx_34SO2)*n(:,patmo_idx_OH) &
        ! - krate(:,112)*n(:,patmo_idx_36SO2)*n(:,patmo_idx_OH) &
        ! + krate(:,117)*n(:,patmo_idx_32H2SO4) &
        ! + krate(:,117)*n(:,patmo_idx_32H2SO4) &
        ! + krate(:,118)*n(:,patmo_idx_33H2SO4) &
        ! + krate(:,118)*n(:,patmo_idx_33H2SO4) &
        ! + krate(:,119)*n(:,patmo_idx_34H2SO4) &
        ! + krate(:,119)*n(:,patmo_idx_34H2SO4) &
        ! + krate(:,120)*n(:,patmo_idx_36H2SO4) &
        ! + krate(:,120)*n(:,patmo_idx_36H2SO4) &
        ! - krate(:,131)*n(:,patmo_idx_32CH3SCH3)*n(:,patmo_idx_OH) &
        ! - krate(:,132)*n(:,patmo_idx_33CH3SCH3)*n(:,patmo_idx_OH) &
        ! - krate(:,133)*n(:,patmo_idx_34CH3SCH3)*n(:,patmo_idx_OH) &
        ! - krate(:,134)*n(:,patmo_idx_36CH3SCH3)*n(:,patmo_idx_OH) &
        ! - krate(:,135)*n(:,patmo_idx_32CH3SCH3)*n(:,patmo_idx_OH) &
        ! - krate(:,136)*n(:,patmo_idx_33CH3SCH3)*n(:,patmo_idx_OH) &
        ! - krate(:,137)*n(:,patmo_idx_34CH3SCH3)*n(:,patmo_idx_OH) &
        ! - krate(:,138)*n(:,patmo_idx_36CH3SCH3)*n(:,patmo_idx_OH) &
        ! + krate(:,209)*n(:,patmo_idx_CO2)*n(:,patmo_idx_32SH) &
        ! + krate(:,210)*n(:,patmo_idx_CO2)*n(:,patmo_idx_33SH) &
        ! + krate(:,211)*n(:,patmo_idx_CO2)*n(:,patmo_idx_34SH) &
        ! + krate(:,212)*n(:,patmo_idx_CO2)*n(:,patmo_idx_36SH) &
        ! + krate(:,217)*n(:,patmo_idx_32SH)*n(:,patmo_idx_32COS) &
        ! + krate(:,218)*n(:,patmo_idx_33SH)*n(:,patmo_idx_33COS) &
        ! + krate(:,219)*n(:,patmo_idx_34SH)*n(:,patmo_idx_34COS) &
        ! + krate(:,220)*n(:,patmo_idx_36SH)*n(:,patmo_idx_36COS) &
        ! + krate(:,237)*n(:,patmo_idx_H2O)*n(:,patmo_idx_32SH) &
        ! + krate(:,238)*n(:,patmo_idx_H2O)*n(:,patmo_idx_33SH) &
        ! + krate(:,239)*n(:,patmo_idx_H2O)*n(:,patmo_idx_34SH) &
        ! + krate(:,240)*n(:,patmo_idx_H2O)*n(:,patmo_idx_36SH) &
        ! - krate(:,241)*n(:,patmo_idx_OH)*n(:,patmo_idx_32SH) &
        ! - krate(:,242)*n(:,patmo_idx_OH)*n(:,patmo_idx_33SH) &
        ! - krate(:,243)*n(:,patmo_idx_OH)*n(:,patmo_idx_34SH) &
        ! - krate(:,244)*n(:,patmo_idx_OH)*n(:,patmo_idx_36SH) &
        ! - krate(:,257)*n(:,patmo_idx_OH)*n(:,patmo_idx_32SO) &
        ! - krate(:,258)*n(:,patmo_idx_OH)*n(:,patmo_idx_33SO) &
        ! - krate(:,259)*n(:,patmo_idx_OH)*n(:,patmo_idx_34SO) &
        ! - krate(:,260)*n(:,patmo_idx_OH)*n(:,patmo_idx_36SO) &
        ! + krate(:,273)*n(:,patmo_idx_32SO2)*n(:,patmo_idx_H) &
        ! + krate(:,274)*n(:,patmo_idx_33SO2)*n(:,patmo_idx_H) &
        ! + krate(:,275)*n(:,patmo_idx_34SO2)*n(:,patmo_idx_H) &
        ! + krate(:,276)*n(:,patmo_idx_36SO2)*n(:,patmo_idx_H) &
        ! + krate(:,285)*n(:,patmo_idx_H)*n(:,patmo_idx_32SO) &
        ! + krate(:,286)*n(:,patmo_idx_H)*n(:,patmo_idx_33SO) &
        ! + krate(:,287)*n(:,patmo_idx_H)*n(:,patmo_idx_34SO) &
        ! + krate(:,288)*n(:,patmo_idx_H)*n(:,patmo_idx_36SO) &
        ! - krate(:,289)*n(:,patmo_idx_OH)*n(:,patmo_idx_32SO3) &
        ! - krate(:,290)*n(:,patmo_idx_OH)*n(:,patmo_idx_33SO3) &
        ! - krate(:,291)*n(:,patmo_idx_OH)*n(:,patmo_idx_34SO3) &
        ! - krate(:,292)*n(:,patmo_idx_OH)*n(:,patmo_idx_36SO3) &
        ! - krate(:,297)*n(:,patmo_idx_32SO2)*n(:,patmo_idx_OH) &
        ! - krate(:,298)*n(:,patmo_idx_33SO2)*n(:,patmo_idx_OH) &
        ! - krate(:,299)*n(:,patmo_idx_34SO2)*n(:,patmo_idx_OH) &
        ! - krate(:,300)*n(:,patmo_idx_36SO2)*n(:,patmo_idx_OH) &
        ! + krate(:,317)*n(:,patmo_idx_32HSO3) &
        ! + krate(:,318)*n(:,patmo_idx_33HSO3) &
        ! + krate(:,319)*n(:,patmo_idx_34HSO3) &
        ! + krate(:,320)*n(:,patmo_idx_36HSO3) &
        ! - krate(:,325)*n(:,patmo_idx_32SO2)*n(:,patmo_idx_OH)*n(:,patmo_idx_OH) &
        ! - krate(:,325)*n(:,patmo_idx_32SO2)*n(:,patmo_idx_OH)*n(:,patmo_idx_OH) &
        ! - krate(:,326)*n(:,patmo_idx_33SO2)*n(:,patmo_idx_OH)*n(:,patmo_idx_OH) &
        ! - krate(:,326)*n(:,patmo_idx_33SO2)*n(:,patmo_idx_OH)*n(:,patmo_idx_OH) &
        ! - krate(:,327)*n(:,patmo_idx_34SO2)*n(:,patmo_idx_OH)*n(:,patmo_idx_OH) &
        ! - krate(:,327)*n(:,patmo_idx_34SO2)*n(:,patmo_idx_OH)*n(:,patmo_idx_OH) &
        ! - krate(:,328)*n(:,patmo_idx_36SO2)*n(:,patmo_idx_OH)*n(:,patmo_idx_OH) &
        ! - krate(:,328)*n(:,patmo_idx_36SO2)*n(:,patmo_idx_OH)*n(:,patmo_idx_OH) &
        ! + krate(:,339)*n(:,patmo_idx_32SO2) &
        ! + krate(:,340)*n(:,patmo_idx_33SO2) &
        ! + krate(:,341)*n(:,patmo_idx_34SO2) &
        ! + krate(:,342)*n(:,patmo_idx_36SO2) &
        ! + krate(:,343)*n(:,patmo_idx_32SO2)*n(:,patmo_idx_CH4O3S) &
        ! + krate(:,344)*n(:,patmo_idx_33SO2)*n(:,patmo_idx_CH4O3S) &
        ! + krate(:,345)*n(:,patmo_idx_34SO2)*n(:,patmo_idx_CH4O3S) &
        ! + krate(:,346)*n(:,patmo_idx_36SO2)*n(:,patmo_idx_CH4O3S)

    dn(:,patmo_idx_O) = 0d0
        ! - krate(:,5)*n(:,patmo_idx_32COS)*n(:,patmo_idx_O) &
        ! - krate(:,6)*n(:,patmo_idx_33COS)*n(:,patmo_idx_O) &
        ! - krate(:,7)*n(:,patmo_idx_34COS)*n(:,patmo_idx_O) &
        ! - krate(:,8)*n(:,patmo_idx_36COS)*n(:,patmo_idx_O) &
        ! - krate(:,13)*n(:,patmo_idx_32CS2)*n(:,patmo_idx_O) &
        ! - krate(:,14)*n(:,patmo_idx_33CS2)*n(:,patmo_idx_O) &
        ! - krate(:,15)*n(:,patmo_idx_34CS2)*n(:,patmo_idx_O) &
        ! - krate(:,16)*n(:,patmo_idx_36CS2)*n(:,patmo_idx_O) &
        ! + krate(:,17)*n(:,patmo_idx_32CS)*n(:,patmo_idx_O2) &
        ! + krate(:,18)*n(:,patmo_idx_33CS)*n(:,patmo_idx_O2) &
        ! + krate(:,19)*n(:,patmo_idx_34CS)*n(:,patmo_idx_O2) &
        ! + krate(:,20)*n(:,patmo_idx_36CS)*n(:,patmo_idx_O2) &
        ! - krate(:,25)*n(:,patmo_idx_32CS)*n(:,patmo_idx_O) &
        ! - krate(:,26)*n(:,patmo_idx_33CS)*n(:,patmo_idx_O) &
        ! - krate(:,27)*n(:,patmo_idx_34CS)*n(:,patmo_idx_O) &
        ! - krate(:,28)*n(:,patmo_idx_36CS)*n(:,patmo_idx_O) &
        ! - krate(:,33)*n(:,patmo_idx_32H2S)*n(:,patmo_idx_O) &
        ! - krate(:,34)*n(:,patmo_idx_33H2S)*n(:,patmo_idx_O) &
        ! - krate(:,35)*n(:,patmo_idx_34H2S)*n(:,patmo_idx_O) &
        ! - krate(:,36)*n(:,patmo_idx_36H2S)*n(:,patmo_idx_O) &
        ! - krate(:,45)*n(:,patmo_idx_32SH)*n(:,patmo_idx_O) &
        ! - krate(:,46)*n(:,patmo_idx_33SH)*n(:,patmo_idx_O) &
        ! - krate(:,47)*n(:,patmo_idx_34SH)*n(:,patmo_idx_O) &
        ! - krate(:,48)*n(:,patmo_idx_36SH)*n(:,patmo_idx_O) &
        ! + krate(:,61)*n(:,patmo_idx_32SO)*n(:,patmo_idx_O2) &
        ! + krate(:,62)*n(:,patmo_idx_33SO)*n(:,patmo_idx_O2) &
        ! + krate(:,63)*n(:,patmo_idx_34SO)*n(:,patmo_idx_O2) &
        ! + krate(:,64)*n(:,patmo_idx_36SO)*n(:,patmo_idx_O2) &
        ! + krate(:,69)*n(:,patmo_idx_32S)*n(:,patmo_idx_O2) &
        ! + krate(:,70)*n(:,patmo_idx_33S)*n(:,patmo_idx_O2) &
        ! + krate(:,71)*n(:,patmo_idx_34S)*n(:,patmo_idx_O2) &
        ! + krate(:,72)*n(:,patmo_idx_36S)*n(:,patmo_idx_O2) &
        ! - krate(:,105)*n(:,patmo_idx_32SO2)*n(:,patmo_idx_O) &
        ! - krate(:,106)*n(:,patmo_idx_33SO2)*n(:,patmo_idx_O) &
        ! - krate(:,107)*n(:,patmo_idx_34SO2)*n(:,patmo_idx_O) &
        ! - krate(:,108)*n(:,patmo_idx_36SO2)*n(:,patmo_idx_O) &
        ! - krate(:,121)*n(:,patmo_idx_O)*n(:,patmo_idx_O2) &
        ! - krate(:,127)*n(:,patmo_idx_32CH3SCH3)*n(:,patmo_idx_O) &
        ! - krate(:,128)*n(:,patmo_idx_33CH3SCH3)*n(:,patmo_idx_O) &
        ! - krate(:,129)*n(:,patmo_idx_34CH3SCH3)*n(:,patmo_idx_O) &
        ! - krate(:,130)*n(:,patmo_idx_36CH3SCH3)*n(:,patmo_idx_O) &
        ! + krate(:,167)*n(:,patmo_idx_32SO2_3)*n(:,patmo_idx_O2) &
        ! + krate(:,168)*n(:,patmo_idx_33SO2_3)*n(:,patmo_idx_O2) &
        ! + krate(:,169)*n(:,patmo_idx_34SO2_3)*n(:,patmo_idx_O2) &
        ! + krate(:,170)*n(:,patmo_idx_36SO2_3)*n(:,patmo_idx_O2) &
        ! + krate(:,175)*n(:,patmo_idx_36SO2) &
        ! + krate(:,176)*n(:,patmo_idx_36SO2) &
        ! + krate(:,177)*n(:,patmo_idx_36SO2) &
        ! + krate(:,178)*n(:,patmo_idx_36SO3) &
        ! + krate(:,180)*n(:,patmo_idx_34SO2) &
        ! + krate(:,181)*n(:,patmo_idx_34SO2) &
        ! + krate(:,182)*n(:,patmo_idx_34SO2) &
        ! + krate(:,184)*n(:,patmo_idx_O3) &
        ! + krate(:,185)*n(:,patmo_idx_O2) &
        ! + krate(:,185)*n(:,patmo_idx_O2) &
        ! + krate(:,188)*n(:,patmo_idx_33SO3) &
        ! + krate(:,190)*n(:,patmo_idx_32SO2) &
        ! + krate(:,191)*n(:,patmo_idx_32SO2) &
        ! + krate(:,192)*n(:,patmo_idx_32SO2) &
        ! + krate(:,194)*n(:,patmo_idx_32SO3) &
        ! + krate(:,195)*n(:,patmo_idx_34SO) &
        ! + krate(:,196)*n(:,patmo_idx_32SO) &
        ! + krate(:,198)*n(:,patmo_idx_33SO2) &
        ! + krate(:,199)*n(:,patmo_idx_33SO2) &
        ! + krate(:,200)*n(:,patmo_idx_33SO2) &
        ! + krate(:,203)*n(:,patmo_idx_33SO) &
        ! + krate(:,204)*n(:,patmo_idx_36SO) &
        ! + krate(:,208)*n(:,patmo_idx_34SO3) &
        ! + krate(:,213)*n(:,patmo_idx_CO)*n(:,patmo_idx_32SO) &
        ! + krate(:,214)*n(:,patmo_idx_CO)*n(:,patmo_idx_33SO) &
        ! + krate(:,215)*n(:,patmo_idx_CO)*n(:,patmo_idx_34SO) &
        ! + krate(:,216)*n(:,patmo_idx_CO)*n(:,patmo_idx_36SO) &
        ! + krate(:,221)*n(:,patmo_idx_32CS)*n(:,patmo_idx_32SO) &
        ! + krate(:,222)*n(:,patmo_idx_33CS)*n(:,patmo_idx_33SO) &
        ! + krate(:,223)*n(:,patmo_idx_34CS)*n(:,patmo_idx_34SO) &
        ! + krate(:,224)*n(:,patmo_idx_36CS)*n(:,patmo_idx_36SO) &
        ! - krate(:,225)*n(:,patmo_idx_32COS)*n(:,patmo_idx_O) &
        ! - krate(:,226)*n(:,patmo_idx_33COS)*n(:,patmo_idx_O) &
        ! - krate(:,227)*n(:,patmo_idx_34COS)*n(:,patmo_idx_O) &
        ! - krate(:,228)*n(:,patmo_idx_36COS)*n(:,patmo_idx_O) &
        ! + krate(:,233)*n(:,patmo_idx_CO)*n(:,patmo_idx_32S) &
        ! + krate(:,234)*n(:,patmo_idx_CO)*n(:,patmo_idx_33S) &
        ! + krate(:,235)*n(:,patmo_idx_CO)*n(:,patmo_idx_34S) &
        ! + krate(:,236)*n(:,patmo_idx_CO)*n(:,patmo_idx_36S) &
        ! + krate(:,241)*n(:,patmo_idx_OH)*n(:,patmo_idx_32SH) &
        ! + krate(:,242)*n(:,patmo_idx_OH)*n(:,patmo_idx_33SH) &
        ! + krate(:,243)*n(:,patmo_idx_OH)*n(:,patmo_idx_34SH) &
        ! + krate(:,244)*n(:,patmo_idx_OH)*n(:,patmo_idx_36SH) &
        ! + krate(:,253)*n(:,patmo_idx_H)*n(:,patmo_idx_32SO) &
        ! + krate(:,254)*n(:,patmo_idx_H)*n(:,patmo_idx_33SO) &
        ! + krate(:,255)*n(:,patmo_idx_H)*n(:,patmo_idx_34SO) &
        ! + krate(:,256)*n(:,patmo_idx_H)*n(:,patmo_idx_36SO) &
        ! - krate(:,269)*n(:,patmo_idx_32SO2)*n(:,patmo_idx_O) &
        ! - krate(:,270)*n(:,patmo_idx_33SO2)*n(:,patmo_idx_O) &
        ! - krate(:,271)*n(:,patmo_idx_34SO2)*n(:,patmo_idx_O) &
        ! - krate(:,272)*n(:,patmo_idx_36SO2)*n(:,patmo_idx_O) &
        ! - krate(:,277)*n(:,patmo_idx_32SO)*n(:,patmo_idx_O) &
        ! - krate(:,278)*n(:,patmo_idx_33SO)*n(:,patmo_idx_O) &
        ! - krate(:,279)*n(:,patmo_idx_34SO)*n(:,patmo_idx_O) &
        ! - krate(:,280)*n(:,patmo_idx_36SO)*n(:,patmo_idx_O) &
        ! + krate(:,313)*n(:,patmo_idx_32SO3) &
        ! + krate(:,314)*n(:,patmo_idx_33SO3) &
        ! + krate(:,315)*n(:,patmo_idx_34SO3) &
        ! + krate(:,316)*n(:,patmo_idx_36SO3) &
        ! + krate(:,329)*n(:,patmo_idx_O3) &
        ! + krate(:,335)*n(:,patmo_idx_32SO2) &
        ! + krate(:,336)*n(:,patmo_idx_33SO2) &
        ! + krate(:,337)*n(:,patmo_idx_34SO2) &
        ! + krate(:,338)*n(:,patmo_idx_36SO2) &
        ! - krate(:,375)*n(:,patmo_idx_32SO3)*n(:,patmo_idx_O) &
        ! - krate(:,376)*n(:,patmo_idx_33SO3)*n(:,patmo_idx_O) &
        ! - krate(:,377)*n(:,patmo_idx_34SO3)*n(:,patmo_idx_O) &
        ! - krate(:,378)*n(:,patmo_idx_36SO3)*n(:,patmo_idx_O)

    if(end_of_run)then
            budget(:,patmo_idx_34H2S,1) = &
        + krate(:,239)*n(:,patmo_idx_H2O)*n(:,patmo_idx_34SH) &
        + krate(:,243)*n(:,patmo_idx_OH)*n(:,patmo_idx_34SH) &
        + krate(:,247)*n(:,patmo_idx_H2)*n(:,patmo_idx_34SH) &
        + krate(:,251)*n(:,patmo_idx_H2O)*n(:,patmo_idx_34HSO)
            budget(:,patmo_idx_34H2S,2) = &
                            - krate(:,31)*n(:,patmo_idx_34H2S)*n(:,patmo_idx_OH) &
        - krate(:,35)*n(:,patmo_idx_34H2S)*n(:,patmo_idx_O) &
        - krate(:,39)*n(:,patmo_idx_34H2S)*n(:,patmo_idx_H) &
        - krate(:,43)*n(:,patmo_idx_34H2S)*n(:,patmo_idx_HO2) 
       endif


    dn(:,patmo_idx_34H2S) = &
        - krate(:,31)*n(:,patmo_idx_34H2S)*n(:,patmo_idx_OH) &
        - krate(:,35)*n(:,patmo_idx_34H2S)*n(:,patmo_idx_O) &
        - krate(:,39)*n(:,patmo_idx_34H2S)*n(:,patmo_idx_H) &
        - krate(:,43)*n(:,patmo_idx_34H2S)*n(:,patmo_idx_HO2) &
        - krate(:,201)*n(:,patmo_idx_34H2S) &
        + krate(:,239)*n(:,patmo_idx_H2O)*n(:,patmo_idx_34SH) &
        + krate(:,243)*n(:,patmo_idx_OH)*n(:,patmo_idx_34SH) &
        + krate(:,247)*n(:,patmo_idx_H2)*n(:,patmo_idx_34SH) &
        + krate(:,251)*n(:,patmo_idx_H2O)*n(:,patmo_idx_34HSO)
!H2SO4 chemistry budget

    if (end_of_run) then 
        budget(:,patmo_idx_32H2SO4,1) = &
        + krate(:,113)*n(:,patmo_idx_32SO3)*n(:,patmo_idx_H2O) &
        + krate(:,325)*n(:,patmo_idx_32SO2)*n(:,patmo_idx_OH)*n(:,patmo_idx_OH)
        
        budget(:,patmo_idx_32H2SO4,8) = &
        + krate(:,325)*n(:,patmo_idx_32SO2)*n(:,patmo_idx_OH)*n(:,patmo_idx_OH)

        budget(:,patmo_idx_32H2SO4,9) = &
        + krate(:,113)*n(:,patmo_idx_32SO3)*n(:,patmo_idx_H2O) 

        budget(:,patmo_idx_32H2SO4,2) = &
        + krate(:,117)*n(:,patmo_idx_32H2SO4) &
        + krate(:,321)*n(:,patmo_idx_32H2SO4) 

        budget(:,patmo_idx_32H2SO4,10) = &
        + krate(:,117)*n(:,patmo_idx_32H2SO4) 

    endif

    dn(:,patmo_idx_32H2SO4) = &
        + krate(:,113)*n(:,patmo_idx_32SO3)*n(:,patmo_idx_H2O) &
        - krate(:,117)*n(:,patmo_idx_32H2SO4) &
        - krate(:,321)*n(:,patmo_idx_32H2SO4) &
        + krate(:,325)*n(:,patmo_idx_32SO2)*n(:,patmo_idx_OH)*n(:,patmo_idx_OH)

    dn(:,patmo_idx_34SO2_1) = &
        - krate(:,141)*n(:,patmo_idx_34SO2_1) &
        - krate(:,145)*n(:,patmo_idx_34SO2_1) &
        - krate(:,149)*n(:,patmo_idx_34SO2_1)*n(:,patmo_idx_32SO2) &
        - krate(:,153)*n(:,patmo_idx_34SO2_1) &
        - krate(:,157)*n(:,patmo_idx_34SO2_1) &
        + krate(:,182)*n(:,patmo_idx_34SO2) &
        + krate(:,349)*n(:,patmo_idx_34SO2)*n(:,patmo_idx_N2) &
        + krate(:,353)*n(:,patmo_idx_34SO2) &
        + krate(:,357)*n(:,patmo_idx_34SO3)*n(:,patmo_idx_32SO) &
        + krate(:,361)*n(:,patmo_idx_34SO2_3)*n(:,patmo_idx_N2) &
        + krate(:,365)*n(:,patmo_idx_34SO2_3)

    dn(:,patmo_idx_33SH) = &
        + krate(:,2)*n(:,patmo_idx_33COS)*n(:,patmo_idx_OH) &
        + krate(:,10)*n(:,patmo_idx_33CS2)*n(:,patmo_idx_OH) &
        + krate(:,30)*n(:,patmo_idx_33H2S)*n(:,patmo_idx_OH) &
        + krate(:,34)*n(:,patmo_idx_33H2S)*n(:,patmo_idx_O) &
        + krate(:,38)*n(:,patmo_idx_33H2S)*n(:,patmo_idx_H) &
        - krate(:,46)*n(:,patmo_idx_33SH)*n(:,patmo_idx_O) &
        - krate(:,50)*n(:,patmo_idx_33SH)*n(:,patmo_idx_O2) &
        - krate(:,54)*n(:,patmo_idx_33SH)*n(:,patmo_idx_O3) &
        + krate(:,94)*n(:,patmo_idx_33HSO)*n(:,patmo_idx_O3) &
        + krate(:,189)*n(:,patmo_idx_33H2S) &
        - krate(:,210)*n(:,patmo_idx_CO2)*n(:,patmo_idx_33SH) &
        - krate(:,218)*n(:,patmo_idx_33SH)*n(:,patmo_idx_33COS) &
        - krate(:,238)*n(:,patmo_idx_H2O)*n(:,patmo_idx_33SH) &
        - krate(:,242)*n(:,patmo_idx_OH)*n(:,patmo_idx_33SH) &
        - krate(:,246)*n(:,patmo_idx_H2)*n(:,patmo_idx_33SH) &
        + krate(:,254)*n(:,patmo_idx_H)*n(:,patmo_idx_33SO) &
        + krate(:,258)*n(:,patmo_idx_OH)*n(:,patmo_idx_33SO) &
        + krate(:,262)*n(:,patmo_idx_33HSO)*n(:,patmo_idx_O2) &
        - krate(:,302)*n(:,patmo_idx_O2)*n(:,patmo_idx_O2)*n(:,patmo_idx_33SH)
     
    !if (end_of_run) then 
    !    budget(:,patmo_idx_33COS,1) = &
    !    + (0.83*(krate(:,10)*n(:,patmo_idx_33CS2)*n(:,patmo_idx_OH))) &
    !    + krate(:,18)*n(:,patmo_idx_33CS)*n(:,patmo_idx_O2) &
    !    + krate(:,22)*n(:,patmo_idx_33CS)*n(:,patmo_idx_O3) &
    !    + krate(:,210)*n(:,patmo_idx_CO2)*n(:,patmo_idx_33SH) &
    !    + krate(:,214)*n(:,patmo_idx_CO)*n(:,patmo_idx_33SO) &
    !    + 0.007*(krate(:,132))*n(:,patmo_idx_33CH3SCH3)*n(:,patmo_idx_OH)
    
    !    budget(:,patmo_idx_33COS,2) = &
    !    - krate(:,2)*n(:,patmo_idx_33COS)*n(:,patmo_idx_OH) &
    !    - krate(:,6)*n(:,patmo_idx_33COS)*n(:,patmo_idx_O) &
    !    - krate(:,197)*n(:,patmo_idx_33COS) &
    !    - krate(:,218)*n(:,patmo_idx_33SH)*n(:,patmo_idx_33COS) &
    !    - krate(:,226)*n(:,patmo_idx_33COS)*n(:,patmo_idx_O) &
    !    - krate(:,230)*n(:,patmo_idx_33COS)*n(:,patmo_idx_O2)
    !endif

    dn(:,patmo_idx_33COS) = &
        - krate(:,2)*n(:,patmo_idx_33COS)*n(:,patmo_idx_OH) &
        - krate(:,6)*n(:,patmo_idx_33COS)*n(:,patmo_idx_O) &
        + (0.83*(krate(:,10)*n(:,patmo_idx_33CS2)*n(:,patmo_idx_OH))) &
        + krate(:,18)*n(:,patmo_idx_33CS)*n(:,patmo_idx_O2) &
        + krate(:,22)*n(:,patmo_idx_33CS)*n(:,patmo_idx_O3) &
        - krate(:,197)*n(:,patmo_idx_33COS) &
        + krate(:,210)*n(:,patmo_idx_CO2)*n(:,patmo_idx_33SH) &
        + krate(:,214)*n(:,patmo_idx_CO)*n(:,patmo_idx_33SO) &
        - krate(:,218)*n(:,patmo_idx_33SH)*n(:,patmo_idx_33COS) &
        - krate(:,226)*n(:,patmo_idx_33COS)*n(:,patmo_idx_O) &
        - krate(:,230)*n(:,patmo_idx_33COS)*n(:,patmo_idx_O2) &
        + 0.007*(krate(:,132)*n(:,patmo_idx_33CH3SCH3)*n(:,patmo_idx_OH)) &
        + 0.007*(krate(:,136)*n(:,patmo_idx_33CH3SCH3)*n(:,patmo_idx_OH)) &
        + 0.007*(krate(:,128)*n(:,patmo_idx_33CH3SCH3)*n(:,patmo_idx_O))

    dn(:,patmo_idx_36CS) = &
        + krate(:,16)*n(:,patmo_idx_36CS2)*n(:,patmo_idx_O) &
        - krate(:,20)*n(:,patmo_idx_36CS)*n(:,patmo_idx_O2) &
        - krate(:,24)*n(:,patmo_idx_36CS)*n(:,patmo_idx_O3) &
        - krate(:,28)*n(:,patmo_idx_36CS)*n(:,patmo_idx_O) &
        + krate(:,205)*n(:,patmo_idx_36CS2) &
        - krate(:,224)*n(:,patmo_idx_36CS)*n(:,patmo_idx_36SO) &
        + krate(:,228)*n(:,patmo_idx_36COS)*n(:,patmo_idx_O) &
        + krate(:,232)*n(:,patmo_idx_36COS)*n(:,patmo_idx_O2) &
        + krate(:,236)*n(:,patmo_idx_CO)*n(:,patmo_idx_36S)

        if (end_of_run) then
                budget(:,patmo_idx_32CH3SCH3,1) = & 
                         + krate(:,335)*n(:,patmo_idx_32SO2) &
                         + krate(:,339)*n(:,patmo_idx_32SO2) &
                         + krate(:,343)*n(:,patmo_idx_32SO2)*n(:,patmo_idx_CH4O3S) 

                budget(:,patmo_idx_32CH3SCH3,2) = & 
                        - krate(:,127)*n(:,patmo_idx_32CH3SCH3)*n(:,patmo_idx_O) &
                        - krate(:,131)*n(:,patmo_idx_32CH3SCH3)*n(:,patmo_idx_OH) &
                        - krate(:,135)*n(:,patmo_idx_32CH3SCH3)*n(:,patmo_idx_OH) 
        endif


    dn(:,patmo_idx_32CH3SCH3) = &
        - krate(:,127)*n(:,patmo_idx_32CH3SCH3)*n(:,patmo_idx_O) &
        - krate(:,131)*n(:,patmo_idx_32CH3SCH3)*n(:,patmo_idx_OH) &
        - krate(:,135)*n(:,patmo_idx_32CH3SCH3)*n(:,patmo_idx_OH) &
        + krate(:,335)*n(:,patmo_idx_32SO2) &
        + krate(:,339)*n(:,patmo_idx_32SO2) &
        + krate(:,343)*n(:,patmo_idx_32SO2)*n(:,patmo_idx_CH4O3S) 
       ! - (7.940d-6)*n(:,patmo_idx_32CH3SCH3)*n(:,patmo_idx_OH)

    dn(:,patmo_idx_33H2S) = &
        - krate(:,30)*n(:,patmo_idx_33H2S)*n(:,patmo_idx_OH) &
        - krate(:,34)*n(:,patmo_idx_33H2S)*n(:,patmo_idx_O) &
        - krate(:,38)*n(:,patmo_idx_33H2S)*n(:,patmo_idx_H) &
        - krate(:,42)*n(:,patmo_idx_33H2S)*n(:,patmo_idx_HO2) &
        - krate(:,189)*n(:,patmo_idx_33H2S) &
        + krate(:,238)*n(:,patmo_idx_H2O)*n(:,patmo_idx_33SH) &
        + krate(:,242)*n(:,patmo_idx_OH)*n(:,patmo_idx_33SH) &
        + krate(:,246)*n(:,patmo_idx_H2)*n(:,patmo_idx_33SH) &
        + krate(:,250)*n(:,patmo_idx_H2O)*n(:,patmo_idx_33HSO)

    dn(:,patmo_idx_H) = 0d0
        ! - krate(:,37)*n(:,patmo_idx_32H2S)*n(:,patmo_idx_H) &
        ! - krate(:,38)*n(:,patmo_idx_33H2S)*n(:,patmo_idx_H) &
        ! - krate(:,39)*n(:,patmo_idx_34H2S)*n(:,patmo_idx_H) &
        ! - krate(:,40)*n(:,patmo_idx_36H2S)*n(:,patmo_idx_H) &
        ! + krate(:,45)*n(:,patmo_idx_32SH)*n(:,patmo_idx_O) &
        ! + krate(:,46)*n(:,patmo_idx_33SH)*n(:,patmo_idx_O) &
        ! + krate(:,47)*n(:,patmo_idx_34SH)*n(:,patmo_idx_O) &
        ! + krate(:,48)*n(:,patmo_idx_36SH)*n(:,patmo_idx_O) &
        ! + krate(:,65)*n(:,patmo_idx_32SO)*n(:,patmo_idx_OH) &
        ! + krate(:,66)*n(:,patmo_idx_33SO)*n(:,patmo_idx_OH) &
        ! + krate(:,67)*n(:,patmo_idx_34SO)*n(:,patmo_idx_OH) &
        ! + krate(:,68)*n(:,patmo_idx_36SO)*n(:,patmo_idx_OH) &
        ! + krate(:,77)*n(:,patmo_idx_32S)*n(:,patmo_idx_OH) &
        ! + krate(:,78)*n(:,patmo_idx_33S)*n(:,patmo_idx_OH) &
        ! + krate(:,79)*n(:,patmo_idx_34S)*n(:,patmo_idx_OH) &
        ! + krate(:,80)*n(:,patmo_idx_36S)*n(:,patmo_idx_OH) &
        ! + krate(:,187)*n(:,patmo_idx_36H2S) &
        ! + krate(:,189)*n(:,patmo_idx_33H2S) &
        ! + krate(:,201)*n(:,patmo_idx_34H2S) &
        ! + krate(:,206)*n(:,patmo_idx_32H2S) &
        ! + krate(:,245)*n(:,patmo_idx_H2)*n(:,patmo_idx_32SH) &
        ! + krate(:,246)*n(:,patmo_idx_H2)*n(:,patmo_idx_33SH) &
        ! + krate(:,247)*n(:,patmo_idx_H2)*n(:,patmo_idx_34SH) &
        ! + krate(:,248)*n(:,patmo_idx_H2)*n(:,patmo_idx_36SH) &
        ! - krate(:,253)*n(:,patmo_idx_H)*n(:,patmo_idx_32SO) &
        ! - krate(:,254)*n(:,patmo_idx_H)*n(:,patmo_idx_33SO) &
        ! - krate(:,255)*n(:,patmo_idx_H)*n(:,patmo_idx_34SO) &
        ! - krate(:,256)*n(:,patmo_idx_H)*n(:,patmo_idx_36SO) &
        ! - krate(:,273)*n(:,patmo_idx_32SO2)*n(:,patmo_idx_H) &
        ! - krate(:,274)*n(:,patmo_idx_33SO2)*n(:,patmo_idx_H) &
        ! - krate(:,275)*n(:,patmo_idx_34SO2)*n(:,patmo_idx_H) &
        ! - krate(:,276)*n(:,patmo_idx_36SO2)*n(:,patmo_idx_H) &
        ! - krate(:,285)*n(:,patmo_idx_H)*n(:,patmo_idx_32SO) &
        ! - krate(:,286)*n(:,patmo_idx_H)*n(:,patmo_idx_33SO) &
        ! - krate(:,287)*n(:,patmo_idx_H)*n(:,patmo_idx_34SO) &
        ! - krate(:,288)*n(:,patmo_idx_H)*n(:,patmo_idx_36SO)

    dn(:,patmo_idx_33SO2_1) = &
        - krate(:,140)*n(:,patmo_idx_33SO2_1) &
        - krate(:,144)*n(:,patmo_idx_33SO2_1) &
        - krate(:,148)*n(:,patmo_idx_33SO2_1)*n(:,patmo_idx_32SO2) &
        - krate(:,152)*n(:,patmo_idx_33SO2_1) &
        - krate(:,156)*n(:,patmo_idx_33SO2_1) &
        + krate(:,200)*n(:,patmo_idx_33SO2) &
        + krate(:,348)*n(:,patmo_idx_33SO2)*n(:,patmo_idx_N2) &
        + krate(:,352)*n(:,patmo_idx_33SO2) &
        + krate(:,356)*n(:,patmo_idx_33SO3)*n(:,patmo_idx_32SO) &
        + krate(:,360)*n(:,patmo_idx_33SO2_3)*n(:,patmo_idx_N2) &
        + krate(:,364)*n(:,patmo_idx_33SO2_3)

    if (end_of_run) then
        budget(:,patmo_idx_32SO2,1) = &
        + krate(:,57)*n(:,patmo_idx_32SO)*n(:,patmo_idx_O3) &
        + krate(:,61)*n(:,patmo_idx_32SO)*n(:,patmo_idx_O2) &
        + krate(:,65)*n(:,patmo_idx_32SO)*n(:,patmo_idx_OH) &
        + krate(:,89)*n(:,patmo_idx_32HSO)*n(:,patmo_idx_O2) &
        + krate(:,97)*n(:,patmo_idx_32HSO2)*n(:,patmo_idx_O2) &
        + krate(:,117)*n(:,patmo_idx_32H2SO4) &
        + 0.99*(krate(:,127)*n(:,patmo_idx_32CH3SCH3)*n(:,patmo_idx_O)) &
        + 0.99*(krate(:,131)*n(:,patmo_idx_32CH3SCH3)*n(:,patmo_idx_OH)) &
        + 0.99*(krate(:,135)*n(:,patmo_idx_32CH3SCH3)*n(:,patmo_idx_OH)) &
        + krate(:,139)*n(:,patmo_idx_32SO2_1) &
        + krate(:,143)*n(:,patmo_idx_32SO2_1) &
        + krate(:,159)*n(:,patmo_idx_32SO2_3) &
        + krate(:,163)*n(:,patmo_idx_32SO2_3) &
        + krate(:,194)*n(:,patmo_idx_32SO3) &
        + krate(:,289)*n(:,patmo_idx_OH)*n(:,patmo_idx_32SO3) &
        + krate(:,293)*n(:,patmo_idx_32SO3)*n(:,patmo_idx_O2) &
        + krate(:,313)*n(:,patmo_idx_32SO3) &
        + krate(:,317)*n(:,patmo_idx_32HSO3) &
        + krate(:,331)*n(:,patmo_idx_32SO4) 

        budget(:,patmo_idx_32SO2,2) = & 
        - krate(:,81)*n(:,patmo_idx_32SO2)*n(:,patmo_idx_HO2) &
        - krate(:,85)*n(:,patmo_idx_32SO2)*n(:,patmo_idx_O3) &
        - krate(:,105)*n(:,patmo_idx_32SO2)*n(:,patmo_idx_O) &
        - krate(:,109)*n(:,patmo_idx_32SO2)*n(:,patmo_idx_OH) &
        - krate(:,123)*n(:,patmo_idx_32SO2) &
        - krate(:,190)*n(:,patmo_idx_32SO2) &
        - krate(:,191)*n(:,patmo_idx_32SO2) &
        - krate(:,192)*n(:,patmo_idx_32SO2) &
        - krate(:,265)*n(:,patmo_idx_32SO2)*n(:,patmo_idx_O2) &
        - krate(:,269)*n(:,patmo_idx_32SO2)*n(:,patmo_idx_O) &
        - krate(:,273)*n(:,patmo_idx_32SO2)*n(:,patmo_idx_H) &
        - krate(:,297)*n(:,patmo_idx_32SO2)*n(:,patmo_idx_OH) &
        - krate(:,305)*n(:,patmo_idx_HO2)*n(:,patmo_idx_32SO2) &
        - krate(:,325)*n(:,patmo_idx_32SO2)*n(:,patmo_idx_OH)*n(:,patmo_idx_OH) &
        - krate(:,335)*n(:,patmo_idx_32SO2) &
        - krate(:,339)*n(:,patmo_idx_32SO2) &
        - krate(:,343)*n(:,patmo_idx_32SO2)*n(:,patmo_idx_CH4O3S) &
        - krate(:,347)*n(:,patmo_idx_32SO2)*n(:,patmo_idx_N2) &
        - krate(:,351)*n(:,patmo_idx_32SO2) &
        - krate(:,367)*n(:,patmo_idx_32SO2) &
        - krate(:,371)*n(:,patmo_idx_32SO2)*n(:,patmo_idx_N2) 
        
        budget(:,patmo_idx_32SO2,9) = 
        + krate(:,355)*n(:,patmo_idx_32SO3)*n(:,patmo_idx_32SO) &
        + krate(:,356)*n(:,patmo_idx_33SO3)*n(:,patmo_idx_32SO) &
        + krate(:,357)*n(:,patmo_idx_34SO3)*n(:,patmo_idx_32SO) &
        + krate(:,358)*n(:,patmo_idx_36SO3)*n(:,patmo_idx_32SO) &
        + krate(:,379)*n(:,patmo_idx_32SO3)*n(:,patmo_idx_32SO) &
        + krate(:,380)*n(:,patmo_idx_33SO3)*n(:,patmo_idx_32SO) &
        + krate(:,381)*n(:,patmo_idx_34SO3)*n(:,patmo_idx_32SO) &
        + krate(:,382)*n(:,patmo_idx_36SO3)*n(:,patmo_idx_32SO) &
        - krate(:,147)*n(:,patmo_idx_32SO2_1)*n(:,patmo_idx_32SO2) &
        - krate(:,148)*n(:,patmo_idx_33SO2_1)*n(:,patmo_idx_32SO2) &
        - krate(:,149)*n(:,patmo_idx_34SO2_1)*n(:,patmo_idx_32SO2) &
        - krate(:,150)*n(:,patmo_idx_36SO2_1)*n(:,patmo_idx_32SO2) &
        - krate(:,171)*n(:,patmo_idx_32SO2_3)*n(:,patmo_idx_32SO2) &
        - krate(:,172)*n(:,patmo_idx_33SO2_3)*n(:,patmo_idx_32SO2) &
        - krate(:,173)*n(:,patmo_idx_34SO2_3)*n(:,patmo_idx_32SO2) &
        - krate(:,174)*n(:,patmo_idx_36SO2_3)*n(:,patmo_idx_32SO2) 


    endif
!
    dn(:,patmo_idx_32SO2) = &
        + krate(:,57)*n(:,patmo_idx_32SO)*n(:,patmo_idx_O3) &
        + krate(:,61)*n(:,patmo_idx_32SO)*n(:,patmo_idx_O2) &
        + krate(:,65)*n(:,patmo_idx_32SO)*n(:,patmo_idx_OH) &
        - krate(:,81)*n(:,patmo_idx_32SO2)*n(:,patmo_idx_HO2) &
        - krate(:,85)*n(:,patmo_idx_32SO2)*n(:,patmo_idx_O3) &
        + krate(:,89)*n(:,patmo_idx_32HSO)*n(:,patmo_idx_O2) &
        + krate(:,97)*n(:,patmo_idx_32HSO2)*n(:,patmo_idx_O2) &
        - krate(:,105)*n(:,patmo_idx_32SO2)*n(:,patmo_idx_O) &
        - krate(:,109)*n(:,patmo_idx_32SO2)*n(:,patmo_idx_OH) &
        + krate(:,117)*n(:,patmo_idx_32H2SO4) &
        - krate(:,123)*n(:,patmo_idx_32SO2) &
        + 0.99*(krate(:,127)*n(:,patmo_idx_32CH3SCH3)*n(:,patmo_idx_O)) &
        + 0.99*(krate(:,131)*n(:,patmo_idx_32CH3SCH3)*n(:,patmo_idx_OH)) &
        + 0.99*(krate(:,135)*n(:,patmo_idx_32CH3SCH3)*n(:,patmo_idx_OH)) &
        + krate(:,139)*n(:,patmo_idx_32SO2_1) &
        + krate(:,143)*n(:,patmo_idx_32SO2_1) &
        - krate(:,147)*n(:,patmo_idx_32SO2_1)*n(:,patmo_idx_32SO2) &
        - krate(:,148)*n(:,patmo_idx_33SO2_1)*n(:,patmo_idx_32SO2) &
        - krate(:,149)*n(:,patmo_idx_34SO2_1)*n(:,patmo_idx_32SO2) &
        - krate(:,150)*n(:,patmo_idx_36SO2_1)*n(:,patmo_idx_32SO2) &
        + krate(:,159)*n(:,patmo_idx_32SO2_3) &
        + krate(:,163)*n(:,patmo_idx_32SO2_3) &
        - krate(:,171)*n(:,patmo_idx_32SO2_3)*n(:,patmo_idx_32SO2) &
        - krate(:,172)*n(:,patmo_idx_33SO2_3)*n(:,patmo_idx_32SO2) &
        - krate(:,173)*n(:,patmo_idx_34SO2_3)*n(:,patmo_idx_32SO2) &
        - krate(:,174)*n(:,patmo_idx_36SO2_3)*n(:,patmo_idx_32SO2) &
        - krate(:,190)*n(:,patmo_idx_32SO2) &
        - krate(:,191)*n(:,patmo_idx_32SO2) &
        - krate(:,192)*n(:,patmo_idx_32SO2) &
        + krate(:,194)*n(:,patmo_idx_32SO3) &
        - krate(:,265)*n(:,patmo_idx_32SO2)*n(:,patmo_idx_O2) &
        - krate(:,269)*n(:,patmo_idx_32SO2)*n(:,patmo_idx_O) &
        - krate(:,273)*n(:,patmo_idx_32SO2)*n(:,patmo_idx_H) &
        + krate(:,289)*n(:,patmo_idx_OH)*n(:,patmo_idx_32SO3) &
        + krate(:,293)*n(:,patmo_idx_32SO3)*n(:,patmo_idx_O2) &
        - krate(:,297)*n(:,patmo_idx_32SO2)*n(:,patmo_idx_OH) &
        - krate(:,305)*n(:,patmo_idx_HO2)*n(:,patmo_idx_32SO2) &
        + krate(:,313)*n(:,patmo_idx_32SO3) &
        + krate(:,317)*n(:,patmo_idx_32HSO3) &
        - krate(:,325)*n(:,patmo_idx_32SO2)*n(:,patmo_idx_OH)*n(:,patmo_idx_OH) &
        + krate(:,331)*n(:,patmo_idx_32SO4) &
        - krate(:,335)*n(:,patmo_idx_32SO2) &
        - krate(:,339)*n(:,patmo_idx_32SO2) &
        - krate(:,343)*n(:,patmo_idx_32SO2)*n(:,patmo_idx_CH4O3S) &
        - krate(:,347)*n(:,patmo_idx_32SO2)*n(:,patmo_idx_N2) &
        - krate(:,351)*n(:,patmo_idx_32SO2) &
        + krate(:,355)*n(:,patmo_idx_32SO3)*n(:,patmo_idx_32SO) &
        + krate(:,356)*n(:,patmo_idx_33SO3)*n(:,patmo_idx_32SO) &
        + krate(:,357)*n(:,patmo_idx_34SO3)*n(:,patmo_idx_32SO) &
        + krate(:,358)*n(:,patmo_idx_36SO3)*n(:,patmo_idx_32SO) &
        - krate(:,367)*n(:,patmo_idx_32SO2) &
        - krate(:,371)*n(:,patmo_idx_32SO2)*n(:,patmo_idx_N2) &
        + krate(:,379)*n(:,patmo_idx_32SO3)*n(:,patmo_idx_32SO) &
        + krate(:,380)*n(:,patmo_idx_33SO3)*n(:,patmo_idx_32SO) &
        + krate(:,381)*n(:,patmo_idx_34SO3)*n(:,patmo_idx_32SO) &
        + krate(:,382)*n(:,patmo_idx_36SO3)*n(:,patmo_idx_32SO)

    dn(:,patmo_idx_33HSO3) = &
        - krate(:,102)*n(:,patmo_idx_33HSO3)*n(:,patmo_idx_O2) &
        + krate(:,110)*n(:,patmo_idx_33SO2)*n(:,patmo_idx_OH) &
        + krate(:,310)*n(:,patmo_idx_HO2)*n(:,patmo_idx_33SO3) &
        - krate(:,318)*n(:,patmo_idx_33HSO3)

    if(end_of_run) then
            budget(:,patmo_idx_34CH3SCH3,1) = & 
                      + krate(:,337)*n(:,patmo_idx_34SO2) &
        + krate(:,341)*n(:,patmo_idx_34SO2) &
        + krate(:,345)*n(:,patmo_idx_34SO2)*n(:,patmo_idx_CH4O3S) 
            budget(:,patmo_idx_34CH3SCH3,2) = & 
                - krate(:,129)*n(:,patmo_idx_34CH3SCH3)*n(:,patmo_idx_O) &
        - krate(:,133)*n(:,patmo_idx_34CH3SCH3)*n(:,patmo_idx_OH) &
        - krate(:,137)*n(:,patmo_idx_34CH3SCH3)*n(:,patmo_idx_OH) 
    endif


    dn(:,patmo_idx_34CH3SCH3) = &
        - krate(:,129)*n(:,patmo_idx_34CH3SCH3)*n(:,patmo_idx_O) &
        - krate(:,133)*n(:,patmo_idx_34CH3SCH3)*n(:,patmo_idx_OH) &
        - krate(:,137)*n(:,patmo_idx_34CH3SCH3)*n(:,patmo_idx_OH) &
        + krate(:,337)*n(:,patmo_idx_34SO2) &
        + krate(:,341)*n(:,patmo_idx_34SO2) &
        + krate(:,345)*n(:,patmo_idx_34SO2)*n(:,patmo_idx_CH4O3S) 
       ! - (7.940d-6)*n(:,patmo_idx_34CH3SCH3)*n(:,patmo_idx_OH)

    if (end_of_run) then
        budget(:,patmo_idx_32SO4,1) = + krate(:,123)*n(:,patmo_idx_32SO2)

        budget(:,patmo_idx_32SO4,2) = - krate(:,331)*n(:,patmo_idx_32SO4)
    endif
!
    dn(:,patmo_idx_32SO4) = &
        + krate(:,123)*n(:,patmo_idx_32SO2) &
        - krate(:,331)*n(:,patmo_idx_32SO4)

    dn(:,patmo_idx_CH4O3S) = 0d0
        ! + krate(:,135)*n(:,patmo_idx_32CH3SCH3)*n(:,patmo_idx_OH) &
        ! + krate(:,136)*n(:,patmo_idx_33CH3SCH3)*n(:,patmo_idx_OH) &
        ! + krate(:,137)*n(:,patmo_idx_34CH3SCH3)*n(:,patmo_idx_OH) &
        ! + krate(:,138)*n(:,patmo_idx_36CH3SCH3)*n(:,patmo_idx_OH) &
        ! - krate(:,343)*n(:,patmo_idx_32SO2)*n(:,patmo_idx_CH4O3S) &
        ! - krate(:,344)*n(:,patmo_idx_33SO2)*n(:,patmo_idx_CH4O3S) &
        ! - krate(:,345)*n(:,patmo_idx_34SO2)*n(:,patmo_idx_CH4O3S) &
        ! - krate(:,346)*n(:,patmo_idx_36SO2)*n(:,patmo_idx_CH4O3S)

    dn(:,patmo_idx_36SH) = &
        + krate(:,4)*n(:,patmo_idx_36COS)*n(:,patmo_idx_OH) &
        + krate(:,12)*n(:,patmo_idx_36CS2)*n(:,patmo_idx_OH) &
        + krate(:,32)*n(:,patmo_idx_36H2S)*n(:,patmo_idx_OH) &
        + krate(:,36)*n(:,patmo_idx_36H2S)*n(:,patmo_idx_O) &
        + krate(:,40)*n(:,patmo_idx_36H2S)*n(:,patmo_idx_H) &
        - krate(:,48)*n(:,patmo_idx_36SH)*n(:,patmo_idx_O) &
        - krate(:,52)*n(:,patmo_idx_36SH)*n(:,patmo_idx_O2) &
        - krate(:,56)*n(:,patmo_idx_36SH)*n(:,patmo_idx_O3) &
        + krate(:,96)*n(:,patmo_idx_36HSO)*n(:,patmo_idx_O3) &
        + krate(:,187)*n(:,patmo_idx_36H2S) &
        - krate(:,212)*n(:,patmo_idx_CO2)*n(:,patmo_idx_36SH) &
        - krate(:,220)*n(:,patmo_idx_36SH)*n(:,patmo_idx_36COS) &
        - krate(:,240)*n(:,patmo_idx_H2O)*n(:,patmo_idx_36SH) &
        - krate(:,244)*n(:,patmo_idx_OH)*n(:,patmo_idx_36SH) &
        - krate(:,248)*n(:,patmo_idx_H2)*n(:,patmo_idx_36SH) &
        + krate(:,256)*n(:,patmo_idx_H)*n(:,patmo_idx_36SO) &
        + krate(:,260)*n(:,patmo_idx_OH)*n(:,patmo_idx_36SO) &
        + krate(:,264)*n(:,patmo_idx_36HSO)*n(:,patmo_idx_O2) &
        - krate(:,304)*n(:,patmo_idx_O2)*n(:,patmo_idx_O2)*n(:,patmo_idx_36SH)

    dn(:,patmo_idx_36SO3) = &
        + krate(:,84)*n(:,patmo_idx_36SO2)*n(:,patmo_idx_HO2) &
        + krate(:,88)*n(:,patmo_idx_36SO2)*n(:,patmo_idx_O3) &
        + krate(:,104)*n(:,patmo_idx_36HSO3)*n(:,patmo_idx_O2) &
        + krate(:,108)*n(:,patmo_idx_36SO2)*n(:,patmo_idx_O) &
        - krate(:,116)*n(:,patmo_idx_36SO3)*n(:,patmo_idx_H2O) &
        + krate(:,150)*n(:,patmo_idx_36SO2_1)*n(:,patmo_idx_32SO2) &
        + krate(:,170)*n(:,patmo_idx_36SO2_3)*n(:,patmo_idx_O2) &
        + krate(:,174)*n(:,patmo_idx_36SO2_3)*n(:,patmo_idx_32SO2) &
        - krate(:,178)*n(:,patmo_idx_36SO3) &
        - krate(:,292)*n(:,patmo_idx_OH)*n(:,patmo_idx_36SO3) &
        - krate(:,296)*n(:,patmo_idx_36SO3)*n(:,patmo_idx_O2) &
        - krate(:,312)*n(:,patmo_idx_HO2)*n(:,patmo_idx_36SO3) &
        - krate(:,316)*n(:,patmo_idx_36SO3) &
        + krate(:,324)*n(:,patmo_idx_36H2SO4) &
        - krate(:,358)*n(:,patmo_idx_36SO3)*n(:,patmo_idx_32SO) &
        - krate(:,378)*n(:,patmo_idx_36SO3)*n(:,patmo_idx_O) &
        - krate(:,382)*n(:,patmo_idx_36SO3)*n(:,patmo_idx_32SO)

    dn(:,patmo_idx_36HSO) = &
        + krate(:,44)*n(:,patmo_idx_36H2S)*n(:,patmo_idx_HO2) &
        + krate(:,56)*n(:,patmo_idx_36SH)*n(:,patmo_idx_O3) &
        - krate(:,92)*n(:,patmo_idx_36HSO)*n(:,patmo_idx_O2) &
        - krate(:,96)*n(:,patmo_idx_36HSO)*n(:,patmo_idx_O3) &
        - krate(:,252)*n(:,patmo_idx_H2O)*n(:,patmo_idx_36HSO) &
        - krate(:,264)*n(:,patmo_idx_36HSO)*n(:,patmo_idx_O2) &
        + krate(:,300)*n(:,patmo_idx_36SO2)*n(:,patmo_idx_OH) &
        + krate(:,304)*n(:,patmo_idx_O2)*n(:,patmo_idx_O2)*n(:,patmo_idx_36SH)

    dn(:,patmo_idx_32CS) = &
        + krate(:,13)*n(:,patmo_idx_32CS2)*n(:,patmo_idx_O) &
        - krate(:,17)*n(:,patmo_idx_32CS)*n(:,patmo_idx_O2) &
        - krate(:,21)*n(:,patmo_idx_32CS)*n(:,patmo_idx_O3) &
        - krate(:,25)*n(:,patmo_idx_32CS)*n(:,patmo_idx_O) &
        + krate(:,202)*n(:,patmo_idx_32CS2) &
        - krate(:,221)*n(:,patmo_idx_32CS)*n(:,patmo_idx_32SO) &
        + krate(:,225)*n(:,patmo_idx_32COS)*n(:,patmo_idx_O) &
        + krate(:,229)*n(:,patmo_idx_32COS)*n(:,patmo_idx_O2) &
        + krate(:,233)*n(:,patmo_idx_CO)*n(:,patmo_idx_32S)

    dn(:,patmo_idx_34SO3) = &
        + krate(:,83)*n(:,patmo_idx_34SO2)*n(:,patmo_idx_HO2) &
        + krate(:,87)*n(:,patmo_idx_34SO2)*n(:,patmo_idx_O3) &
        + krate(:,103)*n(:,patmo_idx_34HSO3)*n(:,patmo_idx_O2) &
        + krate(:,107)*n(:,patmo_idx_34SO2)*n(:,patmo_idx_O) &
        - krate(:,115)*n(:,patmo_idx_34SO3)*n(:,patmo_idx_H2O) &
        + krate(:,149)*n(:,patmo_idx_34SO2_1)*n(:,patmo_idx_32SO2) &
        + krate(:,169)*n(:,patmo_idx_34SO2_3)*n(:,patmo_idx_O2) &
        + krate(:,173)*n(:,patmo_idx_34SO2_3)*n(:,patmo_idx_32SO2) &
        - krate(:,208)*n(:,patmo_idx_34SO3) &
        - krate(:,291)*n(:,patmo_idx_OH)*n(:,patmo_idx_34SO3) &
        - krate(:,295)*n(:,patmo_idx_34SO3)*n(:,patmo_idx_O2) &
        - krate(:,311)*n(:,patmo_idx_HO2)*n(:,patmo_idx_34SO3) &
        - krate(:,315)*n(:,patmo_idx_34SO3) &
        + krate(:,323)*n(:,patmo_idx_34H2SO4) &
        - krate(:,357)*n(:,patmo_idx_34SO3)*n(:,patmo_idx_32SO) &
        - krate(:,377)*n(:,patmo_idx_34SO3)*n(:,patmo_idx_O) &
        - krate(:,381)*n(:,patmo_idx_34SO3)*n(:,patmo_idx_32SO)

    dn(:,patmo_idx_H2) = 0d0
        ! + krate(:,37)*n(:,patmo_idx_32H2S)*n(:,patmo_idx_H) &
        ! + krate(:,38)*n(:,patmo_idx_33H2S)*n(:,patmo_idx_H) &
        ! + krate(:,39)*n(:,patmo_idx_34H2S)*n(:,patmo_idx_H) &
        ! + krate(:,40)*n(:,patmo_idx_36H2S)*n(:,patmo_idx_H) &
        ! - krate(:,245)*n(:,patmo_idx_H2)*n(:,patmo_idx_32SH) &
        ! - krate(:,246)*n(:,patmo_idx_H2)*n(:,patmo_idx_33SH) &
        ! - krate(:,247)*n(:,patmo_idx_H2)*n(:,patmo_idx_34SH) &
        ! - krate(:,248)*n(:,patmo_idx_H2)*n(:,patmo_idx_36SH)

    dn(:,patmo_idx_36CS2) = &
        - krate(:,12)*n(:,patmo_idx_36CS2)*n(:,patmo_idx_OH) &
        - krate(:,16)*n(:,patmo_idx_36CS2)*n(:,patmo_idx_O) &
        - krate(:,205)*n(:,patmo_idx_36CS2) &
        + krate(:,220)*n(:,patmo_idx_36SH)*n(:,patmo_idx_36COS) &
        + krate(:,224)*n(:,patmo_idx_36CS)*n(:,patmo_idx_36SO)

    dn(:,patmo_idx_34HSO) = &
        + krate(:,43)*n(:,patmo_idx_34H2S)*n(:,patmo_idx_HO2) &
        + krate(:,55)*n(:,patmo_idx_34SH)*n(:,patmo_idx_O3) &
        - krate(:,91)*n(:,patmo_idx_34HSO)*n(:,patmo_idx_O2) &
        - krate(:,95)*n(:,patmo_idx_34HSO)*n(:,patmo_idx_O3) &
        - krate(:,251)*n(:,patmo_idx_H2O)*n(:,patmo_idx_34HSO) &
        - krate(:,263)*n(:,patmo_idx_34HSO)*n(:,patmo_idx_O2) &
        + krate(:,299)*n(:,patmo_idx_34SO2)*n(:,patmo_idx_OH) &
        + krate(:,303)*n(:,patmo_idx_O2)*n(:,patmo_idx_O2)*n(:,patmo_idx_34SH)

    dn(:,patmo_idx_33HSO) = &
        + krate(:,42)*n(:,patmo_idx_33H2S)*n(:,patmo_idx_HO2) &
        + krate(:,54)*n(:,patmo_idx_33SH)*n(:,patmo_idx_O3) &
        - krate(:,90)*n(:,patmo_idx_33HSO)*n(:,patmo_idx_O2) &
        - krate(:,94)*n(:,patmo_idx_33HSO)*n(:,patmo_idx_O3) &
        - krate(:,250)*n(:,patmo_idx_H2O)*n(:,patmo_idx_33HSO) &
        - krate(:,262)*n(:,patmo_idx_33HSO)*n(:,patmo_idx_O2) &
        + krate(:,298)*n(:,patmo_idx_33SO2)*n(:,patmo_idx_OH) &
        + krate(:,302)*n(:,patmo_idx_O2)*n(:,patmo_idx_O2)*n(:,patmo_idx_33SH)

    dn(:,patmo_idx_33CS2) = &
        - krate(:,10)*n(:,patmo_idx_33CS2)*n(:,patmo_idx_OH) &
        - krate(:,14)*n(:,patmo_idx_33CS2)*n(:,patmo_idx_O) &
        - krate(:,183)*n(:,patmo_idx_33CS2) &
        + krate(:,218)*n(:,patmo_idx_33SH)*n(:,patmo_idx_33COS) &
        + krate(:,222)*n(:,patmo_idx_33CS)*n(:,patmo_idx_33SO)

    !if (end_of_run) then
    !    budget(:,patmo_idx_36COS,1) = &
    !    + (0.83*(krate(:,12)*n(:,patmo_idx_36CS2)*n(:,patmo_idx_OH))) &
    !    + krate(:,20)*n(:,patmo_idx_36CS)*n(:,patmo_idx_O2) &
    !    + krate(:,24)*n(:,patmo_idx_36CS)*n(:,patmo_idx_O3) &
    !    + krate(:,212)*n(:,patmo_idx_CO2)*n(:,patmo_idx_36SH) &
    !    + krate(:,216)*n(:,patmo_idx_CO)*n(:,patmo_idx_36SO) &
    !    + 0.007*(krate(:,134))*n(:,patmo_idx_36CH3SCH3)*n(:,patmo_idx_OH)

    !    budget(:,patmo_idx_36COS,2) = &
    !    + krate(:,4)*n(:,patmo_idx_36COS)*n(:,patmo_idx_OH) &
    !    + krate(:,8)*n(:,patmo_idx_36COS)*n(:,patmo_idx_O) &
    !    + krate(:,207)*n(:,patmo_idx_36COS) &
    !    + krate(:,220)*n(:,patmo_idx_36SH)*n(:,patmo_idx_36COS) &
    !    + krate(:,228)*n(:,patmo_idx_36COS)*n(:,patmo_idx_O) &
    !    + krate(:,232)*n(:,patmo_idx_36COS)*n(:,patmo_idx_O2)
    !endif

    dn(:,patmo_idx_36COS) = &
        - krate(:,4)*n(:,patmo_idx_36COS)*n(:,patmo_idx_OH) &
        - krate(:,8)*n(:,patmo_idx_36COS)*n(:,patmo_idx_O) &
        + (0.83*(krate(:,12)*n(:,patmo_idx_36CS2)*n(:,patmo_idx_OH))) &
        + krate(:,20)*n(:,patmo_idx_36CS)*n(:,patmo_idx_O2) &
        + krate(:,24)*n(:,patmo_idx_36CS)*n(:,patmo_idx_O3) &
        - krate(:,207)*n(:,patmo_idx_36COS) &
        + krate(:,212)*n(:,patmo_idx_CO2)*n(:,patmo_idx_36SH) &
        + krate(:,216)*n(:,patmo_idx_CO)*n(:,patmo_idx_36SO) &
        - krate(:,220)*n(:,patmo_idx_36SH)*n(:,patmo_idx_36COS) &
        - krate(:,228)*n(:,patmo_idx_36COS)*n(:,patmo_idx_O) &
        - krate(:,232)*n(:,patmo_idx_36COS)*n(:,patmo_idx_O2) &
        + 0.007*(krate(:,134)*n(:,patmo_idx_36CH3SCH3)*n(:,patmo_idx_OH)) &
        + 0.007*(krate(:,138)*n(:,patmo_idx_36CH3SCH3)*n(:,patmo_idx_OH)) &
        + 0.007*(krate(:,130)*n(:,patmo_idx_36CH3SCH3)*n(:,patmo_idx_O))
    if (end_of_run) then
        budget(:,patmo_idx_34SO4,1) = + krate(:,125)*n(:,patmo_idx_34SO2)
!
        budget(:,patmo_idx_34SO4,2) = - krate(:,333)*n(:,patmo_idx_34SO4)
    endif
!
    dn(:,patmo_idx_34SO4) = &
        + krate(:,125)*n(:,patmo_idx_34SO2) &
        - krate(:,333)*n(:,patmo_idx_34SO4)

    dn(:,patmo_idx_32SH) = &
        + krate(:,1)*n(:,patmo_idx_32COS)*n(:,patmo_idx_OH) &
        + krate(:,9)*n(:,patmo_idx_32CS2)*n(:,patmo_idx_OH) &
        + krate(:,29)*n(:,patmo_idx_32H2S)*n(:,patmo_idx_OH) &
        + krate(:,33)*n(:,patmo_idx_32H2S)*n(:,patmo_idx_O) &
        + krate(:,37)*n(:,patmo_idx_32H2S)*n(:,patmo_idx_H) &
        - krate(:,45)*n(:,patmo_idx_32SH)*n(:,patmo_idx_O) &
        - krate(:,49)*n(:,patmo_idx_32SH)*n(:,patmo_idx_O2) &
        - krate(:,53)*n(:,patmo_idx_32SH)*n(:,patmo_idx_O3) &
        + krate(:,93)*n(:,patmo_idx_32HSO)*n(:,patmo_idx_O3) &
        + krate(:,206)*n(:,patmo_idx_32H2S) &
        - krate(:,209)*n(:,patmo_idx_CO2)*n(:,patmo_idx_32SH) &
        - krate(:,217)*n(:,patmo_idx_32SH)*n(:,patmo_idx_32COS) &
        - krate(:,237)*n(:,patmo_idx_H2O)*n(:,patmo_idx_32SH) &
        - krate(:,241)*n(:,patmo_idx_OH)*n(:,patmo_idx_32SH) &
        - krate(:,245)*n(:,patmo_idx_H2)*n(:,patmo_idx_32SH) &
        + krate(:,253)*n(:,patmo_idx_H)*n(:,patmo_idx_32SO) &
        + krate(:,257)*n(:,patmo_idx_OH)*n(:,patmo_idx_32SO) &
        + krate(:,261)*n(:,patmo_idx_32HSO)*n(:,patmo_idx_O2) &
        - krate(:,301)*n(:,patmo_idx_O2)*n(:,patmo_idx_O2)*n(:,patmo_idx_32SH)
 
    if (end_of_run) then
        budget(:,patmo_idx_32COS,1) = &
        + (0.83*(krate(:,9)*n(:,patmo_idx_32CS2)*n(:,patmo_idx_OH))) &
        + krate(:,17)*n(:,patmo_idx_32CS)*n(:,patmo_idx_O2) &
        + krate(:,21)*n(:,patmo_idx_32CS)*n(:,patmo_idx_O3) &
        + krate(:,209)*n(:,patmo_idx_CO2)*n(:,patmo_idx_32SH) &
        + krate(:,213)*n(:,patmo_idx_CO)*n(:,patmo_idx_32SO) &
        + 0.007*(krate(:,131))*n(:,patmo_idx_32CH3SCH3)*n(:,patmo_idx_OH) &
        + 0.007*(krate(:,127)*n(:,patmo_idx_32CH3SCH3)*n(:,patmo_idx_O)) &
        + 0.007*(krate(:,135))*n(:,patmo_idx_32CH3SCH3)*n(:,patmo_idx_OH)
        
        budget(:,patmo_idx_32COS,2) = &
        - krate(:,1)*n(:,patmo_idx_32COS)*n(:,patmo_idx_OH) &
        - krate(:,5)*n(:,patmo_idx_32COS)*n(:,patmo_idx_O) &
        - krate(:,179)*n(:,patmo_idx_32COS) &
        - krate(:,217)*n(:,patmo_idx_32SH)*n(:,patmo_idx_32COS) &
        - krate(:,225)*n(:,patmo_idx_32COS)*n(:,patmo_idx_O) &
        - krate(:,229)*n(:,patmo_idx_32COS)*n(:,patmo_idx_O2) 

        budget(:,patmo_idx_32COS,9) = &
        + krate(:,179)*n(:,patmo_idx_32COS) 
    endif 
    
    dn(:,patmo_idx_32COS) = &
        - krate(:,1)*n(:,patmo_idx_32COS)*n(:,patmo_idx_OH) &
        - krate(:,5)*n(:,patmo_idx_32COS)*n(:,patmo_idx_O) &
        + (0.83*(krate(:,9)*n(:,patmo_idx_32CS2)*n(:,patmo_idx_OH))) &
        + krate(:,17)*n(:,patmo_idx_32CS)*n(:,patmo_idx_O2) &
        + krate(:,21)*n(:,patmo_idx_32CS)*n(:,patmo_idx_O3) &
        - krate(:,179)*n(:,patmo_idx_32COS) &
        + krate(:,209)*n(:,patmo_idx_CO2)*n(:,patmo_idx_32SH) &
        + krate(:,213)*n(:,patmo_idx_CO)*n(:,patmo_idx_32SO) &
        - krate(:,217)*n(:,patmo_idx_32SH)*n(:,patmo_idx_32COS) &
        - krate(:,225)*n(:,patmo_idx_32COS)*n(:,patmo_idx_O) &
        - krate(:,229)*n(:,patmo_idx_32COS)*n(:,patmo_idx_O2) &
        + 0.007*(krate(:,131))*n(:,patmo_idx_32CH3SCH3)*n(:,patmo_idx_OH) &
        + 0.007*(krate(:,127)*n(:,patmo_idx_32CH3SCH3)*n(:,patmo_idx_O)) &
        + 0.007*(krate(:,135))*n(:,patmo_idx_32CH3SCH3)*n(:,patmo_idx_OH) 
        

    dn(:,patmo_idx_O3) = 0d0
        ! - krate(:,21)*n(:,patmo_idx_32CS)*n(:,patmo_idx_O3) &
        ! - krate(:,22)*n(:,patmo_idx_33CS)*n(:,patmo_idx_O3) &
        ! - krate(:,23)*n(:,patmo_idx_34CS)*n(:,patmo_idx_O3) &
        ! - krate(:,24)*n(:,patmo_idx_36CS)*n(:,patmo_idx_O3) &
        ! - krate(:,53)*n(:,patmo_idx_32SH)*n(:,patmo_idx_O3) &
        ! - krate(:,54)*n(:,patmo_idx_33SH)*n(:,patmo_idx_O3) &
        ! - krate(:,55)*n(:,patmo_idx_34SH)*n(:,patmo_idx_O3) &
        ! - krate(:,56)*n(:,patmo_idx_36SH)*n(:,patmo_idx_O3) &
        ! - krate(:,57)*n(:,patmo_idx_32SO)*n(:,patmo_idx_O3) &
        ! - krate(:,58)*n(:,patmo_idx_33SO)*n(:,patmo_idx_O3) &
        ! - krate(:,59)*n(:,patmo_idx_34SO)*n(:,patmo_idx_O3) &
        ! - krate(:,60)*n(:,patmo_idx_36SO)*n(:,patmo_idx_O3) &
        ! - krate(:,73)*n(:,patmo_idx_32S)*n(:,patmo_idx_O3) &
        ! - krate(:,74)*n(:,patmo_idx_33S)*n(:,patmo_idx_O3) &
        ! - krate(:,75)*n(:,patmo_idx_34S)*n(:,patmo_idx_O3) &
        ! - krate(:,76)*n(:,patmo_idx_36S)*n(:,patmo_idx_O3) &
        ! - krate(:,85)*n(:,patmo_idx_32SO2)*n(:,patmo_idx_O3) &
        ! - krate(:,86)*n(:,patmo_idx_33SO2)*n(:,patmo_idx_O3) &
        ! - krate(:,87)*n(:,patmo_idx_34SO2)*n(:,patmo_idx_O3) &
        ! - krate(:,88)*n(:,patmo_idx_36SO2)*n(:,patmo_idx_O3) &
        ! - krate(:,93)*n(:,patmo_idx_32HSO)*n(:,patmo_idx_O3) &
        ! - krate(:,94)*n(:,patmo_idx_33HSO)*n(:,patmo_idx_O3) &
        ! - krate(:,95)*n(:,patmo_idx_34HSO)*n(:,patmo_idx_O3) &
        ! - krate(:,96)*n(:,patmo_idx_36HSO)*n(:,patmo_idx_O3) &
        ! + krate(:,121)*n(:,patmo_idx_O)*n(:,patmo_idx_O2) &
        ! - krate(:,184)*n(:,patmo_idx_O3) &
        ! + krate(:,229)*n(:,patmo_idx_32COS)*n(:,patmo_idx_O2) &
        ! + krate(:,230)*n(:,patmo_idx_33COS)*n(:,patmo_idx_O2) &
        ! + krate(:,231)*n(:,patmo_idx_34COS)*n(:,patmo_idx_O2) &
        ! + krate(:,232)*n(:,patmo_idx_36COS)*n(:,patmo_idx_O2) &
        ! + krate(:,261)*n(:,patmo_idx_32HSO)*n(:,patmo_idx_O2) &
        ! + krate(:,262)*n(:,patmo_idx_33HSO)*n(:,patmo_idx_O2) &
        ! + krate(:,263)*n(:,patmo_idx_34HSO)*n(:,patmo_idx_O2) &
        ! + krate(:,264)*n(:,patmo_idx_36HSO)*n(:,patmo_idx_O2) &
        ! + krate(:,265)*n(:,patmo_idx_32SO2)*n(:,patmo_idx_O2) &
        ! + krate(:,266)*n(:,patmo_idx_33SO2)*n(:,patmo_idx_O2) &
        ! + krate(:,267)*n(:,patmo_idx_34SO2)*n(:,patmo_idx_O2) &
        ! + krate(:,268)*n(:,patmo_idx_36SO2)*n(:,patmo_idx_O2) &
        ! + krate(:,281)*n(:,patmo_idx_O2)*n(:,patmo_idx_32SO) &
        ! + krate(:,282)*n(:,patmo_idx_O2)*n(:,patmo_idx_33SO) &
        ! + krate(:,283)*n(:,patmo_idx_O2)*n(:,patmo_idx_34SO) &
        ! + krate(:,284)*n(:,patmo_idx_O2)*n(:,patmo_idx_36SO) &
        ! + krate(:,293)*n(:,patmo_idx_32SO3)*n(:,patmo_idx_O2) &
        ! + krate(:,294)*n(:,patmo_idx_33SO3)*n(:,patmo_idx_O2) &
        ! + krate(:,295)*n(:,patmo_idx_34SO3)*n(:,patmo_idx_O2) &
        ! + krate(:,296)*n(:,patmo_idx_36SO3)*n(:,patmo_idx_O2) &
        ! + krate(:,301)*n(:,patmo_idx_O2)*n(:,patmo_idx_O2)*n(:,patmo_idx_32SH) &
        ! + krate(:,302)*n(:,patmo_idx_O2)*n(:,patmo_idx_O2)*n(:,patmo_idx_33SH) &
        ! + krate(:,303)*n(:,patmo_idx_O2)*n(:,patmo_idx_O2)*n(:,patmo_idx_34SH) &
        ! + krate(:,304)*n(:,patmo_idx_O2)*n(:,patmo_idx_O2)*n(:,patmo_idx_36SH) &
        ! - krate(:,329)*n(:,patmo_idx_O3)

    dn(:,patmo_idx_32SO3) = &
        + krate(:,81)*n(:,patmo_idx_32SO2)*n(:,patmo_idx_HO2) &
        + krate(:,85)*n(:,patmo_idx_32SO2)*n(:,patmo_idx_O3) &
        + krate(:,101)*n(:,patmo_idx_32HSO3)*n(:,patmo_idx_O2) &
        + krate(:,105)*n(:,patmo_idx_32SO2)*n(:,patmo_idx_O) &
        - krate(:,113)*n(:,patmo_idx_32SO3)*n(:,patmo_idx_H2O) &
        + krate(:,147)*n(:,patmo_idx_32SO2_1)*n(:,patmo_idx_32SO2) &
        + krate(:,167)*n(:,patmo_idx_32SO2_3)*n(:,patmo_idx_O2) &
        + krate(:,171)*n(:,patmo_idx_32SO2_3)*n(:,patmo_idx_32SO2) &
        - krate(:,194)*n(:,patmo_idx_32SO3) &
        - krate(:,289)*n(:,patmo_idx_OH)*n(:,patmo_idx_32SO3) &
        - krate(:,293)*n(:,patmo_idx_32SO3)*n(:,patmo_idx_O2) &
        - krate(:,309)*n(:,patmo_idx_HO2)*n(:,patmo_idx_32SO3) &
        - krate(:,313)*n(:,patmo_idx_32SO3) &
        + krate(:,321)*n(:,patmo_idx_32H2SO4) &
        - krate(:,355)*n(:,patmo_idx_32SO3)*n(:,patmo_idx_32SO) &
        - krate(:,375)*n(:,patmo_idx_32SO3)*n(:,patmo_idx_O) &
        - krate(:,379)*n(:,patmo_idx_32SO3)*n(:,patmo_idx_32SO)

    dn(:,patmo_idx_36SO2_1) = &
        - krate(:,142)*n(:,patmo_idx_36SO2_1) &
        - krate(:,146)*n(:,patmo_idx_36SO2_1) &
        - krate(:,150)*n(:,patmo_idx_36SO2_1)*n(:,patmo_idx_32SO2) &
        - krate(:,154)*n(:,patmo_idx_36SO2_1) &
        - krate(:,158)*n(:,patmo_idx_36SO2_1) &
        + krate(:,176)*n(:,patmo_idx_36SO2) &
        + krate(:,350)*n(:,patmo_idx_36SO2)*n(:,patmo_idx_N2) &
        + krate(:,354)*n(:,patmo_idx_36SO2) &
        + krate(:,358)*n(:,patmo_idx_36SO3)*n(:,patmo_idx_32SO) &
        + krate(:,362)*n(:,patmo_idx_36SO2_3)*n(:,patmo_idx_N2) &
        + krate(:,366)*n(:,patmo_idx_36SO2_3)

    if (end_of_run) then
        budget(:,patmo_idx_34SO2,1) = &
        + krate(:,59)*n(:,patmo_idx_34SO)*n(:,patmo_idx_O3) &
        + krate(:,63)*n(:,patmo_idx_34SO)*n(:,patmo_idx_O2) &
        + krate(:,67)*n(:,patmo_idx_34SO)*n(:,patmo_idx_OH) &
        + krate(:,91)*n(:,patmo_idx_34HSO)*n(:,patmo_idx_O2) &
        + krate(:,99)*n(:,patmo_idx_34HSO2)*n(:,patmo_idx_O2) &
        + krate(:,119)*n(:,patmo_idx_34H2SO4) &
        + 0.99*(krate(:,129)*n(:,patmo_idx_34CH3SCH3)*n(:,patmo_idx_O)) &
        + 0.99*(krate(:,133)*n(:,patmo_idx_34CH3SCH3)*n(:,patmo_idx_OH)) &
        + 0.99*(krate(:,137)*n(:,patmo_idx_34CH3SCH3)*n(:,patmo_idx_OH)) &
        + krate(:,141)*n(:,patmo_idx_34SO2_1) &
        + krate(:,145)*n(:,patmo_idx_34SO2_1) &
        + krate(:,161)*n(:,patmo_idx_34SO2_3) &
        + krate(:,165)*n(:,patmo_idx_34SO2_3) &
        + krate(:,208)*n(:,patmo_idx_34SO3) &
        + krate(:,291)*n(:,patmo_idx_OH)*n(:,patmo_idx_34SO3) &
        + krate(:,295)*n(:,patmo_idx_34SO3)*n(:,patmo_idx_O2) &
        + krate(:,315)*n(:,patmo_idx_34SO3) &
        + krate(:,319)*n(:,patmo_idx_34HSO3) &
        + krate(:,333)*n(:,patmo_idx_34SO4) &
        + krate(:,333)*n(:,patmo_idx_34SO4) 
        



        budget(:,patmo_idx_34SO2,2) = &
        - krate(:,83)*n(:,patmo_idx_34SO2)*n(:,patmo_idx_HO2) &
        - krate(:,87)*n(:,patmo_idx_34SO2)*n(:,patmo_idx_O3) &
        - krate(:,107)*n(:,patmo_idx_34SO2)*n(:,patmo_idx_O) &
        - krate(:,111)*n(:,patmo_idx_34SO2)*n(:,patmo_idx_OH) &
        - krate(:,125)*n(:,patmo_idx_34SO2) &
        - krate(:,180)*n(:,patmo_idx_34SO2) &
        - krate(:,181)*n(:,patmo_idx_34SO2) &
        - krate(:,182)*n(:,patmo_idx_34SO2) &
        - krate(:,267)*n(:,patmo_idx_34SO2)*n(:,patmo_idx_O2) &
        - krate(:,271)*n(:,patmo_idx_34SO2)*n(:,patmo_idx_O) &
        - krate(:,275)*n(:,patmo_idx_34SO2)*n(:,patmo_idx_H) &
        - krate(:,299)*n(:,patmo_idx_34SO2)*n(:,patmo_idx_OH) &
        - krate(:,307)*n(:,patmo_idx_HO2)*n(:,patmo_idx_34SO2) &
        - krate(:,327)*n(:,patmo_idx_34SO2)*n(:,patmo_idx_OH)*n(:,patmo_idx_OH) &
        - krate(:,337)*n(:,patmo_idx_34SO2) &
        - krate(:,341)*n(:,patmo_idx_34SO2) &
        - krate(:,345)*n(:,patmo_idx_34SO2)*n(:,patmo_idx_CH4O3S) &
        - krate(:,349)*n(:,patmo_idx_34SO2)*n(:,patmo_idx_N2) &
        - krate(:,353)*n(:,patmo_idx_34SO2) &
        - krate(:,369)*n(:,patmo_idx_34SO2) &
        - krate(:,373)*n(:,patmo_idx_34SO2)*n(:,patmo_idx_N2)
 
    endif
!
    dn(:,patmo_idx_34SO2) = &
        + krate(:,59)*n(:,patmo_idx_34SO)*n(:,patmo_idx_O3) &
        + krate(:,63)*n(:,patmo_idx_34SO)*n(:,patmo_idx_O2) &
        + krate(:,67)*n(:,patmo_idx_34SO)*n(:,patmo_idx_OH) &
        - krate(:,83)*n(:,patmo_idx_34SO2)*n(:,patmo_idx_HO2) &
        - krate(:,87)*n(:,patmo_idx_34SO2)*n(:,patmo_idx_O3) &
        + krate(:,91)*n(:,patmo_idx_34HSO)*n(:,patmo_idx_O2) &
        + krate(:,99)*n(:,patmo_idx_34HSO2)*n(:,patmo_idx_O2) &
        - krate(:,107)*n(:,patmo_idx_34SO2)*n(:,patmo_idx_O) &
        - krate(:,111)*n(:,patmo_idx_34SO2)*n(:,patmo_idx_OH) &
        + krate(:,119)*n(:,patmo_idx_34H2SO4) &
        - krate(:,125)*n(:,patmo_idx_34SO2) &
        + 0.99*(krate(:,129)*n(:,patmo_idx_34CH3SCH3)*n(:,patmo_idx_O)) &
        + 0.99*(krate(:,133)*n(:,patmo_idx_34CH3SCH3)*n(:,patmo_idx_OH)) &
        + 0.99*(krate(:,137)*n(:,patmo_idx_34CH3SCH3)*n(:,patmo_idx_OH)) &
        + krate(:,141)*n(:,patmo_idx_34SO2_1) &
        + krate(:,145)*n(:,patmo_idx_34SO2_1) &
        + krate(:,161)*n(:,patmo_idx_34SO2_3) &
        + krate(:,165)*n(:,patmo_idx_34SO2_3) &
        - krate(:,180)*n(:,patmo_idx_34SO2) &
        - krate(:,181)*n(:,patmo_idx_34SO2) &
        - krate(:,182)*n(:,patmo_idx_34SO2) &
        + krate(:,208)*n(:,patmo_idx_34SO3) &
        - krate(:,267)*n(:,patmo_idx_34SO2)*n(:,patmo_idx_O2) &
        - krate(:,271)*n(:,patmo_idx_34SO2)*n(:,patmo_idx_O) &
        - krate(:,275)*n(:,patmo_idx_34SO2)*n(:,patmo_idx_H) &
        + krate(:,291)*n(:,patmo_idx_OH)*n(:,patmo_idx_34SO3) &
        + krate(:,295)*n(:,patmo_idx_34SO3)*n(:,patmo_idx_O2) &
        - krate(:,299)*n(:,patmo_idx_34SO2)*n(:,patmo_idx_OH) &
        - krate(:,307)*n(:,patmo_idx_HO2)*n(:,patmo_idx_34SO2) &
        + krate(:,315)*n(:,patmo_idx_34SO3) &
        + krate(:,319)*n(:,patmo_idx_34HSO3) &
        - krate(:,327)*n(:,patmo_idx_34SO2)*n(:,patmo_idx_OH)*n(:,patmo_idx_OH) &
        + krate(:,333)*n(:,patmo_idx_34SO4) &
        - krate(:,337)*n(:,patmo_idx_34SO2) &
        - krate(:,341)*n(:,patmo_idx_34SO2) &
        - krate(:,345)*n(:,patmo_idx_34SO2)*n(:,patmo_idx_CH4O3S) &
        - krate(:,349)*n(:,patmo_idx_34SO2)*n(:,patmo_idx_N2) &
        - krate(:,353)*n(:,patmo_idx_34SO2) &
        - krate(:,369)*n(:,patmo_idx_34SO2) &
        - krate(:,373)*n(:,patmo_idx_34SO2)*n(:,patmo_idx_N2)

!Diffusion code
    ngas_hpp(:) = ngas_hp(:)/ngas_p(:)
    ngas_hpz(:) = ngas_hp(:)/ngas(:)
    ngas_hmm(:) = ngas_hm(:)/ngas_m(:)
    ngas_hmz(:) = ngas_hm(:)/ngas(:) 

    if (end_of_run) then 
        do j=1,60
	    budget(j,patmo_idx_32COS,3) = &
            + (k_hp(j)-d_hp(j,patmo_idx_32COS)) * ngas_hpp(j) * n_p(j,patmo_idx_32COS) &
            - (k_hp(j)+d_hp(j,patmo_idx_32COS)) * ngas_hpz(j) * n(j,patmo_idx_32COS)

            budget(j,patmo_idx_32SO2,3) = &
            + (k_hp(j)-d_hp(j,patmo_idx_32SO2)) * ngas_hpp(j) * n_p(j,patmo_idx_32SO2) &
            - (k_hp(j)+d_hp(j,patmo_idx_32SO2)) * ngas_hpz(j) * n(j,patmo_idx_32SO2)
!
!            budget(j,patmo_idx_33SO2,3) = &
!            + (k_hp(j)-d_hp(j,patmo_idx_33SO2)) * ngas_hpp(j) * n_p(j,patmo_idx_33SO2) &
!            - (k_hp(j)+d_hp(j,patmo_idx_33SO2)) * ngas_hpz(j) * n(j,patmo_idx_33SO2)

            budget(j,patmo_idx_34COS,3) = &
            + (k_hp(j)-d_hp(j,patmo_idx_34COS)) * ngas_hpp(j) * n_p(j,patmo_idx_34COS) &
            - (k_hp(j)+d_hp(j,patmo_idx_34COS)) * ngas_hpz(j) * n(j,patmo_idx_34COS)
!
            budget(j,patmo_idx_34SO2,3) = &
            + (k_hp(j)-d_hp(j,patmo_idx_34SO2)) * ngas_hpp(j) * n_p(j,patmo_idx_34SO2) &
            - (k_hp(j)+d_hp(j,patmo_idx_34SO2)) * ngas_hpz(j) * n(j,patmo_idx_34SO2)
!
!            budget(j,patmo_idx_36SO2,3) = &
!            + (k_hp(j)-d_hp(j,patmo_idx_36SO2)) * ngas_hpp(j) * n_p(j,patmo_idx_36SO2) &
!            - (k_hp(j)+d_hp(j,patmo_idx_36SO2)) * ngas_hpz(j) * n(j,patmo_idx_36SO2)
!
            budget(j,patmo_idx_32COS,4) = &
            + (k_hm(j)-d_hm(j,patmo_idx_32COS)) * ngas_hmm(j) * n_m(j,patmo_idx_32COS) &
            - (k_hm(j)+d_hm(j,patmo_idx_32COS)) * ngas_hmz(j) * n(j,patmo_idx_32COS)

            budget(j,patmo_idx_32SO2,4) = &
            + (k_hm(j)-d_hm(j,patmo_idx_32SO2)) * ngas_hmm(j) * n_m(j,patmo_idx_32SO2) &
            - (k_hm(j)+d_hm(j,patmo_idx_32SO2)) * ngas_hmz(j) * n(j,patmo_idx_32SO2)
!
!            budget(j,patmo_idx_33SO2,4) = &
!            + (k_hm(j)-d_hm(j,patmo_idx_33SO2)) * ngas_hmm(j) * n_m(j,patmo_idx_33SO2) &
!            - (k_hm(j)+d_hm(j,patmo_idx_33SO2)) * ngas_hmz(j) * n(j,patmo_idx_33SO2)
!
            budget(j,patmo_idx_34COS,4) = &
            + (k_hm(j)-d_hm(j,patmo_idx_34COS)) * ngas_hmm(j) * n_m(j,patmo_idx_34COS) &
            - (k_hm(j)+d_hm(j,patmo_idx_34COS)) * ngas_hmz(j) * n(j,patmo_idx_34COS)

            budget(j,patmo_idx_34SO2,4) = &
            + (k_hm(j)-d_hm(j,patmo_idx_34SO2)) * ngas_hmm(j) * n_m(j,patmo_idx_34SO2) &
            - (k_hm(j)+d_hm(j,patmo_idx_34SO2)) * ngas_hmz(j) * n(j,patmo_idx_34SO2)
            
             budget(j,patmo_idx_32CS2,3) = &
            + (k_hp(j)-d_hp(j,patmo_idx_32CS2)) * ngas_hpp(j) * n_p(j,patmo_idx_32CS2) &
            - (k_hp(j)+d_hp(j,patmo_idx_32CS2)) * ngas_hpz(j) * n(j,patmo_idx_32CS2)

            budget(j,patmo_idx_32H2S,3) = &
            + (k_hp(j)-d_hp(j,patmo_idx_32H2S)) * ngas_hpp(j) * n_p(j,patmo_idx_32H2S) &
            - (k_hp(j)+d_hp(j,patmo_idx_32H2S)) * ngas_hpz(j) * n(j,patmo_idx_32H2S)
!
!            budget(j,patmo_idx_33H2S,3) = &
!            + (k_hp(j)-d_hp(j,patmo_idx_33H2S)) * ngas_hpp(j) * n_p(j,patmo_idx_33H2S) &
!            - (k_hp(j)+d_hp(j,patmo_idx_33H2S)) * ngas_hpz(j) * n(j,patmo_idx_33H2S)

            budget(j,patmo_idx_34CS2,3) = &
            + (k_hp(j)-d_hp(j,patmo_idx_34CS2)) * ngas_hpp(j) * n_p(j,patmo_idx_34CS2) &
            - (k_hp(j)+d_hp(j,patmo_idx_34CS2)) * ngas_hpz(j) * n(j,patmo_idx_34CS2)
!
            budget(j,patmo_idx_34H2S,3) = &
            + (k_hp(j)-d_hp(j,patmo_idx_34H2S)) * ngas_hpp(j) * n_p(j,patmo_idx_34H2S) &
            - (k_hp(j)+d_hp(j,patmo_idx_34H2S)) * ngas_hpz(j) * n(j,patmo_idx_34H2S)

            budget(j,patmo_idx_32CH3SCH3,3) = &
            + (k_hp(j)-d_hp(j,patmo_idx_32CH3SCH3)) * ngas_hpp(j) * n_p(j,patmo_idx_32CH3SCH3) &
            - (k_hp(j)+d_hp(j,patmo_idx_32CH3SCH3)) * ngas_hpz(j) * n(j,patmo_idx_32CH3SCH3)
            
!
!            budget(j,patmo_idx_36H2S,3) = &
!            + (k_hp(j)-d_hp(j,patmo_idx_36H2S)) * ngas_hpp(j) * n_p(j,patmo_idx_36H2S) &
!            - (k_hp(j)+d_hp(j,patmo_idx_36H2S)) * ngas_hpz(j) * n(j,patmo_idx_36H2S)
!
            budget(j,patmo_idx_32CS2,4) = &
            + (k_hm(j)-d_hm(j,patmo_idx_32CS2)) * ngas_hmm(j) * n_m(j,patmo_idx_32CS2) &
            - (k_hm(j)+d_hm(j,patmo_idx_32CS2)) * ngas_hmz(j) * n(j,patmo_idx_32CS2)

            budget(j,patmo_idx_32H2S,4) = &
            + (k_hm(j)-d_hm(j,patmo_idx_32H2S)) * ngas_hmm(j) * n_m(j,patmo_idx_32H2S) &
            - (k_hm(j)+d_hm(j,patmo_idx_32H2S)) * ngas_hmz(j) * n(j,patmo_idx_32H2S)
!
!            budget(j,patmo_idx_33H2S,4) = &
!            + (k_hm(j)-d_hm(j,patmo_idx_33H2S)) * ngas_hmm(j) * n_m(j,patmo_idx_33H2S) &
!            - (k_hm(j)+d_hm(j,patmo_idx_33H2S)) * ngas_hmz(j) * n(j,patmo_idx_33H2S)
!
            budget(j,patmo_idx_34CS2,4) = &
            + (k_hm(j)-d_hm(j,patmo_idx_34CS2)) * ngas_hmm(j) * n_m(j,patmo_idx_34CS2) &
            - (k_hm(j)+d_hm(j,patmo_idx_34CS2)) * ngas_hmz(j) * n(j,patmo_idx_34CS2)

            budget(j,patmo_idx_34H2S,4) = &
            + (k_hm(j)-d_hm(j,patmo_idx_34H2S)) * ngas_hmm(j) * n_m(j,patmo_idx_34H2S) &
            - (k_hm(j)+d_hm(j,patmo_idx_34H2S)) * ngas_hmz(j) * n(j,patmo_idx_34H2S)

            budget(j,patmo_idx_32CH3SCH3,4) = &
            + (k_hm(j)-d_hm(j,patmo_idx_32CH3SCH3)) * ngas_hmm(j) * n_m(j,patmo_idx_32CH3SCH3) &
            - (k_hm(j)+d_hm(j,patmo_idx_32CH3SCH3)) * ngas_hmz(j) * n(j,patmo_idx_32CH3SCH3)
    
!
!            budget(j,patmo_idx_36SO2,4) = &
!            + (k_hm(j)-d_hm(j,patmo_idx_36SO2)) * ngas_hmm(j) * n_m(j,patmo_idx_36SO2) &
!            - (k_hm(j)+d_hm(j,patmo_idx_36SO2)) * ngas_hmz(j) * n(j,patmo_idx_36SO2)     
!            
            budget(j,patmo_idx_32H2SO4,3) = &
            + (k_hp(j)-d_hp(j,patmo_idx_32H2SO4)) * ngas_hpp(j) * n_p(j,patmo_idx_32H2SO4) &
            - (k_hp(j)+d_hp(j,patmo_idx_32H2SO4)) * ngas_hpz(j) * n(j,patmo_idx_32H2SO4)
!
            budget(j,patmo_idx_32SO4,3) = &
            + (k_hp(j)-d_hp(j,patmo_idx_32SO4)) * ngas_hpp(j) * n_p(j,patmo_idx_32SO4) &
            - (k_hp(j)+d_hp(j,patmo_idx_32SO4)) * ngas_hpz(j) * n(j,patmo_idx_32SO4)
!
            budget(j,patmo_idx_34H2SO4,3) = &
            + (k_hp(j)-d_hp(j,patmo_idx_34H2SO4)) * ngas_hpp(j) * n_p(j,patmo_idx_34H2SO4) &
            - (k_hp(j)+d_hp(j,patmo_idx_34H2SO4)) * ngas_hpz(j) * n(j,patmo_idx_34H2SO4)
!
            budget(j,patmo_idx_34SO4,3) = &
            + (k_hp(j)-d_hp(j,patmo_idx_34SO4)) * ngas_hpp(j) * n_p(j,patmo_idx_34SO4) &
            - (k_hp(j)+d_hp(j,patmo_idx_34SO4)) * ngas_hpz(j) * n(j,patmo_idx_34SO4)
!
            budget(j,patmo_idx_32H2SO4,4) = &
            + (k_hm(j)-d_hm(j,patmo_idx_32H2SO4)) * ngas_hmm(j) * n_m(j,patmo_idx_32H2SO4) &
            - (k_hm(j)+d_hm(j,patmo_idx_32H2SO4)) * ngas_hmz(j) * n(j,patmo_idx_32H2SO4)
!
            budget(j,patmo_idx_32SO4,4) = &
            + (k_hm(j)-d_hm(j,patmo_idx_32SO4)) * ngas_hmm(j) * n_m(j,patmo_idx_32SO4) &
            - (k_hm(j)+d_hm(j,patmo_idx_32SO4)) * ngas_hmz(j) * n(j,patmo_idx_32SO4)
!
            budget(j,patmo_idx_34H2SO4,4) = &
            + (k_hm(j)-d_hm(j,patmo_idx_34H2SO4)) * ngas_hmm(j) * n_m(j,patmo_idx_34H2SO4) &
            - (k_hm(j)+d_hm(j,patmo_idx_34H2SO4)) * ngas_hmz(j) * n(j,patmo_idx_34H2SO4)
!
            budget(j,patmo_idx_34SO4,4) = &
            + (k_hm(j)-d_hm(j,patmo_idx_34SO4)) * ngas_hmm(j) * n_m(j,patmo_idx_34SO4) &
            - (k_hm(j)+d_hm(j,patmo_idx_34SO4)) * ngas_hmz(j) * n(j,patmo_idx_34SO4)         
        enddo 
    endif
!
!
    do i=1,7
        do j=1,60
            dn(j,i) = dn(j,i) &
            + (k_hp(j)-d_hp(j,i)) * ngas_hpp(j) * n_p(j,i) &
            - ((k_hp(j)+d_hp(j,i)) * ngas_hpz(j) &
            + (k_hm(j)-d_hm(j,i)) * ngas_hmz(j)) * n(j,i) &
            + (k_hm(j)+d_hm(j,i)) * ngas_hmm(j) * n_m(j,i)
        end do
    end do

    do i=9,18
        do j=1,60
            dn(j,i) = dn(j,i) &
            + (k_hp(j)-d_hp(j,i)) * ngas_hpp(j) * n_p(j,i) &
            - ((k_hp(j)+d_hp(j,i)) * ngas_hpz(j) &
            + (k_hm(j)-d_hm(j,i)) * ngas_hmz(j)) * n(j,i) &
            + (k_hm(j)+d_hm(j,i)) * ngas_hmm(j) * n_m(j,i)
        end do
    end do

    do i=21,34
        do j=1,60
            dn(j,i) = dn(j,i) &
            + (k_hp(j)-d_hp(j,i)) * ngas_hpp(j) * n_p(j,i) &
            - ((k_hp(j)+d_hp(j,i)) * ngas_hpz(j) &
            + (k_hm(j)-d_hm(j,i)) * ngas_hmz(j)) * n(j,i) &
            + (k_hm(j)+d_hm(j,i)) * ngas_hmm(j) * n_m(j,i)
        end do
    end do

    do i=36,43
        do j=1,60
            dn(j,i) = dn(j,i) &
            + (k_hp(j)-d_hp(j,i)) * ngas_hpp(j) * n_p(j,i) &
            - ((k_hp(j)+d_hp(j,i)) * ngas_hpz(j) &
            + (k_hm(j)-d_hm(j,i)) * ngas_hmz(j)) * n(j,i) &
            + (k_hm(j)+d_hm(j,i)) * ngas_hmm(j) * n_m(j,i)
        end do
    end do

    do i=45,47
        do j=1,60
            dn(j,i) = dn(j,i) &
            + (k_hp(j)-d_hp(j,i)) * ngas_hpp(j) * n_p(j,i) &
            - ((k_hp(j)+d_hp(j,i)) * ngas_hpz(j) &
            + (k_hm(j)-d_hm(j,i)) * ngas_hmz(j)) * n(j,i) &
            + (k_hm(j)+d_hm(j,i)) * ngas_hmm(j) * n_m(j,i)
        end do
    end do
        
    do i=49,51
        do j=1,60
            dn(j,i) = dn(j,i) &
            + (k_hp(j)-d_hp(j,i)) * ngas_hpp(j) * n_p(j,i) &
            - ((k_hp(j)+d_hp(j,i)) * ngas_hpz(j) &
            + (k_hm(j)-d_hm(j,i)) * ngas_hmz(j)) * n(j,i) &
            + (k_hm(j)+d_hm(j,i)) * ngas_hmm(j) * n_m(j,i)
        end do
    end do

    do i=54,57
        do j=1,60
            dn(j,i) = dn(j,i) &
            + (k_hp(j)-d_hp(j,i)) * ngas_hpp(j) * n_p(j,i) &
            - ((k_hp(j)+d_hp(j,i)) * ngas_hpz(j) &
            + (k_hm(j)-d_hm(j,i)) * ngas_hmz(j)) * n(j,i) &
            + (k_hm(j)+d_hm(j,i)) * ngas_hmm(j) * n_m(j,i)
        end do
    end do

    do i=59,64
        do j=1,60
            dn(j,i) = dn(j,i) &
            + (k_hp(j)-d_hp(j,i)) * ngas_hpp(j) * n_p(j,i) &
            - ((k_hp(j)+d_hp(j,i)) * ngas_hpz(j) &
            + (k_hm(j)-d_hm(j,i)) * ngas_hmz(j)) * n(j,i) &
            + (k_hm(j)+d_hm(j,i)) * ngas_hmm(j) * n_m(j,i)
        end do
    end do

    do i=66,67
        do j=1,60
            dn(j,i) = dn(j,i) &
            + (k_hp(j)-d_hp(j,i)) * ngas_hpp(j) * n_p(j,i) &
            - ((k_hp(j)+d_hp(j,i)) * ngas_hpz(j) &
            + (k_hm(j)-d_hm(j,i)) * ngas_hmz(j)) * n(j,i) &
            + (k_hm(j)+d_hm(j,i)) * ngas_hmm(j) * n_m(j,i)
        end do
    end do

    do i=70,76
        do j=1,60
            dn(j,i) = dn(j,i) &
            + (k_hp(j)-d_hp(j,i)) * ngas_hpp(j) * n_p(j,i) &
            - ((k_hp(j)+d_hp(j,i)) * ngas_hpz(j) &
            + (k_hm(j)-d_hm(j,i)) * ngas_hmz(j)) * n(j,i) &
            + (k_hm(j)+d_hm(j,i)) * ngas_hmm(j) * n_m(j,i)
        end do
    end do

    do i=78,81
        do j=1,60
            dn(j,i) = dn(j,i) &
            + (k_hp(j)-d_hp(j,i)) * ngas_hpp(j) * n_p(j,i) &
            - ((k_hp(j)+d_hp(j,i)) * ngas_hpz(j) &
            + (k_hm(j)-d_hm(j,i)) * ngas_hmz(j)) * n(j,i) &
            + (k_hm(j)+d_hm(j,i)) * ngas_hmm(j) * n_m(j,i)
        end do
    end do
    
    !budget emission
    if(end_of_run) then 
        budget(1,patmo_idx_32COS,5) = + 6.8675359d7/1d5
        budget(1,patmo_idx_32SO2,5) = + 5.612610135d9/1d5
        budget(1,patmo_idx_34COS,5) = +  3.0747112d6/1d5
        budget(1,patmo_idx_34SO2,5) = + 2.499182710d8/1d5 
        budget(1,patmo_idx_32CS2,5) = + 3.743042d7/1d5
        budget(1,patmo_idx_34CS2,5) = + 1.675657d6/1d5
        budget(1,patmo_idx_32H2S,5) = + 1.046224903d9/1d5
        budget(1,patmo_idx_34H2S,5) = + 4.640087823d7/1d5
        budget(1,patmo_idx_32CH3SCH3,5) = + 2.507371315d9/1d5
        budget(1,patmo_idx_34CH3SCH3,5) = + 1.133146064d8/1d5
    endif
!    !budget emission end

    !emission
        !COS emissions total: 723 molec/cm3/s with 10.5 permil for 34S and NO MIF
        dn(1,patmo_idx_32COS) = dn(1,patmo_idx_32COS) + 6.8675359d7/1d5
        dn(1,patmo_idx_33COS) = dn(1,patmo_idx_33COS) + 5.4499096d5/1d5
        dn(1,patmo_idx_34COS) = dn(1,patmo_idx_34COS) + 3.0747112d6/1d5
        dn(1,patmo_idx_36COS) = dn(1,patmo_idx_36COS) + 1.4743303d4/1d5

        !CS2 emissions total: 394.1 molec/cm3/s with 10.4 permil
        dn(1,patmo_idx_32CS2) = dn(1,patmo_idx_32CS2) + 3.743042d7/1d5
        dn(1,patmo_idx_33CS2) = dn(1,patmo_idx_33CS2) + 2.970235d5/1d5
        dn(1,patmo_idx_34CS2) = dn(1,patmo_idx_34CS2) + 1.675657d6/1d5
        dn(1,patmo_idx_36CS2) = dn(1,patmo_idx_36CS2) + 8.034107d3/1d5

        !H2S emissions total: 11011 molec/cm3/s with 1 permil for 34S
        dn(1,patmo_idx_32H2S) = dn(1,patmo_idx_32H2S) + 1.046224903d9/1d5
        dn(1,patmo_idx_33H2S) = dn(1,patmo_idx_33H2S) + 8.258314469d6/1d5
        dn(1,patmo_idx_34H2S) = dn(1,patmo_idx_34H2S) + 4.640087823d7/1d5
        dn(1,patmo_idx_36H2S) = dn(1,patmo_idx_36H2S) + 2.206299159d5/1d5

        !SO2 emissions total: 59081 mole/cm3/s with 5 permil for 34S
        dn(1,patmo_idx_32SO2) = dn(1,patmo_idx_32SO2) + 5.612610135d9/1d5
        dn(1,patmo_idx_33SO2) = dn(1,patmo_idx_33SO2) + 4.431101857d7/1d5
        dn(1,patmo_idx_34SO2) = dn(1,patmo_idx_34SO2) + 2.499182710d8/1d5
        dn(1,patmo_idx_36SO2) = dn(1,patmo_idx_36SO2) + 1.192576285d6/1d5

        !DMS emissions total: 26406 molec/cm3/s with 20 permil for 34S
        dn(1,patmo_idx_32CH3SCH3) = dn(1,patmo_idx_32CH3SCH3) + 2.507371315d9/1d5
        dn(1,patmo_idx_33CH3SCH3) = dn(1,patmo_idx_33CH3SCH3) + 1.999471616d7/1d5
        dn(1,patmo_idx_34CH3SCH3) = dn(1,patmo_idx_34CH3SCH3) + 1.133146064d8/1d5
        dn(1,patmo_idx_36CH3SCH3) = dn(1,patmo_idx_36CH3SCH3) + 5.478112870d5/1d5
    !end emissions

    !budget dry_dep (lifetime tuned by 5.55)
    if(end_of_run) then 
         budget(1,patmo_idx_32SO2,6) = - (5.55*1.355d-7)*n(1,patmo_idx_32SO2)
!        budget(1,patmo_idx_33SO2,6) = -(5.5*1.355d-7)*n(1,patmo_idx_33SO2)
         budget(1,patmo_idx_34SO2,6) = - (5.55*1.355d-7)*n(1,patmo_idx_34SO2)
!        budget(1,patmo_idx_36SO2,6) = - (5.5*1.355d-7)*n(1,patmo_idx_36SO2)
	 budget(1,patmo_idx_32COS,6) = - (5.55*1.585d-8)*n(1,patmo_idx_32COS)
	 budget(1,patmo_idx_34COS,6) = -(5.55*1.5819885d-8)*n(1,patmo_idx_34COS)
         budget(1,patmo_idx_32CS2,6) = - ((5.4d-9*5.55))*n(1,patmo_idx_32CS2)
         budget(1,patmo_idx_34CS2,6) = - ((5.4d-9*5.55))*n(1,patmo_idx_34CS2)

    endif
!    !budget dry dep end
!
  !dry deposition (lifetime tuned by 5.55 to offset high emissions)
        !COS dry deposition (tuned to offset higher emissions + with -1.9 permil fractionation for biogenic uptake for 34S as per Angert et al., 2020)
        dn(1,patmo_idx_32COS) = dn(1,patmo_idx_32COS) - (5.55*1.585d-8)*n(1,patmo_idx_32COS)
        dn(1,patmo_idx_33COS) = dn(1,patmo_idx_33COS) -(5.55*1.585d-8)*n(1,patmo_idx_33COS)
        dn(1,patmo_idx_34COS) = dn(1,patmo_idx_34COS) -(5.55*1.5819885d-8)*n(1,patmo_idx_34COS)
        dn(1,patmo_idx_36COS) = dn(1,patmo_idx_36COS) -(5.55*1.585d-8)*n(1,patmo_idx_36COS)

        !SO2 dry deposition (also tuned), no fractionation
        dn(1,patmo_idx_32SO2) = dn(1,patmo_idx_32SO2) - (5.55*1.355d-7)*n(1,patmo_idx_32SO2)
        dn(1,patmo_idx_33SO2) = dn(1,patmo_idx_33SO2) -(5.55*1.355d-7)*n(1,patmo_idx_33SO2)
        dn(1,patmo_idx_34SO2) = dn(1,patmo_idx_34SO2) - (5.55*1.355d-7)*n(1,patmo_idx_34SO2)
        dn(1,patmo_idx_36SO2) = dn(1,patmo_idx_36SO2) - (5.55*1.355d-7)*n(1,patmo_idx_36SO2)
    
       ! CS2 dry deposition (also tuned), no fractionation
        dn(1,patmo_idx_32CS2) = dn(1,patmo_idx_32CS2) - ((5.4d-9*5.55))*n(1,patmo_idx_32CS2)
        dn(1,patmo_idx_33CS2) = dn(1,patmo_idx_33CS2) - ((5.4d-9*5.55))*n(1,patmo_idx_33CS2)
        dn(1,patmo_idx_34CS2) = dn(1,patmo_idx_34CS2) - ((5.4d-9*5.55))*n(1,patmo_idx_34CS2)
        dn(1,patmo_idx_36CS2) = dn(1,patmo_idx_36CS2) - ((5.4d-9*5.55))*n(1,patmo_idx_36CS2)
    !dry deposition end

!    !aerosol formation budget stored in emission for SO4/H2SO4
    if (end_of_run) then
        do i=13,34
            if (va32(i) <= n(i,15) .and. pa32(i) >= n(i,15)) then
                budget(i,patmo_idx_32H2SO4,5) = - (n(i,15)-va32(i))
                budget(i,patmo_idx_32SO4,5) = (n(i,15)-va32(i))
            endif
            if (va34(i) <= n(i,36) .and. pa34(i) >= n(i,36)) then
                budget(i,patmo_idx_34H2SO4,5) = - (n(i,36)-va34(i))
                budget(i,patmo_idx_34SO4,5) = (n(i,36)-va34(i))
            endif
        enddo
    endif

    !aerosol formation
    do i=13,34
        if (va32(i) <= n(i,15) .and. pa32(i) >= n(i,15)) then
            dn(i,patmo_idx_32H2SO4) = dn(i,patmo_idx_32H2SO4) - (n(i,15)-va32(i))
            dn(i,patmo_idx_32SO4) = dn(i,patmo_idx_32SO4) + (n(i,15)-va32(i))
        end if
    end do
  
    do i=13,34
        if (va33(i) <= n(i,1) .and. pa33(i) >= n(i,1)) then
            dn(i,patmo_idx_33H2SO4) = dn(i,patmo_idx_33H2SO4) - (n(i,1)-va33(i))
            dn(i,patmo_idx_33SO4) = dn(i,patmo_idx_33SO4) + (n(i,1)-va33(i))
        end if
    end do
   
    do i=13,34
        if (va34(i) <= n(i,36) .and. pa34(i) >= n(i,36)) then
            dn(i,patmo_idx_34H2SO4) = dn(i,patmo_idx_34H2SO4) - (n(i,36)-va34(i))
            dn(i,patmo_idx_34SO4) = dn(i,patmo_idx_34SO4) + (n(i,36)-va34(i))
        end if
    end do
   
    do i=13,34
        if (va36(i) <= n(i,17) .and. pa36(i) >= n(i,17)) then
            dn(i,patmo_idx_36H2SO4) = dn(i,patmo_idx_36H2SO4) - (n(i,17)-va36(i))
            dn(i,patmo_idx_36SO4) = dn(i,patmo_idx_36SO4) + (n(i,17)-va36(i))
        end if
    end do
    !aerosol formation

!    ! budget for gravity settling for SO4 stored in dry deposition marker
    if (end_of_run) then
        budget(60,patmo_idx_32SO4,6) = - gd(60)*n(60,patmo_idx_32SO4)
!        budget(60,patmo_idx_33SO4,6) = - gd(60)*n(60,patmo_idx_33SO4)
        budget(60,patmo_idx_34SO4,6) = - gd(60)*n(60,patmo_idx_34SO4)
!        budget(60,patmo_idx_36SO4,6) = - gd(60)*n(60,patmo_idx_36SO4)
        do j=59,1,-1
            budget(j,patmo_idx_32SO4,6) = gd(j+1)*n(j+1,patmo_idx_32SO4) & 
                - gd(j)*n(j,patmo_idx_32SO4)
!            budget(j,patmo_idx_33SO4,6) = gd(j+1)*n(j+1,patmo_idx_33SO4) &
!                - gd(j)*n(j,patmo_idx_33SO4)
            budget(j,patmo_idx_34SO4,6) = gd(j+1)*n(j+1,patmo_idx_34SO4) & 
                - gd(j)*n(j,patmo_idx_34SO4)
!            budget(j,patmo_idx_36SO4,6) = gd(j+1)*n(j+1,patmo_idx_36SO4) &
!                - gd(j)*n(j,patmo_idx_36SO4)
        enddo 
    endif 
!    !gravity settling SO4 Aerosol (JAM-Kasten-1968,r=0.3)
    do j=60,2,-1
        dn(j,patmo_idx_32SO4) = dn(j,patmo_idx_32SO4)-gd(j)*n(j,patmo_idx_32SO4)
        dn(j-1,patmo_idx_32SO4) = dn(j-1,patmo_idx_32SO4)+gd(j)*n(j,patmo_idx_32SO4)
    end do  
    dn(1,patmo_idx_32SO4) = dn(1,patmo_idx_32SO4)-gd(1)*n(1,patmo_idx_32SO4)

    do j=60,2,-1
        dn(j,patmo_idx_33SO4) = dn(j,patmo_idx_33SO4)-gd(j)*n(j,patmo_idx_33SO4)
        dn(j-1,patmo_idx_33SO4) = dn(j-1,patmo_idx_33SO4)+gd(j)*n(j,patmo_idx_33SO4)
    end do  
    dn(1,patmo_idx_33SO4) = dn(1,patmo_idx_33SO4)-gd(1)*n(1,patmo_idx_33SO4)

    do j=60,2,-1
        dn(j,patmo_idx_34SO4) = dn(j,patmo_idx_34SO4)-gd(j)*n(j,patmo_idx_34SO4)
        dn(j-1,patmo_idx_34SO4) = dn(j-1,patmo_idx_34SO4)+gd(j)*n(j,patmo_idx_34SO4)
    end do  
    dn(1,patmo_idx_34SO4) = dn(1,patmo_idx_34SO4)-gd(1)*n(1,patmo_idx_34SO4)
   
    do j=60,2,-1
        dn(j,patmo_idx_36SO4) = dn(j,patmo_idx_36SO4)-gd(j)*n(j,patmo_idx_36SO4)
        dn(j-1,patmo_idx_36SO4) = dn(j-1,patmo_idx_36SO4)+gd(j)*n(j,patmo_idx_36SO4)
    end do  
    dn(1,patmo_idx_36SO4) = dn(1,patmo_idx_36SO4)-gd(1)*n(1,patmo_idx_36SO4)
    !gravity settling

!    !budget wet_dep
     if (end_of_run) then
        budget(12,patmo_idx_32COS,7) = - wetdep(12,patmo_idx_32COS)*n(12,patmo_idx_32COS) 
        budget(12,patmo_idx_33COS,7) = - wetdep(12,patmo_idx_33COS)*n(12,patmo_idx_33COS)
        budget(12,patmo_idx_34COS,7) = - wetdep(12,patmo_idx_34COS)*n(12,patmo_idx_34COS)
        budget(12,patmo_idx_36COS,7) = - wetdep(12,patmo_idx_36COS)*n(12,patmo_idx_36COS)
	do j=11,1,-1
            budget(j,patmo_idx_32COS,7) = wetdep(j+1,patmo_idx_32COS)*n(j+1,patmo_idx_32COS) & 
                - wetdep(j,patmo_idx_32COS)*n(j,patmo_idx_32COS) 
            budget(j,patmo_idx_33COS,7) = wetdep(j+1,patmo_idx_33COS)*n(j+1,patmo_idx_33COS) &
                - wetdep(j,patmo_idx_33COS)*n(j,patmo_idx_33COS)
            budget(j,patmo_idx_34COS,7) = wetdep(j+1,patmo_idx_34COS)*n(j+1,patmo_idx_34COS) & 
                - wetdep(j,patmo_idx_34COS)*n(j,patmo_idx_34COS)
            budget(j,patmo_idx_36COS,7) = wetdep(j+1,patmo_idx_36COS)*n(j+1,patmo_idx_36COS) &
                - wetdep(j,patmo_idx_36COS)*n(j,patmo_idx_36COS)
	enddo
     endif
    if (end_of_run) then
        budget(12,patmo_idx_32SO2,39) = - wetdep(12,patmo_idx_32SO2)*n(12,patmo_idx_32SO2) 
!        budget(12,patmo_idx_33SO2,7) = - wetdep(12,patmo_idx_33SO2)*n(12,patmo_idx_33SO2)
        budget(12,patmo_idx_34SO2,39) = - wetdep(12,patmo_idx_34SO2)*n(12,patmo_idx_34SO2)
!        budget(12,patmo_idx_36SO2,7) = - wetdep(12,patmo_idx_36SO2)*n(12,patmo_idx_36SO2)
!        
        budget(12,patmo_idx_32H2SO4,7) = - wetdep(12,patmo_idx_32H2SO4)*n(12,patmo_idx_32H2SO4) 
        budget(12,patmo_idx_32SO4,7) = - wetdep(12,patmo_idx_32SO4)*n(12,patmo_idx_32SO4)
        budget(12,patmo_idx_34H2SO4,7) = - wetdep(12,patmo_idx_34H2SO4)*n(12,patmo_idx_34H2SO4)
        budget(12,patmo_idx_34SO4,7) = - wetdep(12,patmo_idx_34SO4)*n(12,patmo_idx_34SO4)
!    
        do j=11,1,-1
            budget(j,patmo_idx_32SO2,39) = wetdep(j+1,patmo_idx_32SO2)*n(j+1,patmo_idx_32SO2) & 
                - wetdep(j,patmo_idx_32SO2)*n(j,patmo_idx_32SO2) 
!            budget(j,patmo_idx_33SO2,7) = wetdep(j+1,patmo_idx_33SO2)*n(j+1,patmo_idx_33SO2) &
!                - wetdep(j,patmo_idx_33SO2)*n(j,patmo_idx_33SO2)
            budget(j,patmo_idx_34SO2,39) = wetdep(j+1,patmo_idx_34SO2)*n(j+1,patmo_idx_34SO2) & 
                - wetdep(j,patmo_idx_34SO2)*n(j,patmo_idx_34SO2)
!            budget(j,patmo_idx_36SO2,7) = wetdep(j+1,patmo_idx_36SO2)*n(j+1,patmo_idx_36SO2) &
!                - wetdep(j,patmo_idx_36SO2)*n(j,patmo_idx_36SO2)
!            
            budget(j,patmo_idx_32H2SO4,7) = wetdep(j+1,patmo_idx_32H2SO4)*n(j+1,patmo_idx_32H2SO4) & 
                - wetdep(j,patmo_idx_32H2SO4)*n(j,patmo_idx_32H2SO4) 
            budget(j,patmo_idx_32SO4,7) = wetdep(j+1,patmo_idx_32SO4)*n(j+1,patmo_idx_32SO4) &
                - wetdep(j,patmo_idx_32SO4)*n(j,patmo_idx_32SO4)
            budget(j,patmo_idx_34H2SO4,7) = wetdep(j+1,patmo_idx_34H2SO4)*n(j+1,patmo_idx_34H2SO4) & 
                - wetdep(j,patmo_idx_34H2SO4)*n(j,patmo_idx_34H2SO4)
            budget(j,patmo_idx_34SO4,7) = wetdep(j+1,patmo_idx_34SO4)*n(j+1,patmo_idx_34SO4) &
                - wetdep(j,patmo_idx_34SO4)*n(j,patmo_idx_34SO4)
        enddo 
    endif 
!     
    !end budget_wetdep
 
    !wet deposition
    do j=12,2,-1
        do i=1,chemSpeciesNumber
            dn(j,i) = dn(j,i)-wetdep(j,i)*n(j,i)
            dn(j-1,i) = dn(j-1,i)+wetdep(j,i)*n(j,i)
        end do   
    end do
    
    do i=1,chemSpeciesNumber
        dn(1,i) = dn(1,i)-wetdep(1,i)*n(1,i)
    end do  

    
    !unroll chemistry
    dy(:) = 0d0
    do i=1,speciesNumber
        dy((i-1)*cellsNumber+1:(i*cellsNumber)) = dn(:,i)
    end do

  end subroutine fex
end module patmo_ode
