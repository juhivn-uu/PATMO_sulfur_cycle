module patmo_ode
contains
  subroutine fex(neq,tt,nin,dy)
    use patmo_commons
    use patmo_constants
    use patmo_parameters
    use patmo_utils
    use patmo_budget, only: budget!budget module added JV
    implicit none
    integer,intent(in)::neq
    real*8,intent(in)::tt,nin(neqAll)
    real*8,intent(out)::dy(neqAll)
    INTEGER iout
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
    integer::i,j
    real*8::contribution(cellsNumber,speciesNumber)

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

    k_hp(:) = 100*(kzz_hp(:)+dzz_hp(:))*idh2(:)
    k_hm(:) = 100*(kzz_hm(:)+dzz_hm(:))*idh2(:)

! write(*,*) 
! STOP

    dn(:,:) = 0d0
    if (end_of_run) then
      budget(:,patmo_idx_COS,1) = &
        + 0.83*krate(:,5)*n(:,patmo_idx_CS)*n(:,patmo_idx_O2) &
        + krate(:,6)*n(:,patmo_idx_CS)*n(:,patmo_idx_O3) &
        + krate(:,46)*n(:,patmo_idx_CO2)*n(:,patmo_idx_SH) &
        + krate(:,47)*n(:,patmo_idx_CO)*n(:,patmo_idx_SO) &
        + 0.007*(krate(:,35)*n(:,patmo_idx_CH3SCH3)*n(:,patmo_idx_O)) &
        + 0.007*(krate(:,36)*n(:,patmo_idx_CH3SCH3)*n(:,patmo_idx_OH)) &        
        + 0.007*(krate(:,37)*n(:,patmo_idx_CH3SCH3)*n(:,patmo_idx_OH)) 



      budget(:,patmo_idx_COS,2) = &
        + krate(:,1)*n(:,patmo_idx_COS)*n(:,patmo_idx_OH) &
        + krate(:,2)*n(:,patmo_idx_COS)*n(:,patmo_idx_O) &
        + krate(:,38)*n(:,patmo_idx_COS) &
        + krate(:,48)*n(:,patmo_idx_SH)*n(:,patmo_idx_COS) &
        + krate(:,50)*n(:,patmo_idx_COS)*n(:,patmo_idx_O) &
        + krate(:,51)*n(:,patmo_idx_COS)*n(:,patmo_idx_O2)
      
      budget(:,patmo_idx_COS,8) = &
         + krate(:,38)*n(:,patmo_idx_COS)
      endif
              
    dn(:,patmo_idx_COS) = &
        - krate(:,1)*n(:,patmo_idx_COS)*n(:,patmo_idx_OH) &
        - krate(:,2)*n(:,patmo_idx_COS)*n(:,patmo_idx_O) &
        + 0.83*(krate(:,3)*n(:,patmo_idx_CS2)*n(:,patmo_idx_OH)) &
        + krate(:,5)*n(:,patmo_idx_CS)*n(:,patmo_idx_O2) &
        + krate(:,6)*n(:,patmo_idx_CS)*n(:,patmo_idx_O3) &
        - krate(:,38)*n(:,patmo_idx_COS) &
        + krate(:,46)*n(:,patmo_idx_CO2)*n(:,patmo_idx_SH) &
        + krate(:,47)*n(:,patmo_idx_CO)*n(:,patmo_idx_SO) &
        - krate(:,48)*n(:,patmo_idx_SH)*n(:,patmo_idx_COS) &
        - krate(:,50)*n(:,patmo_idx_COS)*n(:,patmo_idx_O) &
        - krate(:,51)*n(:,patmo_idx_COS)*n(:,patmo_idx_O2) 
        + 0.007*(krate(:,35)*n(:,patmo_idx_CH3SCH3)*n(:,patmo_idx_O)) &        
        + 0.007*((krate(:,36)*n(:,patmo_idx_CH3SCH3)*n(:,patmo_idx_OH)) &
        + 0.007*((krate(:,37)*n(:,patmo_idx_CH3SCH3)*n(:,patmo_idx_OH)) 

    dn(:,patmo_idx_HO2) = 0d0
    !    - krate(:,11)*n(:,patmo_idx_H2S)*n(:,patmo_idx_HO2) &
    !    - krate(:,22)*n(:,patmo_idx_SO2)*n(:,patmo_idx_HO2) &
    !    + krate(:,26)*n(:,patmo_idx_HSO2)*n(:,patmo_idx_O2) &
    !    + krate(:,27)*n(:,patmo_idx_HSO3)*n(:,patmo_idx_O2) &
    !    + krate(:,56)*n(:,patmo_idx_H2O)*n(:,patmo_idx_HSO) &
    !    + krate(:,67)*n(:,patmo_idx_OH)*n(:,patmo_idx_SO3) &
    !    - krate(:,71)*n(:,patmo_idx_HO2)*n(:,patmo_idx_SO2) &
    !    - krate(:,72)*n(:,patmo_idx_HO2)*n(:,patmo_idx_SO3)

    dn(:,patmo_idx_NO) = 0d0
    !    + krate(:,18)*n(:,patmo_idx_SO)*n(:,patmo_idx_NO2) &
    !    - krate(:,63)*n(:,patmo_idx_SO2)*n(:,patmo_idx_NO)

    dn(:,patmo_idx_N) = 0d0
    !    + krate(:,33)*n(:,patmo_idx_N2) &
    !    + krate(:,33)*n(:,patmo_idx_N2) &
    !    - krate(:,78)*n(:,patmo_idx_N)*n(:,patmo_idx_N) &
    !    - krate(:,78)*n(:,patmo_idx_N)*n(:,patmo_idx_N)

    dn(:,patmo_idx_HSO) = &
        + krate(:,11)*n(:,patmo_idx_H2S)*n(:,patmo_idx_HO2) &
        + krate(:,14)*n(:,patmo_idx_SH)*n(:,patmo_idx_O3) &
        - krate(:,24)*n(:,patmo_idx_HSO)*n(:,patmo_idx_O2) &
        - krate(:,25)*n(:,patmo_idx_HSO)*n(:,patmo_idx_O3) &
        - krate(:,56)*n(:,patmo_idx_H2O)*n(:,patmo_idx_HSO) &
        - krate(:,59)*n(:,patmo_idx_HSO)*n(:,patmo_idx_O2) &
        + krate(:,69)*n(:,patmo_idx_SO2)*n(:,patmo_idx_OH) &
        + krate(:,70)*n(:,patmo_idx_O2)*n(:,patmo_idx_O2)*n(:,patmo_idx_SH)

    dn(:,patmo_idx_CO2) = 0d0
     !   + krate(:,1)*n(:,patmo_idx_COS)*n(:,patmo_idx_OH) &
     !   - krate(:,46)*n(:,patmo_idx_CO2)*n(:,patmo_idx_SH)

    dn(:,patmo_idx_SO3) = &
        + krate(:,22)*n(:,patmo_idx_SO2)*n(:,patmo_idx_HO2) &
        + krate(:,23)*n(:,patmo_idx_SO2)*n(:,patmo_idx_O3) &
        + krate(:,27)*n(:,patmo_idx_HSO3)*n(:,patmo_idx_O2) &
        + krate(:,28)*n(:,patmo_idx_SO2)*n(:,patmo_idx_O) &
        - krate(:,30)*n(:,patmo_idx_SO3)*n(:,patmo_idx_H2O) &
        - krate(:,42)*n(:,patmo_idx_SO3) &
        - krate(:,67)*n(:,patmo_idx_OH)*n(:,patmo_idx_SO3) &
        - krate(:,68)*n(:,patmo_idx_SO3)*n(:,patmo_idx_O2) &
        - krate(:,72)*n(:,patmo_idx_HO2)*n(:,patmo_idx_SO3) &
        - krate(:,73)*n(:,patmo_idx_SO3) &
        + krate(:,75)*n(:,patmo_idx_H2SO4)

    dn(:,patmo_idx_H2O) = 0d0
     !   + krate(:,8)*n(:,patmo_idx_H2S)*n(:,patmo_idx_OH) &
     !   + krate(:,11)*n(:,patmo_idx_H2S)*n(:,patmo_idx_HO2) &
     !   - krate(:,30)*n(:,patmo_idx_SO3)*n(:,patmo_idx_H2O) &
     !   - krate(:,53)*n(:,patmo_idx_H2O)*n(:,patmo_idx_SH) &
     !   - krate(:,56)*n(:,patmo_idx_H2O)*n(:,patmo_idx_HSO) &
     !   + krate(:,75)*n(:,patmo_idx_H2SO4)

    dn(:,patmo_idx_HSO2) = &
        - krate(:,26)*n(:,patmo_idx_HSO2)*n(:,patmo_idx_O2) &
        + krate(:,71)*n(:,patmo_idx_HO2)*n(:,patmo_idx_SO2)
    if (end_of_run) then
            budget(:,patmo_idx_HSO2,1) = krate(:,71)*n(:,patmo_idx_HO2)*n(:,patmo_idx_SO2)
            budget(:,patmo_idx_HSO2,2) = krate(:,26)*n(:,patmo_idx_HSO2)*n(:,patmo_idx_O2)
    endif
    dn(:,patmo_idx_CO) = 0d0
     !   + krate(:,2)*n(:,patmo_idx_COS)*n(:,patmo_idx_O) &
     !   + krate(:,7)*n(:,patmo_idx_CS)*n(:,patmo_idx_O) &
     !   + krate(:,38)*n(:,patmo_idx_COS) &
     !   - krate(:,47)*n(:,patmo_idx_CO)*n(:,patmo_idx_SO) &
     !   - krate(:,52)*n(:,patmo_idx_CO)*n(:,patmo_idx_S)

    dn(:,patmo_idx_O2) = 0d0
     !   - krate(:,5)*n(:,patmo_idx_CS)*n(:,patmo_idx_O2) &
     !   + krate(:,6)*n(:,patmo_idx_CS)*n(:,patmo_idx_O3) &
     !   - krate(:,13)*n(:,patmo_idx_SH)*n(:,patmo_idx_O2) &
     !   + krate(:,14)*n(:,patmo_idx_SH)*n(:,patmo_idx_O3) &
     !   + krate(:,15)*n(:,patmo_idx_SO)*n(:,patmo_idx_O3) &
     !   - krate(:,16)*n(:,patmo_idx_SO)*n(:,patmo_idx_O2) &
     !   - krate(:,19)*n(:,patmo_idx_S)*n(:,patmo_idx_O2) &
     !   + krate(:,20)*n(:,patmo_idx_S)*n(:,patmo_idx_O3) &
     !   + krate(:,23)*n(:,patmo_idx_SO2)*n(:,patmo_idx_O3) &
     !   - krate(:,24)*n(:,patmo_idx_HSO)*n(:,patmo_idx_O2) &
     !   + krate(:,25)*n(:,patmo_idx_HSO)*n(:,patmo_idx_O3) &
     !   + krate(:,25)*n(:,patmo_idx_HSO)*n(:,patmo_idx_O3) &
     !   - krate(:,26)*n(:,patmo_idx_HSO2)*n(:,patmo_idx_O2) &
     !   - krate(:,27)*n(:,patmo_idx_HSO3)*n(:,patmo_idx_O2) &
     !   - krate(:,32)*n(:,patmo_idx_O)*n(:,patmo_idx_O2) &
     !   + krate(:,39)*n(:,patmo_idx_O3) &
     !   - krate(:,40)*n(:,patmo_idx_O2) &
     !   + krate(:,50)*n(:,patmo_idx_COS)*n(:,patmo_idx_O) &
     !   - krate(:,51)*n(:,patmo_idx_COS)*n(:,patmo_idx_O2) &
     !   + krate(:,58)*n(:,patmo_idx_OH)*n(:,patmo_idx_SO) &
     !   - krate(:,59)*n(:,patmo_idx_HSO)*n(:,patmo_idx_O2) &
     !   - krate(:,60)*n(:,patmo_idx_SO2)*n(:,patmo_idx_O2) &
     !   + krate(:,61)*n(:,patmo_idx_SO2)*n(:,patmo_idx_O) &
     !   + krate(:,64)*n(:,patmo_idx_SO)*n(:,patmo_idx_O) &
     !   - krate(:,65)*n(:,patmo_idx_O2)*n(:,patmo_idx_SO) &
     !   - krate(:,68)*n(:,patmo_idx_SO3)*n(:,patmo_idx_O2) &
     !   + krate(:,69)*n(:,patmo_idx_SO2)*n(:,patmo_idx_OH) &
     !   - krate(:,70)*n(:,patmo_idx_O2)*n(:,patmo_idx_O2)*n(:,patmo_idx_SH) &
     !   - krate(:,70)*n(:,patmo_idx_O2)*n(:,patmo_idx_O2)*n(:,patmo_idx_SH) &
     !   + krate(:,71)*n(:,patmo_idx_HO2)*n(:,patmo_idx_SO2) &
     !   + krate(:,72)*n(:,patmo_idx_HO2)*n(:,patmo_idx_SO3) &
     !   + krate(:,77)*n(:,patmo_idx_O3)

    dn(:,patmo_idx_N2) = 0d0
     !   - krate(:,33)*n(:,patmo_idx_N2) &
     !   + krate(:,78)*n(:,patmo_idx_N)*n(:,patmo_idx_N)
    
    if (end_of_run) then
      budget(:,patmo_idx_CS2,1) = &
        + krate(:,48)*n(:,patmo_idx_SH)*n(:,patmo_idx_COS) &
        + krate(:,49)*n(:,patmo_idx_CS)*n(:,patmo_idx_SO)

    budget(:,patmo_idx_CS2,2) = &
        + krate(:,3)*n(:,patmo_idx_CS2)*n(:,patmo_idx_OH) &
        + krate(:,4)*n(:,patmo_idx_CS2)*n(:,patmo_idx_O) &
        + krate(:,41)*n(:,patmo_idx_CS2)
    budget(:,patmo_idx_CS2,8) = &
        + krate(:,41)*n(:,patmo_idx_CS2)
        
    endif

    dn(:,patmo_idx_CS2) = &
        - krate(:,3)*n(:,patmo_idx_CS2)*n(:,patmo_idx_OH) &
        - krate(:,4)*n(:,patmo_idx_CS2)*n(:,patmo_idx_O) &
        - krate(:,41)*n(:,patmo_idx_CS2) &
        + krate(:,48)*n(:,patmo_idx_SH)*n(:,patmo_idx_COS) &
        + krate(:,49)*n(:,patmo_idx_CS)*n(:,patmo_idx_SO)

    dn(:,patmo_idx_SO) = &
        + krate(:,2)*n(:,patmo_idx_COS)*n(:,patmo_idx_O) &
        + krate(:,4)*n(:,patmo_idx_CS2)*n(:,patmo_idx_O) &
        + krate(:,12)*n(:,patmo_idx_SH)*n(:,patmo_idx_O) &
        + krate(:,13)*n(:,patmo_idx_SH)*n(:,patmo_idx_O2) &
        - krate(:,15)*n(:,patmo_idx_SO)*n(:,patmo_idx_O3) &
        - krate(:,16)*n(:,patmo_idx_SO)*n(:,patmo_idx_O2) &
        - krate(:,17)*n(:,patmo_idx_SO)*n(:,patmo_idx_OH) &
        - krate(:,18)*n(:,patmo_idx_SO)*n(:,patmo_idx_NO2) &
        + krate(:,19)*n(:,patmo_idx_S)*n(:,patmo_idx_O2) &
        + krate(:,20)*n(:,patmo_idx_S)*n(:,patmo_idx_O3) &
        + krate(:,21)*n(:,patmo_idx_S)*n(:,patmo_idx_OH) &
        + krate(:,43)*n(:,patmo_idx_SO2) &
        - krate(:,45)*n(:,patmo_idx_SO) &
        - krate(:,47)*n(:,patmo_idx_CO)*n(:,patmo_idx_SO) &
        - krate(:,49)*n(:,patmo_idx_CS)*n(:,patmo_idx_SO) &
        - krate(:,57)*n(:,patmo_idx_H)*n(:,patmo_idx_SO) &
        - krate(:,58)*n(:,patmo_idx_OH)*n(:,patmo_idx_SO) &
        + krate(:,60)*n(:,patmo_idx_SO2)*n(:,patmo_idx_O2) &
        + krate(:,61)*n(:,patmo_idx_SO2)*n(:,patmo_idx_O) &
        + krate(:,62)*n(:,patmo_idx_SO2)*n(:,patmo_idx_H) &
        + krate(:,63)*n(:,patmo_idx_SO2)*n(:,patmo_idx_NO) &
        - krate(:,64)*n(:,patmo_idx_SO)*n(:,patmo_idx_O) &
        - krate(:,65)*n(:,patmo_idx_O2)*n(:,patmo_idx_SO) &
        - krate(:,66)*n(:,patmo_idx_H)*n(:,patmo_idx_SO)
    if (end_of_run) then
            budget(:,patmo_idx_SO,1) =  + krate(:,2)*n(:,patmo_idx_COS)*n(:,patmo_idx_O) 
            budget(:,patmo_idx_SO,2) =  + krate(:,43)*n(:,patmo_idx_SO2) 
            budget(:,patmo_idx_SO,3) =  + krate(:,60)*n(:,patmo_idx_SO2)*n(:,patmo_idx_O2) &
        + krate(:,61)*n(:,patmo_idx_SO2)*n(:,patmo_idx_O) &
        + krate(:,62)*n(:,patmo_idx_SO2)*n(:,patmo_idx_H) &
        + krate(:,63)*n(:,patmo_idx_SO2)*n(:,patmo_idx_NO) 
            budget(:,patmo_idx_SO,4) =   - krate(:,15)*n(:,patmo_idx_SO)*n(:,patmo_idx_O3) &
        - krate(:,16)*n(:,patmo_idx_SO)*n(:,patmo_idx_O2) &
        - krate(:,17)*n(:,patmo_idx_SO)*n(:,patmo_idx_OH) &
        - krate(:,18)*n(:,patmo_idx_SO)*n(:,patmo_idx_NO2) 
            budget(:,patmo_idx_SO,5) =  + krate(:,12)*n(:,patmo_idx_SH)*n(:,patmo_idx_O) &
        + krate(:,13)*n(:,patmo_idx_SH)*n(:,patmo_idx_O2) 
     endif

     
            dn(:,patmo_idx_OH) = 0d0
     !   - krate(:,1)*n(:,patmo_idx_COS)*n(:,patmo_idx_OH) &
     !   - krate(:,3)*n(:,patmo_idx_CS2)*n(:,patmo_idx_OH) &
     !   - krate(:,8)*n(:,patmo_idx_H2S)*n(:,patmo_idx_OH) &
     !   + krate(:,9)*n(:,patmo_idx_H2S)*n(:,patmo_idx_O) &
     !   + krate(:,13)*n(:,patmo_idx_SH)*n(:,patmo_idx_O2) &
     !   - krate(:,17)*n(:,patmo_idx_SO)*n(:,patmo_idx_OH) &
     !   - krate(:,21)*n(:,patmo_idx_S)*n(:,patmo_idx_OH) &
     !   + krate(:,22)*n(:,patmo_idx_SO2)*n(:,patmo_idx_HO2) &
     !   + krate(:,24)*n(:,patmo_idx_HSO)*n(:,patmo_idx_O2) &
     !   - krate(:,29)*n(:,patmo_idx_SO2)*n(:,patmo_idx_OH) &
     !   + krate(:,31)*n(:,patmo_idx_H2SO4) &
     !   + krate(:,31)*n(:,patmo_idx_H2SO4) &
     !   - krate(:,36)*n(:,patmo_idx_CH3SCH3)*n(:,patmo_idx_OH) &
     !   - krate(:,37)*n(:,patmo_idx_CH3SCH3)*n(:,patmo_idx_OH) &
     !   + krate(:,46)*n(:,patmo_idx_CO2)*n(:,patmo_idx_SH) &
     !   + krate(:,48)*n(:,patmo_idx_SH)*n(:,patmo_idx_COS) &
     !   + krate(:,53)*n(:,patmo_idx_H2O)*n(:,patmo_idx_SH) &
     !   - krate(:,54)*n(:,patmo_idx_OH)*n(:,patmo_idx_SH) &
     !   - krate(:,58)*n(:,patmo_idx_OH)*n(:,patmo_idx_SO) &
     !   + krate(:,62)*n(:,patmo_idx_SO2)*n(:,patmo_idx_H) &
     !   + krate(:,66)*n(:,patmo_idx_H)*n(:,patmo_idx_SO) &
     !   - krate(:,67)*n(:,patmo_idx_OH)*n(:,patmo_idx_SO3) &
     !   - krate(:,69)*n(:,patmo_idx_SO2)*n(:,patmo_idx_OH) &
     !   + krate(:,74)*n(:,patmo_idx_HSO3) &
     !   - krate(:,76)*n(:,patmo_idx_SO2)*n(:,patmo_idx_OH)*n(:,patmo_idx_OH) &
     !   - krate(:,76)*n(:,patmo_idx_SO2)*n(:,patmo_idx_OH)*n(:,patmo_idx_OH) &
     !   + krate(:,81)*n(:,patmo_idx_SO2) &
     !   + krate(:,82)*n(:,patmo_idx_SO2)*n(:,patmo_idx_CH4O3S)

    dn(:,patmo_idx_O) = 0d0
     !   - krate(:,2)*n(:,patmo_idx_COS)*n(:,patmo_idx_O) &
     !   - krate(:,4)*n(:,patmo_idx_CS2)*n(:,patmo_idx_O) &
     !   + krate(:,5)*n(:,patmo_idx_CS)*n(:,patmo_idx_O2) &
     !   - krate(:,7)*n(:,patmo_idx_CS)*n(:,patmo_idx_O) &
     !   - krate(:,9)*n(:,patmo_idx_H2S)*n(:,patmo_idx_O) &
     !   - krate(:,12)*n(:,patmo_idx_SH)*n(:,patmo_idx_O) &
     !   + krate(:,16)*n(:,patmo_idx_SO)*n(:,patmo_idx_O2) &
     !   + krate(:,19)*n(:,patmo_idx_S)*n(:,patmo_idx_O2) &
     !   - krate(:,28)*n(:,patmo_idx_SO2)*n(:,patmo_idx_O) &
     !   - krate(:,32)*n(:,patmo_idx_O)*n(:,patmo_idx_O2) &
     !   - krate(:,35)*n(:,patmo_idx_CH3SCH3)*n(:,patmo_idx_O) &
     !   + krate(:,39)*n(:,patmo_idx_O3) &
     !   + krate(:,40)*n(:,patmo_idx_O2) &
     !   + krate(:,40)*n(:,patmo_idx_O2) &
     !   + krate(:,42)*n(:,patmo_idx_SO3) &
     !   + krate(:,43)*n(:,patmo_idx_SO2) &
     !   + krate(:,45)*n(:,patmo_idx_SO) &
     !   + krate(:,47)*n(:,patmo_idx_CO)*n(:,patmo_idx_SO) &
     !   + krate(:,49)*n(:,patmo_idx_CS)*n(:,patmo_idx_SO) &
     !   - krate(:,50)*n(:,patmo_idx_COS)*n(:,patmo_idx_O) &
     !   + krate(:,52)*n(:,patmo_idx_CO)*n(:,patmo_idx_S) &
     !   + krate(:,54)*n(:,patmo_idx_OH)*n(:,patmo_idx_SH) &
     !   + krate(:,57)*n(:,patmo_idx_H)*n(:,patmo_idx_SO) &
     !   - krate(:,61)*n(:,patmo_idx_SO2)*n(:,patmo_idx_O) &
     !   - krate(:,64)*n(:,patmo_idx_SO)*n(:,patmo_idx_O) &
     !   + krate(:,73)*n(:,patmo_idx_SO3) &
     !   + krate(:,77)*n(:,patmo_idx_O3) &
     !   + krate(:,80)*n(:,patmo_idx_SO2)
    
    if (end_of_run) then
        budget(:,patmo_idx_H2SO4,1) = &
                    + krate(:,30)*n(:,patmo_idx_SO3)*n(:,patmo_idx_H2O) &
                    + krate(:,76)*n(:,patmo_idx_SO2)*n(:,patmo_idx_OH)*n(:,patmo_idx_OH)
            
        budget(:,patmo_idx_H2SO4,2) = &     
                     + krate(:,31)*n(:,patmo_idx_H2SO4) &
                     + krate(:,75)*n(:,patmo_idx_H2SO4) 
     endif

    dn(:,patmo_idx_H2SO4) = &
        + krate(:,30)*n(:,patmo_idx_SO3)*n(:,patmo_idx_H2O) &
        - krate(:,31)*n(:,patmo_idx_H2SO4) &
        - krate(:,75)*n(:,patmo_idx_H2SO4) &
        + krate(:,76)*n(:,patmo_idx_SO2)*n(:,patmo_idx_OH)*n(:,patmo_idx_OH)

    dn(:,patmo_idx_NO2) = 0d0
     !   - krate(:,18)*n(:,patmo_idx_SO)*n(:,patmo_idx_NO2) &
     !   + krate(:,63)*n(:,patmo_idx_SO2)*n(:,patmo_idx_NO)

    dn(:,patmo_idx_SO4) = &
        + krate(:,34)*n(:,patmo_idx_SO2) &
        - krate(:,79)*n(:,patmo_idx_SO4)

    dn(:,patmo_idx_S) = &
        + krate(:,7)*n(:,patmo_idx_CS)*n(:,patmo_idx_O) &
        - krate(:,19)*n(:,patmo_idx_S)*n(:,patmo_idx_O2) &
        - krate(:,20)*n(:,patmo_idx_S)*n(:,patmo_idx_O3) &
        - krate(:,21)*n(:,patmo_idx_S)*n(:,patmo_idx_OH) &
        + krate(:,38)*n(:,patmo_idx_COS) &
        + krate(:,41)*n(:,patmo_idx_CS2) &
        + krate(:,45)*n(:,patmo_idx_SO) &
        - krate(:,52)*n(:,patmo_idx_CO)*n(:,patmo_idx_S) &
        + krate(:,64)*n(:,patmo_idx_SO)*n(:,patmo_idx_O) &
        + krate(:,65)*n(:,patmo_idx_O2)*n(:,patmo_idx_SO) &
        + krate(:,66)*n(:,patmo_idx_H)*n(:,patmo_idx_SO)

    if (end_of_run) then 
            budget(:,patmo_idx_CH3SCH3,1) = &
                 + krate(:,80)*n(:,patmo_idx_SO2) &
                 + krate(:,81)*n(:,patmo_idx_SO2) &
                 + krate(:,82)*n(:,patmo_idx_SO2)*n(:,patmo_idx_CH4O3S)
            
            budget(:,patmo_idx_CH3SCH3,2) = & 
                 + krate(:,35)*n(:,patmo_idx_CH3SCH3)*n(:,patmo_idx_O) &
                 + krate(:,36)*n(:,patmo_idx_CH3SCH3)*n(:,patmo_idx_OH) &
                 + krate(:,37)*n(:,patmo_idx_CH3SCH3)*n(:,patmo_idx_OH)  

    endif
 

    dn(:,patmo_idx_CH3SCH3) = &
        - krate(:,35)*n(:,patmo_idx_CH3SCH3)*n(:,patmo_idx_O) &
        - krate(:,36)*n(:,patmo_idx_CH3SCH3)*n(:,patmo_idx_OH) &
        - krate(:,37)*n(:,patmo_idx_CH3SCH3)*n(:,patmo_idx_OH) 
        + krate(:,80)*n(:,patmo_idx_SO2) &
        + krate(:,81)*n(:,patmo_idx_SO2) &
        + krate(:,82)*n(:,patmo_idx_SO2)*n(:,patmo_idx_CH4O3S)

    if (end_of_run) then
      budget(:,patmo_idx_SO2,1) =  &
        + krate(:,15)*n(:,patmo_idx_SO)*n(:,patmo_idx_O3) &
        + krate(:,16)*n(:,patmo_idx_SO)*n(:,patmo_idx_O2) &
        + krate(:,17)*n(:,patmo_idx_SO)*n(:,patmo_idx_OH) &
        + krate(:,18)*n(:,patmo_idx_SO)*n(:,patmo_idx_NO2)& 
        + krate(:,24)*n(:,patmo_idx_HSO)*n(:,patmo_idx_O2) &
        + krate(:,26)*n(:,patmo_idx_HSO2)*n(:,patmo_idx_O2) &
        + krate(:,31)*n(:,patmo_idx_H2SO4) &
        + (krate(:,35)*n(:,patmo_idx_CH3SCH3)*n(:,patmo_idx_O)) &
        + (krate(:,36)*n(:,patmo_idx_CH3SCH3)*n(:,patmo_idx_OH)) &
        + (krate(:,37)*n(:,patmo_idx_CH3SCH3)*n(:,patmo_idx_OH)) &
        + krate(:,67)*n(:,patmo_idx_OH)*n(:,patmo_idx_SO3) &
        + krate(:,68)*n(:,patmo_idx_SO3)*n(:,patmo_idx_O2) &
        + krate(:,73)*n(:,patmo_idx_SO3) &
        + krate(:,74)*n(:,patmo_idx_HSO3) &
        + krate(:,79)*n(:,patmo_idx_SO4) 
        



        !+ krate(:,42)*n(:,patmo_idx_SO3) &  
    
      budget(:,patmo_idx_SO2,2) = &
         - krate(:,22)*n(:,patmo_idx_SO2)*n(:,patmo_idx_HO2) &
        - krate(:,23)*n(:,patmo_idx_SO2)*n(:,patmo_idx_O3) &
        - krate(:,28)*n(:,patmo_idx_SO2)*n(:,patmo_idx_O) &
        - krate(:,29)*n(:,patmo_idx_SO2)*n(:,patmo_idx_OH) &
        - krate(:,43)*n(:,patmo_idx_SO2) &
        - krate(:,60)*n(:,patmo_idx_SO2)*n(:,patmo_idx_O2) &
        - krate(:,61)*n(:,patmo_idx_SO2)*n(:,patmo_idx_O) &
        - krate(:,62)*n(:,patmo_idx_SO2)*n(:,patmo_idx_H) &
        - krate(:,63)*n(:,patmo_idx_SO2)*n(:,patmo_idx_NO) &
        - krate(:,69)*n(:,patmo_idx_SO2)*n(:,patmo_idx_OH) &
        - krate(:,71)*n(:,patmo_idx_HO2)*n(:,patmo_idx_SO2) &
        - krate(:,76)*n(:,patmo_idx_SO2)*n(:,patmo_idx_OH)*n(:,patmo_idx_OH) &
        - krate(:,80)*n(:,patmo_idx_SO2) &
        - krate(:,81)*n(:,patmo_idx_SO2) &
        - krate(:,82)*n(:,patmo_idx_SO2)*n(:,patmo_idx_CH4O3S)


      budget(:,patmo_idx_SO2,8) = &
              - krate(:,43)*n(:,patmo_idx_SO2) &
         - krate(:,80)*n(:,patmo_idx_SO2) &
        - krate(:,81)*n(:,patmo_idx_SO2)

       

      
      !  + krate(:,24)*n(:,patmo_idx_HSO)*n(:,patmo_idx_O2) &
      !  + krate(:,69)*n(:,patmo_idx_SO2)*n(:,patmo_idx_OH)
      !  - krate(:,22)*n(:,patmo_idx_SO2)*n(:,patmo_idx_HO2) &
      !  - krate(:,23)*n(:,patmo_idx_SO2)*n(:,patmo_idx_O3) &
      !  - krate(:,28)*n(:,patmo_idx_SO2)*n(:,patmo_idx_O) &
      !  - krate(:,29)*n(:,patmo_idx_SO2)*n(:,patmo_idx_OH) &
            !  - krate(:,43)*n(:,patmo_idx_SO2) &
      !  - krate(:,60)*n(:,patmo_idx_SO2)*n(:,patmo_idx_O2) &
      !  - krate(:,61)*n(:,patmo_idx_SO2)*n(:,patmo_idx_O) &
      !  - krate(:,62)*n(:,patmo_idx_SO2)*n(:,patmo_idx_H) &
      !  - krate(:,63)*n(:,patmo_idx_SO2)*n(:,patmo_idx_NO) &
      !  - krate(:,69)*n(:,patmo_idx_SO2)*n(:,patmo_idx_OH) &
      !  - krate(:,71)*n(:,patmo_idx_HO2)*n(:,patmo_idx_SO2) &
      !  - krate(:,76)*n(:,patmo_idx_SO2)*n(:,patmo_idx_OH)*n(:,patmo_idx_OH) &
      !  - krate(:,80)*n(:,patmo_idx_SO2) &
      !  - krate(:,81)*n(:,patmo_idx_SO2) &
      !  - krate(:,82)*n(:,patmo_idx_SO2)*n(:,patmo_idx_CH4O3S)
        

     endif
    
    
    dn(:,patmo_idx_SO2) = &
        + krate(:,15)*n(:,patmo_idx_SO)*n(:,patmo_idx_O3) &
        + krate(:,16)*n(:,patmo_idx_SO)*n(:,patmo_idx_O2) &
        + krate(:,17)*n(:,patmo_idx_SO)*n(:,patmo_idx_OH) &
        + krate(:,18)*n(:,patmo_idx_SO)*n(:,patmo_idx_NO2) &
        - krate(:,22)*n(:,patmo_idx_SO2)*n(:,patmo_idx_HO2) &
        - krate(:,23)*n(:,patmo_idx_SO2)*n(:,patmo_idx_O3) &
        + krate(:,24)*n(:,patmo_idx_HSO)*n(:,patmo_idx_O2) &
        + krate(:,26)*n(:,patmo_idx_HSO2)*n(:,patmo_idx_O2) &
        - krate(:,28)*n(:,patmo_idx_SO2)*n(:,patmo_idx_O) &
        - krate(:,29)*n(:,patmo_idx_SO2)*n(:,patmo_idx_OH) &
        + krate(:,31)*n(:,patmo_idx_H2SO4) &
        - krate(:,34)*n(:,patmo_idx_SO2) &
        + (krate(:,35)*n(:,patmo_idx_CH3SCH3)*n(:,patmo_idx_O)) &
        + (krate(:,36)*n(:,patmo_idx_CH3SCH3)*n(:,patmo_idx_OH)) &
        + (krate(:,37)*n(:,patmo_idx_CH3SCH3)*n(:,patmo_idx_OH)) &
        + krate(:,42)*n(:,patmo_idx_SO3) &
        - krate(:,43)*n(:,patmo_idx_SO2) &
        - krate(:,60)*n(:,patmo_idx_SO2)*n(:,patmo_idx_O2) &
        - krate(:,61)*n(:,patmo_idx_SO2)*n(:,patmo_idx_O) &
        - krate(:,62)*n(:,patmo_idx_SO2)*n(:,patmo_idx_H) &
        - krate(:,63)*n(:,patmo_idx_SO2)*n(:,patmo_idx_NO) &
        + krate(:,67)*n(:,patmo_idx_OH)*n(:,patmo_idx_SO3) &
        + krate(:,68)*n(:,patmo_idx_SO3)*n(:,patmo_idx_O2) &
        - krate(:,69)*n(:,patmo_idx_SO2)*n(:,patmo_idx_OH) &
        - krate(:,71)*n(:,patmo_idx_HO2)*n(:,patmo_idx_SO2) &
        + krate(:,73)*n(:,patmo_idx_SO3) &
        + krate(:,74)*n(:,patmo_idx_HSO3) &
        - krate(:,76)*n(:,patmo_idx_SO2)*n(:,patmo_idx_OH)*n(:,patmo_idx_OH) &
        + krate(:,79)*n(:,patmo_idx_SO4) 
        - krate(:,80)*n(:,patmo_idx_SO2) &
        - krate(:,81)*n(:,patmo_idx_SO2) &
        - krate(:,82)*n(:,patmo_idx_SO2)*n(:,patmo_idx_CH4O3S)

    dn(:,patmo_idx_CH4O3S) = 0d0
     !   + krate(:,37)*n(:,patmo_idx_CH3SCH3)*n(:,patmo_idx_OH) &
     !   - krate(:,82)*n(:,patmo_idx_SO2)*n(:,patmo_idx_CH4O3S)

    dn(:,patmo_idx_HSO3) = &
        - krate(:,27)*n(:,patmo_idx_HSO3)*n(:,patmo_idx_O2) &
        + krate(:,29)*n(:,patmo_idx_SO2)*n(:,patmo_idx_OH) &
        + krate(:,72)*n(:,patmo_idx_HO2)*n(:,patmo_idx_SO3) &
        - krate(:,74)*n(:,patmo_idx_HSO3)
     if (end_of_run) then
       budget(:,patmo_idx_HSO3,1) =  krate(:,29)*n(:,patmo_idx_SO2)*n(:,patmo_idx_OH) 
       budget(:,patmo_idx_HSO3,2) =  n(:,patmo_idx_HSO3)
       budget(:,patmo_idx_HSO3,3) =  krate(:,27)*n(:,patmo_idx_O2)!*n(:,patmo_idx_HSO3)
       budget(:,patmo_idx_HSO3,4) = krate(:,29)

     endif
    dn(:,patmo_idx_H2) = 0d0
     !   + krate(:,10)*n(:,patmo_idx_H2S)*n(:,patmo_idx_H) &
     !   - krate(:,55)*n(:,patmo_idx_H2)*n(:,patmo_idx_SH)
    
     if (end_of_run) then
            budget(:,patmo_idx_H2S,1) = &
                + krate(:,53)*n(:,patmo_idx_H2O)*n(:,patmo_idx_SH) &
                + krate(:,54)*n(:,patmo_idx_OH)*n(:,patmo_idx_SH) &
                + krate(:,55)*n(:,patmo_idx_H2)*n(:,patmo_idx_SH) &
                + krate(:,56)*n(:,patmo_idx_H2O)*n(:,patmo_idx_HSO)
            
            budget(:,patmo_idx_H2S,2) = &
                + krate(:,8)*n(:,patmo_idx_H2S)*n(:,patmo_idx_OH) &
                + krate(:,9)*n(:,patmo_idx_H2S)*n(:,patmo_idx_O) &
                + krate(:,10)*n(:,patmo_idx_H2S)*n(:,patmo_idx_H) &
                + krate(:,11)*n(:,patmo_idx_H2S)*n(:,patmo_idx_HO2) &
                + krate(:,44)*n(:,patmo_idx_H2S)
      
            budget(:,patmo_idx_H2S,8) = & 
                + krate(:,44)*n(:,patmo_idx_H2S)
                
      endif

    dn(:,patmo_idx_H2S) = &
        - krate(:,8)*n(:,patmo_idx_H2S)*n(:,patmo_idx_OH) &
        - krate(:,9)*n(:,patmo_idx_H2S)*n(:,patmo_idx_O) &
        - krate(:,10)*n(:,patmo_idx_H2S)*n(:,patmo_idx_H) &
        - krate(:,11)*n(:,patmo_idx_H2S)*n(:,patmo_idx_HO2) &
        - krate(:,44)*n(:,patmo_idx_H2S) &
        + krate(:,53)*n(:,patmo_idx_H2O)*n(:,patmo_idx_SH) &
        + krate(:,54)*n(:,patmo_idx_OH)*n(:,patmo_idx_SH) &
        + krate(:,55)*n(:,patmo_idx_H2)*n(:,patmo_idx_SH) &
        + krate(:,56)*n(:,patmo_idx_H2O)*n(:,patmo_idx_HSO)

    dn(:,patmo_idx_SH) = &
        + krate(:,1)*n(:,patmo_idx_COS)*n(:,patmo_idx_OH) &
        + 1.17*(krate(:,3)*n(:,patmo_idx_CS2)*n(:,patmo_idx_OH)) &
        + krate(:,8)*n(:,patmo_idx_H2S)*n(:,patmo_idx_OH) &
        + krate(:,9)*n(:,patmo_idx_H2S)*n(:,patmo_idx_O) &
        + krate(:,10)*n(:,patmo_idx_H2S)*n(:,patmo_idx_H) &
        - krate(:,12)*n(:,patmo_idx_SH)*n(:,patmo_idx_O) &
        - krate(:,13)*n(:,patmo_idx_SH)*n(:,patmo_idx_O2) &
        - krate(:,14)*n(:,patmo_idx_SH)*n(:,patmo_idx_O3) &
        + krate(:,25)*n(:,patmo_idx_HSO)*n(:,patmo_idx_O3) &
        + krate(:,44)*n(:,patmo_idx_H2S) &
        - krate(:,46)*n(:,patmo_idx_CO2)*n(:,patmo_idx_SH) &
        - krate(:,48)*n(:,patmo_idx_SH)*n(:,patmo_idx_COS) &
        - krate(:,53)*n(:,patmo_idx_H2O)*n(:,patmo_idx_SH) &
        - krate(:,54)*n(:,patmo_idx_OH)*n(:,patmo_idx_SH) &
        - krate(:,55)*n(:,patmo_idx_H2)*n(:,patmo_idx_SH) &
        + krate(:,57)*n(:,patmo_idx_H)*n(:,patmo_idx_SO) &
        + krate(:,58)*n(:,patmo_idx_OH)*n(:,patmo_idx_SO) &
        + krate(:,59)*n(:,patmo_idx_HSO)*n(:,patmo_idx_O2) &
        - krate(:,70)*n(:,patmo_idx_O2)*n(:,patmo_idx_O2)*n(:,patmo_idx_SH)

    dn(:,patmo_idx_CS) = &
        + krate(:,4)*n(:,patmo_idx_CS2)*n(:,patmo_idx_O) &
        - krate(:,5)*n(:,patmo_idx_CS)*n(:,patmo_idx_O2) &
        - krate(:,6)*n(:,patmo_idx_CS)*n(:,patmo_idx_O3) &
        - krate(:,7)*n(:,patmo_idx_CS)*n(:,patmo_idx_O) &
        + krate(:,41)*n(:,patmo_idx_CS2) &
        - krate(:,49)*n(:,patmo_idx_CS)*n(:,patmo_idx_SO) &
        + krate(:,50)*n(:,patmo_idx_COS)*n(:,patmo_idx_O) &
        + krate(:,51)*n(:,patmo_idx_COS)*n(:,patmo_idx_O2) &
        + krate(:,52)*n(:,patmo_idx_CO)*n(:,patmo_idx_S)

    dn(:,patmo_idx_H) = 0d0
     !   - krate(:,10)*n(:,patmo_idx_H2S)*n(:,patmo_idx_H) &
     !   + krate(:,12)*n(:,patmo_idx_SH)*n(:,patmo_idx_O) &
     !   + krate(:,17)*n(:,patmo_idx_SO)*n(:,patmo_idx_OH) &
     !   + krate(:,21)*n(:,patmo_idx_S)*n(:,patmo_idx_OH) &
     !   + krate(:,44)*n(:,patmo_idx_H2S) &
     !   + krate(:,55)*n(:,patmo_idx_H2)*n(:,patmo_idx_SH) &
     !   - krate(:,57)*n(:,patmo_idx_H)*n(:,patmo_idx_SO) &
     !   - krate(:,62)*n(:,patmo_idx_SO2)*n(:,patmo_idx_H) &
     !   - krate(:,66)*n(:,patmo_idx_H)*n(:,patmo_idx_SO)

    dn(:,patmo_idx_O3) = 0d0
     !   - krate(:,6)*n(:,patmo_idx_CS)*n(:,patmo_idx_O3) &
     !   - krate(:,14)*n(:,patmo_idx_SH)*n(:,patmo_idx_O3) &
     !   - krate(:,15)*n(:,patmo_idx_SO)*n(:,patmo_idx_O3) &
     !   - krate(:,20)*n(:,patmo_idx_S)*n(:,patmo_idx_O3) &
     !   - krate(:,23)*n(:,patmo_idx_SO2)*n(:,patmo_idx_O3) &
     !   - krate(:,25)*n(:,patmo_idx_HSO)*n(:,patmo_idx_O3) &
     !   + krate(:,32)*n(:,patmo_idx_O)*n(:,patmo_idx_O2) &
     !   - krate(:,39)*n(:,patmo_idx_O3) &
     !   + krate(:,51)*n(:,patmo_idx_COS)*n(:,patmo_idx_O2) &
     !   + krate(:,59)*n(:,patmo_idx_HSO)*n(:,patmo_idx_O2) &
     !   + krate(:,60)*n(:,patmo_idx_SO2)*n(:,patmo_idx_O2) &
     !   + krate(:,65)*n(:,patmo_idx_O2)*n(:,patmo_idx_SO) &
     !   + krate(:,68)*n(:,patmo_idx_SO3)*n(:,patmo_idx_O2) &
     !   + krate(:,70)*n(:,patmo_idx_O2)*n(:,patmo_idx_O2)*n(:,patmo_idx_SH) &
     !   - krate(:,77)*n(:,patmo_idx_O3)

    ngas_hpp(:) = ngas_hp(:)/ngas_p(:)
    ngas_hpz(:) = ngas_hp(:)/ngas(:)
    ngas_hmm(:) = ngas_hm(:)/ngas_m(:)
    ngas_hmz(:) = ngas_hm(:)/ngas(:)
    
    if (end_of_run) then
   
      do j=1,60


        budget(j,patmo_idx_H2SO4,3) = &
          + (k_hp(j)-d_hp(j,patmo_idx_H2SO4)) * ngas_hpp(j) * n_p(j,patmo_idx_H2SO4) &
          - (k_hp(j)+d_hp(j,patmo_idx_H2SO4)) * ngas_hpz(j) * n(j,patmo_idx_H2SO4)
        
        budget(j,patmo_idx_SO4,3) = &
          + (k_hp(j)-d_hp(j,patmo_idx_SO4)) * ngas_hpp(j) * n_p(j,patmo_idx_SO4) &
          - (k_hp(j)+d_hp(j,patmo_idx_SO4)) * ngas_hpz(j) * n(j,patmo_idx_SO4)
        

        
        !budget(j,patmo_idx_COS,4) = &
        !  - (k_hm(j)-d_hm(j,patmo_idx_COS)) * ngas_hmz(j) * n(j,patmo_idx_COS) &
        !  + (k_hm(j)+d_hm(j,patmo_idx_COS)) * ngas_hmm(j) * n_m(j,patmo_idx_COS)
        
          budget(j,patmo_idx_H2SO4,4) = &
          - (k_hm(j)-d_hm(j,patmo_idx_H2SO4)) * ngas_hmz(j) * n(j,patmo_idx_H2SO4) &
          + (k_hm(j)+d_hm(j,patmo_idx_H2SO4)) * ngas_hmm(j) * n_m(j,patmo_idx_H2SO4)
        
          budget(j,patmo_idx_SO4,4) = &
          - (k_hm(j)-d_hm(j,patmo_idx_SO4)) * ngas_hmz(j) * n(j,patmo_idx_SO4) &
          + (k_hm(j)+d_hm(j,patmo_idx_SO4)) * ngas_hmm(j) * n_m(j,patmo_idx_SO4)

      enddo
    endif

   ! if (i==1) then 
   !   do j=1,60
   !     dn(j,i) = dn(j,i) &
   !     + (k_hp(j)-d_hp(j,i)) * ngas_hpp(j) * n_p(j,i) &
   !     - ((k_hp(j)+d_hp(j,i)) * ngas_hpz(j) &
   !     + (k_hm(j)-d_hm(j,i)) * ngas_hmz(j)) * n(j,i) &
   !     + (k_hm(j)+d_hm(j,i)) * ngas_hmm(j) * n_m(j,i)
   !   end do
   ! endif

    do i=8,10
      do j=1,60
        dn(j,i) = dn(j,i) &
          + (k_hp(j)-d_hp(j,i)) * ngas_hpp(j) * n_p(j,i) &
          - ((k_hp(j)+d_hp(j,i)) * ngas_hpz(j) &
          + (k_hm(j)-d_hm(j,i)) * ngas_hmz(j)) * n(j,i) &
          + (k_hm(j)+d_hm(j,i)) * ngas_hmm(j) * n_m(j,i)
      end do
    end do
!     
    do i=16,19
      do j=1,60
        dn(j,i) = dn(j,i) &
          + (k_hp(j)-d_hp(j,i)) * ngas_hpp(j) * n_p(j,i) &
          - ((k_hp(j)+d_hp(j,i)) * ngas_hpz(j) &
          + (k_hm(j)-d_hm(j,i)) * ngas_hmz(j)) * n(j,i) &
          + (k_hm(j)+d_hm(j,i)) * ngas_hmm(j) * n_m(j,i)
      end do
    end do
!     
    do i=24,26
      do j=1,60
        dn(j,i) = dn(j,i) &
          + (k_hp(j)-d_hp(j,i)) * ngas_hpp(j) * n_p(j,i) &
          - ((k_hp(j)+d_hp(j,i)) * ngas_hpz(j) &
          + (k_hm(j)-d_hm(j,i)) * ngas_hmz(j)) * n(j,i) &
          + (k_hm(j)+d_hm(j,i)) * ngas_hmm(j) * n_m(j,i)
      end do
    end do
!     
    do i=28,30
      do j=1,60
        dn(j,i) = dn(j,i) &
          + (k_hp(j)-d_hp(j,i)) * ngas_hpp(j) * n_p(j,i) &
          - ((k_hp(j)+d_hp(j,i)) * ngas_hpz(j) &
          + (k_hm(j)-d_hm(j,i)) * ngas_hmz(j)) * n(j,i) &
          + (k_hm(j)+d_hm(j,i)) * ngas_hmm(j) * n_m(j,i)
      end do
    end do
!     
    !diffusion	 
    do i=1,60
      dn(i,1) = dn(i,1) &
        + (k_hp(i)-d_hp(i,1)) * ngas_hpp(i) * n_p(i,1) &
        - ((k_hp(i)+d_hp(i,1)) * ngas_hpz(i) &
        + (k_hm(i)-d_hm(i,1)) * ngas_hmz(i)) * n(i,1) &
        + (k_hm(i)+d_hm(i,1)) * ngas_hmm(i) * n_m(i,1)
    end do
     
    do  i=1,60  
    dn(i,6) = dn(i,6) &
        + (k_hp(i)-d_hp(i,6)) * ngas_hpp(i) * n_p(i,6) &
        - ((k_hp(i)+d_hp(i,6)) * ngas_hpz(i) &
        + (k_hm(i)-d_hm(i,6)) * ngas_hmz(i)) * n(i,6) &
        + (k_hm(i)+d_hm(i,6)) * ngas_hmm(i) * n_m(i,6)
    end do

    
 
    !emission
    dn(1,patmo_idx_COS) = dn(1,patmo_idx_COS) + 7.231324246d7/1d5           !OCS 0.618 TgS/y
    dn(1,patmo_idx_CS2) = dn(1,patmo_idx_CS2) + 3.940809712d7/1d5           !CS2 1.3Tg/y
    dn(1,patmo_idx_H2S) = dn(1,patmo_idx_H2S) + 1.10d9/1d5           !H2S 903.5-1252.6
    dn(1,patmo_idx_SO2) = dn(1,patmo_idx_SO2) + 5.90d9/1d5           !SO2 7196-7715
    dn(1,patmo_idx_CH3SCH3) = dn(1,patmo_idx_CH3SCH3) + 2.6406d9/1d5   !DMS 22TgS

    if (end_of_run) then  
      !budget(1,patmo_idx_COS,5) = 7.231324246d7/1d5    
      budget(1,patmo_idx_SO2,5) = 5.90d9/1d5
      budget(1,patmo_idx_CS2,5) = 3.940809712d7/1d5 
      budget(1,patmo_idx_H2S,5) = 1.10d9/1d5 
      budget(1,patmo_idx_CH3SCH3,5) = 2.6406d9/1d5  
       
      !budget(1,patmo_idx_COS,6) = (1.585d-8*5.55)*n(1,patmo_idx_COS)   
      budget(1,patmo_idx_SO2,6) = (1.355d-7*5.55)*n(1,patmo_idx_SO2)
      budget(1,patmo_idx_H2S,6) = 0.0*n(1,patmo_idx_H2S)   
      budget(1,patmo_idx_CS2,6) = (5.4d-9*5.55)*n(1,patmo_idx_CS2)   
      budget(1,patmo_idx_CH3SCH3,6) = (0.0)*n(1,patmo_idx_CH3SCH3)
    endif     

    !(Tuning mechanism *5.5 for the lowest layer for recalculated OCS emissions!
    !dry deposition
    dn(1,patmo_idx_COS) = dn(1,patmo_idx_COS) - (1.585d-8*5.55)*n(1,patmo_idx_COS)             !OCS  lifetime 2.0y
    dn(1,patmo_idx_SO2) = dn(1,patmo_idx_SO2) - (1.355d-7*5.55)*n(1,patmo_idx_SO2)            !SO2  lifetime 5.0d
    dn(1,patmo_idx_CS2) = dn(1,patmo_idx_CS2) - (5.4d-9*5.55)*n(1,patmo_idx_CS2)
    !aerosol formation
       do i=13,34
      
      if (va(i) <= n(i,24) .and. pa(i) >= n(i,24)) then
         ! if (end_of_run) then 
         !   WRITE(93,*) i, va(i), pa(i), n(i,24),  n(i,24)-va(i),  dn(i,patmo_idx_H2SO4) 
        !endif
      	if (end_of_run) then 
             budget(j,patmo_idx_H2SO4,5) = - (n(i,24)-va(i))
             budget(j,patmo_idx_SO4,5)   = + (n(i,24)-va(i))
        endif

        dn(i,patmo_idx_H2SO4) = dn(i,patmo_idx_H2SO4) - (n(i,24)-va(i))
        dn(i,patmo_idx_SO4)   = dn(i,patmo_idx_SO4  ) + (n(i,24)-va(i))
      end if
    end do     

    !gravity settling SO4 Aerosol (JAM-Kasten-1968,r=0.3)
    if (end_of_run) then
      budget(60,  patmo_idx_SO4, 6) =  -gd(60)*n(60,patmo_idx_SO4)
      
      do j=59,1,-1
        budget(j,  patmo_idx_SO4, 6) =  gd(j+1)*n(j+1, patmo_idx_SO4) -gd(j)*n(j,patmo_idx_SO4)
      end do  
    endif
    
 !   do j=60,2,-1
 !     dn(j,  patmo_idx_SO4) = dn(j,  patmo_idx_SO4) - gd(j)*n(j,patmo_idx_SO4)
 !     dn(j-1,patmo_idx_SO4) = dn(j-1,patmo_idx_SO4) + gd(j)*n(j,patmo_idx_SO4)
 !     
 !   end do  
 !   dn(1,patmo_idx_SO4) = dn(1,patmo_idx_SO4)-gd(1)*n(1,patmo_idx_SO4)

    !wet deposition
    if (end_of_run) then
     ! budget(12,patmo_idx_COS,7) = -wetdep(12,patmo_idx_COS)*n(12,patmo_idx_COS)
      budget(12,patmo_idx_SO2,7) = -wetdep(12,patmo_idx_SO2)*n(12,patmo_idx_SO2)
      budget(12,patmo_idx_CS2,7) = -wetdep(12,patmo_idx_CS2)*n(12,patmo_idx_CS2)
      budget(12,patmo_idx_H2SO4,7) = -wetdep(12,patmo_idx_H2SO4)*n(12,patmo_idx_H2SO4)
      budget(12,patmo_idx_SO4,7) = -wetdep(12,patmo_idx_SO4)*n(12,patmo_idx_SO4)
      budget(12,patmo_idx_H2S,7) =  -wetdep(12,patmo_idx_H2S)*n(12,patmo_idx_H2S)
      budget(12,patmo_idx_CH3SCH3,7) =  -wetdep(12,patmo_idx_CH3SCH3)*n(12,patmo_idx_CH3SCH3)
      


      do j=11,1,-1
      !  budget(j,patmo_idx_COS,7) = wetdep(j+1,patmo_idx_COS)*n(j+1,patmo_idx_COS) - wetdep(j,patmo_idx_COS)*n(j,patmo_idx_COS)
        budget(j,patmo_idx_SO2,7) = wetdep(j+1,patmo_idx_SO2)*n(j+1,patmo_idx_SO2) - wetdep(j,patmo_idx_SO2)*n(j,patmo_idx_SO2)
        budget(j,patmo_idx_CS2,7) = wetdep(j+1,patmo_idx_CS2)*n(j+1,patmo_idx_CS2) - wetdep(j,patmo_idx_CS2)*n(j,patmo_idx_CS2)
        budget(j,patmo_idx_H2S,7) =  wetdep(j+1,patmo_idx_H2S)*n(j+1,patmo_idx_H2S)  - wetdep(j,patmo_idx_H2S)*n(j,patmo_idx_H2S)

        budget(j,patmo_idx_CH3SCH3,7) = wetdep(j+1,patmo_idx_CH3SCH3)*n(j+1,patmo_idx_CH3SCH3) - wetdep(j,patmo_idx_CH3SCH3)*n(j,patmo_idx_CH3SCH3)

        budget(j,patmo_idx_SO4,7) = wetdep(j+1,patmo_idx_SO4)*n(j+1,patmo_idx_SO4) - wetdep(j,patmo_idx_SO4)*n(j,patmo_idx_SO4)
        budget(j,patmo_idx_H2SO4,7) = wetdep(j+1,patmo_idx_H2SO4)*n(j+1,patmo_idx_H2SO4)-wetdep(j,patmo_idx_H2SO4)*n(j,patmo_idx_H2SO4)
      enddo
    endif
    
    do j=12,2,-1
      do i=1,chemSpeciesNumber
        dn(j,i) = dn(j,i)-wetdep(j,i)*n(j,i)
        dn(j-1,i) = dn(j-1,i)+wetdep(j,i)*n(j,i)
      end do   
    end do

    do i=1,chemSpeciesNumber
      dn(1,i) = dn(1,i)-wetdep(1,i)*n(1,i)
    end do  

    !DMS → SO2 (96%)
    ! do i=1,60
    !     dn(i,patmo_idx_CH3SCH3) = dn(i,patmo_idx_CH3SCH3)- (n(i,30)*0.96d0)
    !     dn(i,patmo_idx_SO2) = dn(i,patmo_idx_SO2)+ (n(i,30)*0.96d0)
    ! end do
    ! end if

    !unroll chemistry
    dy(:) = 0d0
    do i=1,speciesNumber
      dy((i-1)*cellsNumber+1:(i*cellsNumber)) = dn(:,i)
    end do

  end subroutine fex
end module patmo_ode
