      subroutine pup
      
!!    ~ ~ ~ PURPOSE ~ ~ ~
!!    this subroutine is called by grow.f and calculates plant pesticide uptake

!!    ~ ~ ~ INCOMING VARIABLES ~ ~ ~
!!    name        |units          |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    curyr       |none           |current year of simulation
!!    nyskip      |none           |number of years to skip output summarization/
!!                                |printing
!!    ihru        |none           |HRU number
!!    hru_dafr(:) |km**2/km**2    |fraction of watershed area in HRU
!!    idplt(:)    |none           |land cover code from crop.dat
!!    npmx        |none           |number of different pesticides used in
!!                                |the simulation
!!    npno(:)     |none           |array of unique pesticides used in watershed
!!    pstuptk(:)  |none           |plant uptake coefficient per pesticide
!!                                |This parameter controls the amount of
!!                                |pesticide removed from the soil pore water
!!                                |by the plant and which accumulates in the
!!                                |plant  !! JK: Not yet active
!!    wuse(:)     |mm H2O         |water uptake by plants in each soil layer
!!    sol_pst(:,:,:)|kg/ha        |amount of pesticide in layer
!!    sol_pst_conc(:,:,:)|kg/mm-ha|concentration of pesticide in each soil layer
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ OUTGOING VARIABLES ~ ~ ~
!!    name        |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    plt_pst2(:,:) |kg/ha       |amount of pesticide stored in plant (hru)
!!                               |(plt_pst is amount of pesticide stored ON plant)
!!    plt_upst_d(:,:)|kg/ha      |plant uptake of pesticide in HRU for the 
!!                               |time step
!!    plt_upst_ly(:,:)|kg/ha     |plant uptake of pesticide in layer for the 
!!                               |current HRU and time step
!!    solpst_tot(:,:)|kg/ha      |amount of pesticide stored in soil column (hru)
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
     
!!    ~ ~ ~ LOCAL DEFINITIONS ~ ~ ~
!!    name        |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    j           |none          |HRU number
!!    ly          |none          |counter (soil layers)
!!    k           |none          |counter (pesticides)
!!    kk          |none          |pesticide number from pest.dat
!!    yy          |kg/ha         |amount of pesticide removed from layer
!!    icrop       |none          |land cover code
!!    pltpstuptk_lst|none        |list with plants that take-up pesticide 
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ SUBROUTINES/FUNCTIONS CALLED ~ ~ ~

!!    ~ ~ ~ ~ ~ ~ END SPECIFICATIONS ~ ~ ~ ~ ~ ~

      use parm
      implicit none 

      integer :: j, icrop, ly, k, kk
      real*8 :: yy
      integer :: pltpstuptk_lst(19) !! Hard coded plant codes to which plant uptake applies.
      integer :: tfapst_ids(6) !! Hard coded pesticide codes to which plant uptake applies

      j = 0
      j = ihru

      !if (subnum(j)=='00001' .and. hruno(j)=='0092' .and. iyr==2014 .and. iida==147) then
      !  write(*,*) iyr, iida;
      !endif

      pltpstuptk_lst = (/3,2,77,31,19,72,20,84,32,63,70,128,30,21,69,83,27,99,28/)
      tfapst_ids = (/236,238,239,240,241,242/)  !! can be removed if read from pest.dat

      if (hrupest(j) /= 0) then
        !! define agricultural crops for which pesticide uptake takes place
        icrop = 0.
        icrop = idplt(j)
        
        if (ANY(pltpstuptk_lst == icrop)) then !!check if plant uptake applies to crop
          do k = 1, npmx
            kk = 0
            kk = npno(k)
            !pstuptk(k) = 0.305  !! to be removed once read from pest.dat            
            if (ANY(tfapst_ids == kk)) then  !! only run if pesticide is TFA (can be removed if read from pest.dat)
              !pstuptk = 0.305  !!JK: I think we need to restrict this to TFA?
              if (kk > 0) then  !! only run if pesticide simulated
                plt_upst_d(k,j) = 0.
                do ly = 1, sol_nly(j)            
                  plt_upst_ly(k,ly) = 0.
                  !if (subnum(j)=='00001' .and. hruno(j)=='0092' .and. iyr==2014 .and. iida==147) then
                  !  write(*,*) k,j,ly, sol_pst(k,j,ly), sol_pst_conc(k,j,ly)
                  !endif
                  if (sol_pst(k,j,ly) >= 1.e-9 .and. sol_pst_conc(k,j,ly) >= 1.e-9) then 
                      !! calculate pesticide lost in layer due to plant uptake
                      yy = 0.
                      yy = pstuptk * wuse(ly) * sol_pst_conc(k,j,ly)
                      if (yy > sol_pst(k,j,ly)) yy = sol_pst(k,j,ly)
                      !! update soil layer, soil column and plant pesticide mass
                      sol_pst(k,j,ly) = sol_pst(k,j,ly) - yy
                      plt_upst_d(k,j) = plt_upst_d(k,j) + yy
                      plt_upst_ly(k,ly) = yy
                  end if
                end do
                if (plt_upst_d(k,j) < 0.) plt_upst_d(k,j) = 0.
                !!update amount of pesticide stored in plant (hru)
                plt_pst2(k,j) = plt_pst2(k,j) + plt_upst_d(k,j)  !!Don't reset if plant is killed/harvested. Just show accumulation
              end if !! only run if pesticide simulated
            end if !! !! only run if pesticide is TFA
          end do
        end if !!check if plant uptake applies to crop
      end if

      return
      end

                        !!JK: to be moved to writed.f
                  !!    wshd_pstup(:)|kg/ha        |average annual amount of pesticide 
                  !!                               |plant uptake 
                  !!if (curyr > nyskip) then
                    !! add wshd_pstup(k) = wshd_pstup(k) / yrs to writeaa.f
                    !!wshd_pstup(k) = wshd_pstup(k) + pltpst(j,k) * hru_dafr(j)  
                  !!end if