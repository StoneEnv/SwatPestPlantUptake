      subroutine rtpest3
      
!!     ~ ~ ~ PURPOSE ~ ~ ~
!!     this subroutine computes the daily stream pesticide balance (3)
!!     (soluble and sorbed)     

!!    ~ ~ ~ INCOMING VARIABLES ~ ~ ~
!!    name          |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    ch_l2(:)      |km            |length of main channel
!!    ch_w(2,:)     |m             |average width of main channel !HR potentially use widht at bottom phi(6,k)
!!    chpst_conc(:) |mg/(m**3)     |initial pesticide concentration in reach
!!    chpst_koc(:)  |m**3/g        |pesticide partition coefficient between
!!                                 |water and sediment in reach
!!    chpst_mix(:)  |m/day         |mixing velocity (diffusion/dispersion) for
!!                                 |pesticide in reach
!!    chpst_rea(:)  |1/day         |pesticide reaction coefficient in reach
!!    chpst_rsp(:)  |m/day         |resuspension velocity in reach for pesticide
!!                                 |sorbed to sediment
!!    chpst_stl(:)  |m/day         |settling velocity in reach for pesticide
!!                                 |sorbed to sediment
!!    chpst_vol(:)  |m/day         |pesticide volatilization coefficient in 
!!                                 |reach
!!    drift(:)      |kg            |amount of pesticide drifting onto main
!!                                 |channel in subbasin
!!    hru_sub(:)    |none          |subbasin number where reach is located
!!    inum1         |none          |reach number
!!    inum2         |none          |inflow hydrograph storage location number
!!    rchdep        |m             |depth of flow on day
!!    rchwtr        |m^3 H2O       |water stored in reach at beginning of day
!!    rnum1         |none          |fraction of overland flow
!!    rtwtr         |m^3 H2O       |water leaving reach on day
!!    sedpst_act(:) |m             |depth of active sediment layer in reach for
!!                                 |pesticide
!!    sedpst_bry(:) |m/day         |pesticide burial velocity in river bed
!!                                 |sediment
!!    sedpst_conc(:)|mg/(m**3)     |inital pesticide concentration in river bed
!!                                 |sediment
!!    sedpst_rea(:) |1/day         |pesticide reaction coefficient in river bed
!!                                 |sediment
!!    varoute(11,:) |mg pst        |pesticide in solution
!!    varoute(12,:) |mg pst        |pesticide sorbed to sediment
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ OUTGOING VARIABLES ~ ~ ~
!!    name        |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    bury        |mg pst        |loss of pesticide from active sediment layer
!!                               |by burial
!!    difus       |mg pst        |diffusion of pesticide from sediment to reach
!!    reactb      |mg pst        |amount of pesticide in sediment that is lost
!!                               |through reactions
!!    reactw      |mg pst        |amount of pesticide in reach that is lost
!!                               |through reactions
!!    resuspst    |mg pst        |amount of pesticide moving from sediment to
!!                               |reach due to resuspension
!!    setlpst     |mg pst        |amount of pesticide moving from water to
!!                               |sediment due to settling
!!    solpesto    |mg pst/m^3    |soluble pesticide concentration in outflow
!!                               |on day
!!    sorpesto    |mg pst/m^3    |sorbed pesticide concentration in outflow
!!                               |on day
!!    volatpst    |mg pst        |amount of pesticide in reach lost by
!!                               |volatilization
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ LOCAL DEFINITIONS ~ ~ ~
!!    name        |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    bedvol      |m^3           |volume of river bed sediment
!!    chpstmass   |mg pst        |mass of pesticide in reach
!!    depth       |m             |depth of water in reach
!!    fd2         |
!!    frsol       |none          |fraction of pesticide in reach that is soluble
!!    frsrb       |none          |fraction of pesticide in reach that is sorbed
!!    jrch        |none          |reach number
!!    pstin       |mg pst        |total pesticide transported into reach
!!                               |during time step
!!    sedcon      |g/m^3         |sediment concentration
!!    sedpstmass  |mg pst        |mass of pesticide in bed sediment
!!    solpstin    |mg pst        |soluble pesticide entering reach during 
!!                               |time step
!!    sorpstin    |mg pst        |sorbed pesticide entering reach during
!!                               |time step
!!    tday        |days          |flow duration
!!    wtrin       |m^3 H2O       |volume of water entering reach during time
!!                               |step
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ LOCAL DEFINITIONS, VVWM Methodology ~ ~ ~
!!    name        |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    foc         |unitless      |fraction of organic carbon in sdeiment
!!    Bur         |kg/s          |burial rate of sediment
!!    c1          |kg/m3         |aqueous concentration in water column
!!    c2          |kg/m3         |aqueous concentration in water column
!!    mSed_1      |kg            |mass of sedment in water column
!!    mDoc_1      |kg            |mass of DOC in water column
!!    mBio_1      |kg            |mass of biota in water column
!!    mSed_2      |kg            |mass of sedment in benthic
!!    mDoc_2      |kg            |mass of DOC in benthic
!!    mBio_2      |kg            |mass of biota in benthic
!!    KdSed_1     |m3/kg         |linear isother partition coefficient for sediment in water column
!!    KdDoc_1     |m3/kg         |linear isother partition coefficient for DOC in water column
!!    KdBio_1     |m3/kg         |linear isother partition coefficient for biota in water column
!!    KdSed_2     |m3/kg         |linear isother partition coefficient for sediment in benthic
!!    KdDoc_2     |m3/kg         |linear isother partition coefficient for DOC in benthic
!!    KdBio_2     |m3/kg         |linear isother partition coefficient for biota in benthic
!!    v1          |m3            |volume of water in water column
!!    v2          |m3            |volume of water in benthic
!!    Q           |m3.s          |volumetric flow rate of water out of water column
!!    uHydr       |1/s           |1st order hydrolysis rate
!!    uPhoto      |1/s           |1st order photolysis rate
!!    uVol        |1/s           |1st order effective volatilization rate
!!    uBio_1      |1/s           |1st order biological metabolism rate in water column
!!    uBio_2      |1/s           |1st order biological metabolism rate in benthic
!!    g1          |1/s           |effective loss rate in water column
!!    g2          |1/s           |effective loss rate in benthic
!!    D           |m2/s          |overall water column to benthic dispesion coefficient; VVWM default = 8.33 E-09
!!    dx          |m             |boundary layer thinkess; VVWM deafult = 1.02 m
!!    omega       |1/s           |mass transfer coefficient betwwen benthic and water coluumn
!!    thta        |unitless      |ratio of solute holding capacity in benthic to that of water column
!!    Doc1        |mg/l          |DOC concentration in water column; VVWM default = 5.0 mg/l
!!    Doc2        |mg/l          |DOC concentration in water benthic; VVWM default = 5.0 mg/l
!!    Plmas       |mg/L          |biosolids concentration in water column; VVWM default = 0.4 mg/l
!!    bd2         |g/cm3 or g/mL |bulk density of solids in benthid; VVWM default = 1.35
!!    por2        |unitless      |prorosity of benthic; VVWM default = 0.5
!!    bnmas       |g/m2          |areal concentration of biosolids in benthic; VVWM default = 0.006
!!    bedarea     |m2            |benthic surface area
!!    cap_1       |kg            |solute holding capacity in water column
!!    cap_2       |kg            |solute holding capacity in benthic
!!      fw1       |kg            |fraction of solute in aqueous phase in water column
!!      fw2       |kg            |fraction of solute in aqueous phase in benthic 
!!        A       |              |analytical solution constant A
!!        B       |              |analytical solution constant B
!!        E       |              |analytical solution constant E
!!        F       |              |analytical solution constant F
!!      C1o       |kg/m3         |initial concentration in water column
!!      C2o       |kg/m3         |initial concentration in benthic
!!       Xd       |unitless      |fraction of water column sorbed pesticide transferred to benthic via deposition
!!       Xr       |unitless      |fraction of benthic sorbed pesticide transferred to water column via resuspension
!!       af       |              |variable in analytical solution method
!!      fxa       |              |variable in analytical solution method
!!      bxe       |              |variable in analytical solution method
!!      dif       |              |variable in analytical solution method
!!      bbb       |              |variable in analytical solution method
!!    root1       |              |variable in analytical solution method
!!    root2       |              |variable in analytical solution method
!!       DD       |              |variable in analytical solution method
!!      EE1       |              |variable in analytical solution method
!!       FF       |              |variable in analytical solution method
!!       m1       |kg/m3         |initial aqueous concentration in water column
!!       m2       |kg/m3         |initial aqueous concentration in benthic
!!       X1       |              |variable in analytical solution method
!!       Y1       |              |variable in analytical solution method
!!      rt1       |              |variable in analytical solution method
!!      rt2       |              |variable in analytical solution method
!!    exrt1       |              |variable in analytical solution method
!!    exrt2       |              |variable in analytical solution method
!!      ccc       |              |variable in analytical solution method
!!      ddd       |              |variable in analytical solution method
!!      mn1       |kg/m3         |end of day aqueous concentration in water column
!!      mn2       |kg/m3         |end of day aqueous concentration in benthic
!!       gx       |              |variable in analytical solution method
!!       hx       |              |variable in analytical solution method
!!    term1       |              |variable in analytical solution method
!!    term2       |              |variable in analytical solution method
!!    term3       |              |variable in analytical solution method
!!    term4       |              |variable in analytical solution method
!!    term1       |              |variable in analytical solution method
!!    T_end       |s             |averaging period for concentrations (default is 1 day)
!!    mavg1       |kg/m3         |average aqueous concentration in water column
!!    mavg2       |kg/m3         |average aqueous concentration in benthic
!!   SedRes       |kg/d          |sediment resuspended in reach
!!   SedDep       |kg/d          |sediment deposited in reach 
!!  SedInit       |kg            |inital beginning of timestep sediment in reach
!! SedRteIn       |kg            |toal sediment routed into reach during daily timestep
!!    jrch        |none          |reach number
!!    pstin       |mg pst        |total pesticide transported into reach
!!                               |during time step
!!    sedpstmass  |mg pst        |mass of pesticide in bed sediment
!!    solpstin    |mg pst        |soluble pesticide entering reach during 
!!                               |time step
!!    sorpstin    |mg pst        |sorbed pesticide entering reach during
!!                               |time step
!!    tday        |days          |flow duration
!!    wtrin       |m^3 H2O       |volume of water entering reach during time
!!                               |step
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ SUBROUTINES/FUNCTIONS CALLED ~ ~ ~
!!    Intrinsic: abs

!!    ~ ~ ~ ~ ~ ~ END SPECIFICATIONS ~ ~ ~ ~ ~ ~

      use parm

      implicit none

      integer :: jrch, hr
      real*8 :: solpstin, sorpstin, pstin, depth, chpstmass, frsol, frsrb
      real*8 :: sedpstmass, bedvol, fd2, wtrin, solmax, sedcon, tday
      real*8 :: fr_stored, fr_routed, reactw1, reactw2
    
      real*8 :: foc, Bur, c1, c2, mSed_1, mDoc_1, mBio_1, mSed_2, mDoc_2, mBio_2
      real*8 :: KdSed_1, KdDoc_1, KdBio_1, KdSed_2, KdDoc_2, KdBio_2, v1, v2, Q
      real*8 :: uHydr, uPhoto, uVol, uBio_1, uBio_2, g1, g2, D, dx, omega, thta
      real*8 :: Doc1, Doc2, Plmas, bd2, por2, bnmas, bedarea, cap1, cap2, fw1, fw2
      real*8 :: A, B, E, F, C1o, C2o, Xd, Xr
      real*8 :: af, fxa, bxe, dif, bbb, root1, root2, DD, EE1, FF !!variables used in analytical solution to Sim DiffEq (VVWM method)
      real*8 :: X1, Y1, m1, m2, rt1, rt2, exrt1, exrt2, ccc, ddd, mn1, mn2 !!variables used in analytical solution to Sim DiffEq (VVWM method)
      real*8 :: gx, hx, term1, term2, term3, term4, T_end, mavg1, mavg2 !!variables used in analytical solution to Sim DiffEq (VVWM method)
      real*8 :: SedRes, SedDep, SedInit, SedRteIn

      jrch = 0
      jrch = inum1

!! calculate volume of water entering reach, 1
      wtrin = 0.
      wtrin = varoute(2,inum2) * (1. - rnum1)

!! Calculate sediment flux terms
      SedRes = 0.
      SedDep = 0.
      SedInit = 0.
      SedRes = ((sedrch / (rtwtr / (rchwtr + wtrin))) - sedinorg + rchdy(57,jrch)) * 1000. !!resuspended sediment in kg
      SedDep = rchdy(57,jrch) * 1000. !!despoited sediment in kg
      SedInit = sedinorg * 1000.   !!Convert metric tons to kg
!! Sediment initialization for hourly inflows:
      SedRteIn = (varoute(3,inum2) * (1. - rnum1)) * 1000.  !!sediment routed in
      SedInit = SedInit - SedRteIn  !!set initial sedument to before incoming routed added

!! calculate flow duration
      tday = 0.
      tday = rttime / 24.0
      if (tday > 1.0) tday = 1.0
         
!! pesticide transported into reach during day
      solpstin = 0.
      sorpstin = 0.
      solpstin = varoute(11,inum2) * (1. - rnum1)
      sorpstin = varoute(12,inum2) * (1. - rnum1)

!! MW: Only look at RCH 1
      if (jrch <= 1) then
         jrch = jrch
      endif

!! initialize pesticide out
      solpesto = 0.
      sorpesto = 0.

!! initialize pesticide in reach from previous timestep concentration and beginning of timestep water stored
      chpstmass = 0.
      chpstmass = chpst_conc(jrch) * rchwtr
      
!! initialize pesticide in benthic sediment from previous timestep concentration
      bedarea = 0.
      bedarea = ch_w(2,jrch) * ch_l2(jrch) * 1000.
      bedvol = 0.
      bedvol = bedarea * sedpst_act(jrch)
      sedpstmass = 0.
      sedpstmass = sedpst_conc(jrch) * bedvol

!! Check if any pesticide exists in system
      if (solpstin + sorpstin + chpstmass + sedpstmass < 1.e-12) then  
        chpst_conc(jrch) = 0.
        sedpst_conc(jrch) = 0.
        return
      end if

!! Issues to still consider
!! 1.) travel time in reach less than 1 day; for now think that constant 1-day time for reach processes is appropriate
!! 2.) issue of washout term in simulatneous differential equation solution, versus mass trasported out of reach based on average timestep concentration time or routed water

      if (rtwtr / 86400. > 1.e-6) then

!! initialize constants
          foc = 0.04 !VVWM default
          Doc1 = 5. !VVWM default
          Doc2 = 5. !VVWM default
          Plmas = 0.4 !VVWM default
          por2 = 0.5 !VVWM default
          bd2 = 1.35 !VVWM default
          bnmas = 0.006 !VVWM default
          !!T_end = 86400. !!1-day timestep
          T_end = 3600. !!1-hr timestep (for sub-daily interbal timestep in subroutine)

!! iniitalize Kd based on incoming Koc, convert from mL/g to m3/kg
          KdSed_1 = foc * chpst_koc(jrch) * 0.001
          KdSed_2 = KdSed_1
          KdDoc_1 = 0.2114 * chpst_koc(jrch) * 0.001
          KdDoc_2 = chpst_koc(jrch) * 0.001
          KdBio_1 = 0.436 * ((chpst_koc(jrch) / 0.35)**0.907) * 0.001
          KdBio_2 = KdBio_1

!! initialize volume and mass in water column at start of timestep
!!Initialize in hourly sub-loop now
          !!v1 = rchwtr + (wtrin / 24.)
          !!mSed_1 = SedInit
          !!mDoc_1 = Doc1 * v1 * 0.001
          !!mBio_1 = Plmas * v1 * 0.001

!! initialize volume and mass in benthic at start of timestep
          v2 = bedvol * por2
          mSed_2 = bedvol * bd2 * 1000. 
          mDoc_2 = Doc2 * v2 * 0.001
          mBio_2 = bnmas * bedarea * 0.001

!! Calculate mass transfer coefficient, omega
          D = 0.00000000833
          dx = 1.02   !!May want to vary with depth
          omega = (bedarea * D) / (dx * bedvol)

!! Define water column and benthic disspation rates
          uHydr = 0. !assume no hydrolysis (would need additional SWAT inputs)
          uPhoto = 0. !assume no photolysis (would need additional SWAT inputs)
          uVol = 0. !assume no volatilization (would need additional SWAT inputs)
          uBio_1 = chpst_rea(jrch) / 86400. !convert from 1/d to 1/s
          uBio_2 = sedpst_rea(jrch) / 86400. !convert from 1/d to 1/s
          
!! Can't have burial sediment mass greater then benthic sediment mass
          if (SedDep < mSed_2) then
             Bur = SedDep / 86400. ! convert deposition from kg/day to kg/s
          else
             Bur = mSed_2 / 86400.
          end if
          Bur = 0. !!Set burial equal to zero

!! Start hourly sub-loop here
          Do  hr = 1, 24, 1 

!! Adjust hourly timestep capacities based
          if (hr == 1) then
             v1 = rchwtr + (wtrin / 24.) !!- (rtwtr / 24.)
             mDoc_1 = Doc1 * v1 * 0.001
             mBio_1 = Plmas * v1 * 0.001
             mSed_1 = SedInit + (SedRteIn / 24.) - ((SedDep - SedRes)/ 24.) !!- ((sedrch / 24.) * 1000.)
          else 
             v1 = v1 + (wtrin / 24.) !!- (rtwtr / 24.)
             mDoc_1 = Doc1 * v1 * 0.001
             mBio_1 = Plmas * v1 * 0.001
             mSed_1 = mSed_1 + (SedRteIn / 24.) - ((SedDep - SedRes)/ 24.) !!- ((sedrch / 24.) * 1000.)
          end if 
          
!! Calculate solute holding capacities
          cap1 = KdSed_1 * mSed_1 + KdBio_1 * mBio_1 + KdDoc_1 * mDoc_1 + v1
          cap2 = KdSed_2 * mSed_2 + KdBio_2 * mBio_2 + KdDoc_2 * mDoc_2 + v2
          fw1 = v1 / cap1
          fw2 = v2 / cap2
          thta = cap2 / cap1

!! Calculate effective dissipation rates in water column and benthic
          g1 = (uHydr + uPhoto + uVol) * fw1 + uBio_1
          g2 = (fw2 * uHydr) + uBio_2 + (Bur * KdSed_2 / cap2)

!! Solve analytical solution to differential equations
          A = (-1. * g1) - (omega * thta)
          B = omega * thta
          E = omega
          F = (-1.* g2) - omega

!! Calculate boundary conditions
!! Assumes instantaneous transfer of sorbed pesticide in water column to benthic based on fraction of initial conditions sediment mass that is deposited in reach
!! Likewise, resuspension of benthic into water column occur instantaneously at beginning of timestep, and sorbed pesticide from benthic to water column occurs proportionaly
          if (SedDep > 0.) then
             !!Xd = min(SedDep / SedInit, 1.0)
             Xd = min((SedDep / 24.0) / SedInit, 1.0)
          else
             Xd = 0.
          end if
          if (SedRes > 0.) then
             !!Xr = min(SedRes / mSed_2, 1.0)
             Xr = min((SedRes / 24.0) / mSed_2, 1.0)
          else
             Xr = 0.
          end if
          !!Xd = 0.
          !! Put in some option "prben" like fractions to allocate sorbed to benthic
          !!Xd = 0.5
          !!Xr = 0.0
         
          if (hr <= 24)  then 
       
              C1o = (fw1/v1) * ((solpstin / 24.) + ((1. - Xd) * (sorpstin / 24.)) + chpstmass) * 0.000001 !! water column aqueous conc at beginning of TS, convert mg to kg
              C2o = (fw2/v2) * ((Xd * (sorpstin / 24.)) + sedpstmass) * 0.000001 !! benthic aqueous conc at beginning of TS, convert mg to kg
              !!Now, adjust mass in water column and benthic by the resuspension pesticide mass
              C1o = C1o + (fw1/v1) * (Xr * (1. -  fw2) * sedpstmass * 0.000001)
              C2o = C2o - (fw2/v2) * (Xr * (1. -  fw2) * sedpstmass * 0.000001)
          else
              C1o = (fw1/v1) * (chpstmass - (Xd * (1. - fw1) * chpstmass)) * 0.000001 !! water column aqueous conc at beginning of TS, convert mg to kg
              C2o = (fw2/v2) * (sedpstmass + (Xd * (1. - fw1) * chpstmass)) * 0.000001 !! benthic aqueous conc at beginning of TS, convert mg to kg
              !!Now, adjust mass in water column and benthic by the resuspension pesticide mass
              C1o = C1o + (fw1/v1) * (Xr * (1. -  fw2) * sedpstmass * 0.000001)
              C2o = C2o - (fw2/v2) * (Xr * (1. -  fw2) * sedpstmass * 0.000001)
          end if 

!! Solve simultaneous diff equations
!! Follow VVWM method to for breaking down the equations
          af=A+F
          fxa=F*A
          bxe=B*E
          dif=4*(fxa-bxe)
          bbb=sqrt(af*af-dif)
          root1 = (af+bbb)/2.
          root2 = (af-bbb)/2.    
          DD = (root1-A)/B
          EE1 = (root2-A)/B
          FF = EE1-DD
          m1 = C1o    !! m1 in kg/m3
          m2 = C2o    !! m2 in kg/m3
          X1 = (EE1*m1-m2)/FF
          Y1 = (m2-DD*m1)/FF

!! Calculate new concentrations for next step (continuation of VVWM method)
          rt1 = root1*T_end
          rt2 = root2*T_end
          exrt1 = exp(rt1)
          exrt2 = exp(rt2)
          ccc = X1*exrt1
          ddd = Y1*exrt2
          mn1 = ccc+ddd  !! End of timestep aqeous concentration in water column (kg/m3)
          mn2= DD*ccc+EE1*ddd !! End of timestep aqeous concentration in benthic (kg/m3)

!! Average daily (hourly) concentration calculations, assumes 1 day (hour) average (continuation of VVWM method) 
          gx=X1/root1
          hx=Y1/root2
          term1 = gx*exrt1                    
          term2 = hx*exrt2                    
          term3 = -gx
          term4 = -hx
          mavg1=(term1+term2+term3+term4)/T_end  !! Average of timestep aqueous concentration in water column (kg/m3)
          mavg2=(term1*DD+term2*EE1+term3*DD+term4*EE1)/T_end  !! Average of timestep aqueous concentration in benthic (kg/m3)

!! Convert end of time step concentrations back to total compatment mass, and from kg to mg
          !!Route out the water and sediment mass for end of timestep
          v1 = v1 - (rtwtr / 24.)
          mSed_1 = mSed_1 - ((sedrch / 24.) * 1000.)
          !!chpstmass = mn1 * 1000000. * ((v1 - (rtwtr / 24.)) / fw1)
          chpstmass = mn1 * 1000000. * (v1 / fw1) !! Water column mass is based on end of timestep volume (not subtracting routed water here)
          sedpstmass = mn2 * 1000000. * (v2 / fw2)  !!for benthic, end of timestep mass

!! verify that water concentration is at or below solubility, transfer mass to benthic if above solubility in water column
          solmax = 0.
          solmax = pest_sol * v1
          !!if (solmax < chpstmass * fw1) then
          !!   chpstmass = chpstmass - ((chpstmass * fw1) - solmax)
          !!   sedpstmass = sedpstmass + ((chpstmass * fw1) - solmax)
          !!end if

!! Calculate end of timestep total compartment concentration, based on end of timestep volumes
          if (v1 <= 0.) then
             chpst_conc(jrch) = 0. !! this would set concentration to 0 if there is no water left in channel after routing
          else
             chpst_conc(jrch) = chpstmass / v1
          end if
          sedpst_conc(jrch) = sedpstmass / bedvol

!! Set pesticide concentation transported out of reach based on timestep average concentration, converted from kg to mg, sum for all hours
!! MW Note, 4/28/20: Am thinking that we consider modifying SWAT outputs to print average in channel concentration (not the routed out concentration)

          solpesto = solpesto + (mavg1 * 1000000.)
          sorpesto = sorpesto + ((1. - fw1) * (mavg1 * 1000000. / fw1))

          End Do
          !! solpesto and sorpesto used in "reachout.f" for routing; gets multiplied by rtwtr to determine routed mass
          !! The average concentration in routed water (solpesto / rtwtr) should be equivalent to average concentration in reach
          solpesto = solpesto / 24.  !!Represents avgerage water column dissolved pesticide conc. over the day (mg/m3)
          sorpesto = sorpesto / 24.  !!Represents avgerage water column sorbed pesticide conc. over the day (mg/m3)
      else

!!Insignificant flow
          sedpstmass = sedpstmass + chpstmass
          chpstmass = 0.
          chpst_conc(jrch) = 0.
          sedpst_conc(jrch) = sedpstmass / bedvol
          solpesto = 0.
          sorpesto = 0.

      end if

      return
      end
