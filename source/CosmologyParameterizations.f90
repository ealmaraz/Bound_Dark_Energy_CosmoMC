    !Default parameterization using theta = r_s/D_a instead of H_0, and tau instead of z_re
    !and log(A_s) instead of A_s
    !Less general, but should give better performance
    !
    !The well-determined parameter A_s exp(-2tau) should be found by the covariance matrix
    !parameter 3 is 100*theta, parameter 4 is tau, others same as params_H except A->log(A)
    !Theta is much better constrained than H_0
    !
    !Also a background-only parameterization, e.g. for use with just supernoave etc

    module CosmologyParameterizations
    use CosmologyTypes
    use CosmoTheory
    use Calculator_Cosmology
    use bbn
    implicit none
    private

    !erick30.11.19 We introduce these variables to solve the BDE background
    real(mcp) cpot_param, alphaBDE_param, wDGc_param, Lambdac_param, ac_param, acLc_param
    real(mcp) grc_param, gDGc_param, OmegaDGc_param
    real(mcp) rhog0_param, rhon0_param, rhor0_param, rhom0_param, rhoDGc_param

    Type, extends(TCosmologyParameterization) :: ThetaParameterization
        real(mcp) :: H0_min = 40, H0_max = 100
        real(mcp) :: H0_prior_mean = 0._mcp, H0_prior_std = 0._mcp
        real(mcp) :: sterile_mphys_max = 10 !maximum allowed physical mass of thermal sterile neutrino in eV
        real(mcp) :: use_min_zre = 0._mcp
        real(mcp) :: zre_prior_mean = 0._mcp, zre_prior_std = 0._mcp
        integer :: num_derived = 0
    contains
    procedure :: ParamArrayToTheoryParams => TP_ParamArrayToTheoryParams
    procedure :: NonBaseParameterPriors => TP_NonBaseParameterPriors
    procedure :: CalcDerivedParams => TP_CalcDerivedParams
    procedure :: InitWithSetNames => TP_Init
    end type ThetaParameterization

    Type, extends(TCosmologyParameterization) :: BackgroundParameterization
    contains
    procedure :: ParamArrayToTheoryParams => BK_ParamArrayToTheoryParams
    procedure :: CalcDerivedParams => BK_CalcDerivedParams
    procedure :: InitWithSetNames => BK_Init
    end type BackgroundParameterization

    public BackgroundParameterization,ThetaParameterization

    contains


      !erick30.11.19 In this subroutine we solve the background equations (after ac) in BDE. When this
      !subroutine is called, ac, Lambdac, cpot, alpha, wDGc & OmegaDGc have been already initialized
      !in SetForH. The structure of this subroutine is very similar than 'bde_background' in
      !camb/equations_bde_v3.f90. However, here we don't export our results to an output file and
      !we don't use dverk to solve the differential equations.
      subroutine cosmomc_bde_background
         use Quint
         use constants

         !parameters for DLSODA
         integer neq, itol, itask, istate, iopt, lrw, liw, jt
         real(dl) y(3), atol(3), rwork(70), iwork(23)
         real(dl) t, tout, rtol
         real(dl) tc, t0

         integer i
         real(dl) splZero

         !input parameters for using DLSODA (see odepack.f for full documentation)
         neq = 3
         itol = 2
         rtol = 1.e-4
         atol(1) = 1.e-8
         atol(2) = 1.e-8
         atol(3) = 1.e-8
         itask = 1
         istate = 1
         iopt = 0
         lrw = 70
         liw = 23
         jt = 1

         !ac_param and others have been already initialized when calling SetForH
         tc = log(ac_param/ac_param)
         t0 = log(1._dl/ac_param)

         !initialization of dN
         dN = log(1._dl/ac_param)/(NumPoints-1)

         y(1) = sqrt(OmegaDGc_param*(1._mcp+wDGc_param)/2._mcp)
         y(2) = sqrt(OmegaDGc_param*(1._mcp-wDGc_param)/2._mcp)
         y(3) = sqrt(Mpl2)*alphaBDE_param/Lambdac_param

         !first point: initial conditions at ac
         Nend(1) = tc
         abde(1) = ac_param

         xbde(1) = y(1)
         ybde(1) = y(2)
         lbde(1) = y(3)

         phi_a(1) = alphaBDE_param/lbde(1)
         phidot_a(1) = abde(1)*sqrt(2._mcp*cosmomc_bde_v(phi_a(1)))*xbde(1)/ybde(1)

         !first call to DLSODA: this is just to break the initial stiffness of the problem
         t    = 0
         tout = 1.e-12
         CALL DLSODA(cosmomc_bde_equations,neq,y,t,tout,itol,rtol,atol,itask,istate,iopt,rwork,lrw,iwork,liw,cosmomc_jac,jt)
         if (istate /= 2) write(*,*) 'DLSODA failed!'

         !subsequent calls to DLSODA: we evolve the system from tout up to the present time
         do i=1, NumPointsEx-1
            tout = dN*i
            CALL DLSODA(cosmomc_bde_equations,neq,y,t,tout,itol,rtol,atol,itask,istate,iopt,rwork,lrw,iwork,liw,cosmomc_jac,jt)
            Nbde(i+1) = tout
            abde(i+1) = ac_param*exp(Nbde(i+1))

            xbde(i+1) = y(1)
            ybde(i+1) = y(2)
            lbde(i+1) = y(3)

            phi_a(i+1) = alphaBDE_param/lbde(i+1)
            phidot_a(i+1) = abde(i+1)*sqrt(2._mcp*cosmomc_bde_v(phi_a(i+1)))*xbde(i+1)/ybde(i+1)
         end do

         !Once the system is solved, we make some needed interpolations
         call spline(Nbde,xbde,NumPoints,splZero,splZero,ddxbde)
         call spline(Nbde,ybde,NumPoints,splZero,splZero,ddybde)
         call spline(Nbde,lbde,NumPoints,splZero,splZero,ddlbde)
         call spline(Nbde,phi_a,NumPointsEx,splZero,splZero,ddphi_a)
         call spline(Nbde,phidot_a,NumPointsEx,splZero,splZero,ddphidot_a)
       end subroutine cosmomc_bde_background


       !erick11.04.19. Scalar field potential NORMALIZED to Mpl^2:
       !                bde_v = V_std/Mpl^2 = cpot*Lambda^(4+alpha)/(Mpl2*phi^alpha)
       !This function has the same structure than 'Vofphi(phi,deriv=0)' in camb/equations_bde_v3.f90
       function cosmomc_bde_v(phi)
         use constants
         implicit none
         real(mcp) phi,cosmomc_bde_v
         real(mcp) norm

         norm = cpot_param*(Lambdac_param/sqrt(Mpl2))**(2._mcp+alphaBDE_param)*Lambdac_param**2*eV2toMpcm2
         cosmomc_bde_v = norm/(phi**alphaBDE_param)

       end  function cosmomc_bde_v


       !erick30.11.19. Dynamical system of background equations. The independent variable is t = ln(a/ac).
       !The dynamical variables are (see 'Dynamics of Dark Energy' by Copeland, Sami & Tsujikawa):
       !     x = phidot_std/(Mpl*\sqrt{6}*Hconf)  ,  y=a*sqrt(V_std)/(Mpl*\sqrt{3}*Hconf)
       !     l = -Mpl*V_std'/V_std                , Gamma = V_std*V_std''/V_std'**2
       !For the IPL potential V_std = cpot*Lambdac^(4+alpha)/phi_std^alpha & Gamma = 1 + 1/alpha.
       !This subroutine is essentially the same as 'bde_dlsoda_equations' in camb/equations_bde_v3.f90
       subroutine cosmomc_bde_equations(neq,t,y,ydot)
         double precision t, y, ydot
         dimension y(3), ydot(3)
         double precision ws
         integer neq

         ws = (1._mcp/3._mcp)*rhor0_param/(rhor0_param+ac_param*exp(t)*rhom0_param)

         ydot(1) = -3._mcp*y(1) + sqrt(3._mcp/2._mcp)*y(3)*y(2)**2+(3._mcp/2._mcp)*y(1)*(2._mcp*y(1)**2+(1._mcp+ws)*(1._mcp-y(1)**2-y(2)**2))
         ydot(2) = -sqrt(3._mcp/2._mcp)*y(1)*y(2)*y(3)+(3._mcp/2._mcp)*y(2)*(2._mcp*y(1)**2+(1._mcp+ws)*(1._mcp-y(1)**2-y(2)**2))
         ydot(3) = -sqrt(6._mcp)*y(1)*y(3)**2/alphaBDE_param
         return
       end subroutine cosmomc_bde_equations


       !erick30.11.19. Jacobian matrix df/dy for the dynamical system in subroutine 'cosmomc_bde_equations'.
       !               f is each one of ydot; y is x, y or l
       !               This subroutine is called by DLSODA
       !This subroutine is essentially the same as 'jac' in camb/equations_bde_v3.f90
       subroutine cosmomc_jac(neq,t,y,ml,mu,pd,nrowpd)
         integer neq, ml, mu, nrowpd
         double precision t, y, pd
         dimension y(neq), pd(nrowpd,neq)
         double precision ws

         ws = (1._mcp/3._mcp)*rhor0_param/(rhor0_param+ac_param*exp(t)*rhom0_param)

         pd(1,1) = -3._mcp+3._mcp/2._mcp*(6._mcp*y(1)*y(1)+(1._mcp+ws)*(1._mcp-3._mcp*y(1)*y(1)-y(2)*y(2)))
         pd(1,2) = sqrt(6._mcp)*y(3)*y(2)-3._mcp*y(1)*y(2)*(1._mcp+ws)
         pd(1,3) = sqrt(6._mcp)/2._mcp*y(2)*y(2)

         pd(2,1) = -sqrt(6._mcp)/2._mcp*y(3)*y(2)+3._mcp/2._mcp*y(2)*(4._mcp*y(1)-2._mcp*y(1)*(1._mcp+ws))
         pd(2,2) = -sqrt(6._mcp)/2._mcp*y(3)*y(1)+3._mcp/2._mcp*(2._mcp*y(1)*y(1)+(1._mcp+ws)*(1._mcp-y(1)*y(1)-3._mcp*y(2)*y(2)))
         pd(2,3) = -sqrt(6._mcp)/2._mcp*y(1)*y(2)

         pd(3,1) = -sqrt(6._mcp)/alphaBDE_param*y(3)*y(3)
         pd(3,2) = 0.
         pd(3,3) = -2._mcp*sqrt(6._mcp)/alphaBDE_param*y(1)*y(3)

         return
       end subroutine cosmomc_jac


    subroutine TP_Init(this, Ini, Names, Config)
    class(ThetaParameterization) :: this
    class(TSettingIni) :: Ini
    class(TParamNames) :: Names
    class(TGeneralConfig), target :: Config
    character(LEN=:), pointer :: prior

    call Ini%Read('H0_min',this%H0_min)
    call Ini%Read('H0_max',this%H0_max)
    call Ini%Read('use_min_zre',this%use_min_zre)
    call Ini%Read('sterile_mphys_max',this%sterile_mphys_max)
    prior => Ini%Read_String('H0_prior',NotFoundFail=.false.)
    if (prior/='') then
        read(prior,*) this%H0_prior_mean, this%H0_prior_std
    end if
    prior => Ini%Read_String('zre_prior',NotFoundFail=.false.)
    if (prior/='') then
        read(prior,*) this%zre_prior_mean, this%zre_prior_std
    end if

    !erick16.10.16. We use paramnames/params_CMB_BDE.paramnames for BDE
    call this%Initialize(Ini,Names, 'paramnames/params_CMB_BDE.paramnames', Config)

    if (CosmoSettings%bbn_consistency) call Names%Add('paramnames/derived_bbn.paramnames')
    call Names%Add('paramnames/derived_theory.paramnames')
    if (CosmoSettings%use_LSS) call Names%Add('paramnames/derived_LSS.paramnames')
    if (CosmoSettings%compute_tensors) call Names%Add('paramnames/derived_tensors.paramnames')

    !Add output ranges to match priors
    call Names%AddDerivedRange('zrei', mn=this%use_min_zre)
    call Names%AddDerivedRange('H0', this%H0_min, this%H0_max)
    this%num_derived = Names%num_derived

    !set number of hard parameters, number of initial power spectrum parameters
    !erick30.11.19 BDE has one more hard parameter than LCDM (to wit, cpot)
    !Therefore, we set 17 below
    call this%SetTheoryParameterNumbers(17,last_power_index)
    end subroutine TP_Init


    function TP_NonBaseParameterPriors(this,CMB)
    class(ThetaParameterization) :: this
    class(TTheoryParams) :: CMB
    real(mcp):: TP_NonBaseParameterPriors

    select type (CMB)
    class is (CMBParams)
        TP_NonBaseParameterPriors = logZero
        if (CMB%H0 < this%H0_min .or. CMB%H0 > this%H0_max) return
        if (CMB%zre < this%Use_min_zre) return
        if (CMB%omnuh2_sterile > 0 .and. CMB%nnu > standard_neutrino_neff) then
            !Check if physical mass of thermal massive sterile too big (look like CDM, so don't need to model separately)
            if (CMB%omnuh2_sterile*neutrino_mass_fac/(CMB%nnu-standard_neutrino_neff)**0.75_mcp > this%sterile_mphys_max) return
        end if
        TP_NonBaseParameterPriors = 0
        if (this%H0_prior_mean/=0._mcp) then
            TP_NonBaseParameterPriors = ((CMB%H0 - this%H0_prior_mean)/this%H0_prior_std)**2/2
        end if
        if (this%zre_prior_mean/=0._mcp) then
            TP_NonBaseParameterPriors = TP_NonBaseParameterPriors + ((CMB%zre - this%zre_prior_mean)/this%zre_prior_std)**2/2
        end if
    end select
    end function TP_NonBaseParameterPriors



    subroutine TP_ParamArrayToTheoryParams(this, Params, CMB)
    !erick30.11.19 camb/constants.f90 and some other stuff from the Quint module are needed here
    use constants
    use Quint
    use IFPORT
    class(ThetaParameterization) :: this
    real(mcp) Params(:)
    integer, parameter :: ncache =2
    Class(TTheoryParams), target :: CMB
    Type(CMBParams), save :: LastCMB(ncache)
    real(mcp) DA
    !real(mcp)  D_b,D_t,D_try,try_b,try_t, lasttry  !original: this variables are used in the binary search of H0

    !erick21.09.15 H0, wDE & theta_{MC} in BDE
    real(mcp) H0_BDE, wDE0_BDE, theta_BDE
    real(mcp) x0, y0, omde0, omrh2, ommh2

    !erick21.01.18. We introduce these variables to interpolate x, y and l at z_outputs.
    real(mcp) aaux, xaux,yaux,laux
    integer idx

    integer, save :: cache=1
    integer i
    Type(CMBParams), pointer :: CP2
    integer error

    select type(CosmoCalc=>this%Config%Calculator)
    class is (TCosmologyCalculator)
        select type (CMB)
        class is (CMBParams)
            do i=1, ncache
                !want to save two slow positions for some fast-slow methods
                if (all(Params(1:num_hard) == LastCMB(i)%BaseParams(1:num_hard))) then
                    CP2 => CMB !needed to make next line work for some odd reason CMB=LastCMB(i) does not work
                    CP2 = LastCMB(i)
                    call this%TCosmologyParameterization%ParamArrayToTheoryParams(Params, CMB)
                    call SetFast(Params,CMB)
                    return
                end if
            end do
            call this%TCosmologyParameterization%ParamArrayToTheoryParams(Params, CMB)
            error = 0   !JD to prevent stops when using bbn_consistency or m_sterile

            !!!!!!!!!!!!! erick30.11.19. Here we get H0, wx0 & theta by solving the background equations  !!!!!!!!!!!!!!!!
            !The procedure is as follows:
            !i) We firstly call SetForH just to initialize the hard parameters. At this point
            !H0 = H0_max, wx0=-1 & 100*theta=1.0411. This is just a harmless initialization used to
            !keep the program running
            call SetForH(Params,CMB,this%H0_max,-1._mcp,1.0411_mcp,.true.,error)

            !ii)Once alphaBDE & others have been initialized, we solve the background
            !Note also that the temporarily values of H0, wx0 & theta are irrelevant
            call cosmomc_bde_background

            !Once xbde & ybde are filled we retrieve the third last elements corresponding to the present epoch
            x0 = xbde(NumPoints)
            y0 = ybde(NumPoints)

            !iii) Here we compute the correct value of wx0 & H0. In this respect, see the code in 'bde_h0'
            !in camb/equations_bde_v3.f90
            wDE0_BDE = (x0**2-y0**2)/(x0**2+y0**2)
            omde0 = x0**2+y0**2
            omrh2 = c8piG_eVm2*CMB%rhor0*eV2toMpcm2/3._mcp
            ommh2 = 1.e10*(CMB%ombh2+CMB%omch2)/(c**2)
            H0_BDE = c/(1.e3)*sqrt((omrh2+ommh2)/(1-omde0))

            !iv)We now call SetForH for the second time to update the info of wx0, H0 and
            !the omegas. With this info we compute the correct value of theta in BDE
            call SetForH(Params,CMB,H0_BDE,wDE0_BDE,1.0411_mcp,.false.,error)
            theta_BDE = 100._mcp*CosmoCalc%CMBToTheta(CMB)

            !v) Finally we call SetForH once again to update the value of theta. After this,
            !everything is fine and the program continues running
            call SetForH(Params,CMB,H0_BDE,wDE0_BDE,theta_BDE,.false.,error)
            !write(*,*) 'wDE0:', wDE0_BDE, 'H0:', H0_BDE, 'theta:', theta_BDE

            !erick30.11.19 We fill in the arrays CMB%xdyn & CMB%ydyn
            do idx=1, size(CosmoSettings%z_outputs)
               aaux=1./(1.+CosmoSettings%z_outputs(idx))
               !cosmomc_bde_xylAta is implemented in camb/equations_bde_v3
               call cosmomc_bde_xylAta(aaux,xaux,yaux,laux)
               CMB%xdyn(idx) = xaux
               CMB%ydyn(idx) = yaux
            end do

            !!!!!!!!!!! erick - comment. Calculation of H0 from the value of \theta_{*} in LCDM !!!!!!!!!!!!!!!!
            !DA = Params(3)/100
            !try_b = this%H0_min
            !!erick - comment. H0 is initialized after CosmologyParametrizations/SetForH is called for the first time
            !call SetForH(Params,CMB,try_b, .true.,error)  !JD for bbn related errors
            !if(error/=0)then
            !    cmb%H0=0
            !    return
            !end if
            !D_b = CosmoCalc%CMBToTheta(CMB)

            !try_t = this%H0_max
            !call SetForH(Params,CMB,try_t, .false.)
            !D_t = CosmoCalc%CMBToTheta(CMB)

            !erick - comment - Info displaying
            !write(*,321) 'H0min', '=', this%H0_min, 'H0max' ,'=',this%H0_max
            !321 FORMAT ('',A5,T7,A1,ES12.4,T24,A5,T31,A1,ES12.4)
            !write(*,654) 't_inf', '=', D_b, 'Pa(3)', '=',DA, 't_sup', '=', D_t_test
            !654 FORMAT('',A5,T7,A1,ES12.4,T24,A5,T31,A1,ES12.4,T49,A5,T56,A1,ES12.4)

            !if (DA < D_b .or. DA > D_t) then
            !    if (Feedback>1) write(*,*) instance, 'Out of range finding H0: ', real(Params(3))
            !    cmb%H0=0 !Reject it
            !else
            !    lasttry = -1
            !    !erick - comment. Binary search of H0 consistent with the input value of 100\Theta_{*}
            !    do
            !        call SetForH(Params,CMB,(try_b+try_t)/2, .false.)
            !        D_try = CosmoCalc%CMBToTheta(CMB)
            !        if (D_try < DA) then
            !            try_b = (try_b+try_t)/2
            !        else
            !            try_t = (try_b+try_t)/2
            !        end if
            !        if (abs(D_try - lasttry)< 1e-7) exit
            !        lasttry = D_try
            !    end do

                !erick - comment. After the search is done, the program continues here
                !!call InitCAMB(CMB,error)
            !    if (CMB%tau==0._mcp) then
            !        CMB%zre=0
            !    else
            !        CMB%zre = CosmoCalc%GetZreFromTau(CMB, CMB%tau)
            !    end if

            !    LastCMB(cache) = CMB
            !    cache = mod(cache,ncache)+1
            !end if
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

            !erick30.11.19 These lines are the same ones written above
            if (CMB%tau==0._mcp) then
                CMB%zre=0
            else
                CMB%zre = CosmoCalc%GetZreFromTau(CMB, CMB%tau)
            end if

            LastCMB(cache) = CMB
            cache = mod(cache,ncache)+1

            !erick30.11.19 Here we omit 'end if' since we are not within the branch
            !                  if (DA < D_b .or. DA > D_t) then
        end select
        class default
        call MpiStop('CosmologyParameterizations: Calculator is not TCosmologyCalculator')
    end select

    end subroutine TP_ParamArrayToTheoryParams


    function GetYPBBN(Yhe)
    !Convert yhe defined as mass fraction (CMB codes), to nucleon ratio definition
    real(mcp), intent(in) :: Yhe
    real(mcp) GetYPBBN
    real(mcp), parameter :: m_proton = 1.672621637e-27
    real(mcp), parameter :: m_H = 1.673575e-27
    real(mcp), parameter :: not4 = 3.9715
    real(mcp), parameter :: m_He = m_H * not4

    GetYPBBN =  4 * m_H * Yhe / (m_He - Yhe * (m_He - 4*m_H))

    end function GetYPBBN


    subroutine TP_CalcDerivedParams(this, P, Theory, derived)
    class(ThetaParameterization) :: this
    real(mcp), allocatable :: derived(:)
    class(TTheoryPredictions), allocatable :: Theory
    real(mcp) :: P(:)
    Type(CMBParams) CMB
    real(mcp) :: lograt
    integer ix,i
    real(mcp) z
    !erick16.10.16 Extra contribution to N_eff by the Dark Group in the early universe
    real(mcp) NextBDE
    integer, parameter :: derivedCL(5) = [40, 220, 810, 1420, 2000]

    !erick30.11.19 Used to fill in the arrays corresponding to xbde & ybde
    integer idx

    !erick - comment
    !write(*,*) 'TP_CalcDerivedParams has been called; a new sample is added'

    if (.not. allocated(Theory)) call MpiStop('Not allocated theory!!!')
    select type (Theory)
    class is (TCosmoTheoryPredictions)
        allocate(Derived(this%num_derived), source=0._mcp)

        call this%ParamArrayToTheoryParams(P,CMB)

        derived(1) = CMB%H0
        derived(2) = CMB%omv
        derived(3) = CMB%omdm+CMB%omb
        derived(4) = CMB%omdmh2 + CMB%ombh2
        derived(5) = CMB%omnuh2
        derived(6) = (CMB%omdmh2 + CMB%ombh2)*CMB%h

        derived(7) = Theory%Sigma_8
        derived(8) = Theory%Sigma_8*((CMB%omdm+CMB%omb))**0.5_mcp
        derived(9) = Theory%Sigma_8*((CMB%omdm+CMB%omb))**0.25_mcp
        derived(10)= Theory%Sigma_8/CMB%h**0.5_mcp

        derived(11) = Theory%Lensing_rms_deflect
        derived(12) = CMB%zre

        !erick30.11.19. We export the info of ac, OmegaDGc, wx0 and theta
        derived(13) = CMB%ac
        derived(14) = CMB%OmegaDGc
        derived(15) = CMB%wx0
        derived(16) = CMB%theta

        !erick30.11.19 We export the info of x & y at z_outputs
        do idx = 1, size(CosmoSettings%z_outputs)
           derived(16+idx)=CMB%xdyn(idx)
           derived(16+size(CosmoSettings%z_outputs)+idx)=CMB%ydyn(idx)
        end do

        !erick30.11.19 Since size(z_outputs)=37 and there is x & y for each z,
        !we set ix from 17 to 91
        ix=91
        derived(ix) = cl_norm*CMB%InitPower(As_index)*1e9
        derived(ix+1) = derived(ix)*exp(-2*CMB%tau)  !A e^{-2 tau}
        ix = ix+2

        if(CosmoSettings%use_CMB .and. allocated(Theory%Cls(1,1)%CL)) then
            !L(L+1)C_L/2pi at various places
            derived(ix:ix+size(DerivedCL)-1) = Theory%Cls(1,1)%CL(derivedCL)
        end if
        ix = ix+size(derivedCL)

        lograt = log(0.002_mcp/CosmoSettings%pivot_k)   !get ns at k=0.002
        derived(ix) = CMB%InitPower(ns_index) +CMB%InitPower(nrun_index)*lograt +&
            CMB%InitPower(nrunrun_index)*lograt**2/2
        ix=ix+1

        derived(ix)= CMB%Yhe !value actually used, may be set from bbn consistency
        derived(ix+1)= GetYpBBN(CMB%Yhe) !same, as nucleon ratio definition
        ix = ix+2

        !erick30.11.19 We compute the fraction of primordial deuterium using NextBDE
        NextBDE = CMB%Nex
        if (CosmoSettings%bbn_consistency) then
            !derived(ix) = 1d5*BBN_DH%Value(CMB%ombh2,CMB%nnu - standard_neutrino_neff) !original
            derived(ix) = 1d5*BBN_DH%Value(CMB%ombh2,NextBDE)
            ix =ix + 1
        end if

        derived(ix:ix + Theory%numderived-1) = Theory%derived_parameters(1: Theory%numderived)
        ix = ix + Theory%numderived

        if (CosmoSettings%Use_LSS) then
            ! f sigma_8 at specified redshift
            do i=1,size(CosmoSettings%z_outputs)
                z =  CosmoSettings%z_outputs(i)
                !erick - comment. Here we can verify which z's we are analyzing
                !write(*,*) 'fsigma8 & sigma8 at z = ', z
                derived(ix) = Theory%growth_z%Value(z)
                derived(ix+1) = Theory%sigma8_z%Value(z)
                ix = ix + 2
            end do
        end if

        if (CosmoSettings%Compute_tensors) then
            derived(ix:ix+5) = [Theory%tensor_ratio_02, Theory%tensor_ratio_BB, log(Theory%tensor_AT*1e10), &
                Theory%tensor_ratio_C10, Theory%tensor_AT*1e9, Theory%tensor_AT*1e9*exp(-2*CMB%tau) ]
            ix=ix+6
        end if

        if (ix - 1 /= this%num_derived) then
            write(*,*) 'num_derived =', this%num_derived, '; ix, Theory%numderived = ', ix, Theory%numderived
            call MpiStop('TP_CalcDerivedParams error in derived parameter numbers')
        end if
    end select

    end subroutine TP_CalcDerivedParams


    subroutine SetFast(Params,CMB)
    real(mcp) Params(num_Params)
    Type(CMBParams) CMB

    CMB%InitPower(1:num_initpower) = Params(index_initpower:index_initpower+num_initpower-1)
    CMB%InitPower(As_index) = exp(CMB%InitPower(As_index))

    end subroutine SetFast


    !erick21.09.15. We extend the number of input parameters of this subroutine.
    !This is to compute wx0, H0 & theta in our model (see TP_ParamArrayToTheoryParams above)
    subroutine SetForH(Params,CMB,H0,wx0,thetaq,firsttime,error)
    !erick30.11.15 constants and settings are required
    use constants
    use settings
    use bbn

    real(mcp) Params(num_Params)
    logical, intent(in) :: firsttime
    Type(CMBParams) CMB
    real(mcp) h2,H0

    !erick16.10.16
    real(mcp) wx0, thetaq, NextBDE

    integer, optional :: error

    !erick - comment. In LCDM CMB%H0 is assigned at the beginning
    !CMB%H0=H0

    if (firsttime) then
      !erick30.11.19 Initialization of the hard parameters and other calculations are done here

      !rho of radiation today (in eV**4) and acLc
      rhog0_param = const_pi**2/15._mcp*COBE_CMBTemp**4*(k_B/eV)**4
      rhon0_param = 7._mcp/8._mcp*(4._mcp/11._mcp)**(4._mcp/3._mcp)*default_nnu*rhog0_param
      rhor0_param = rhog0_param+rhon0_param
      CMB%rhor0 = rhor0_param
      !degrees of freedom of photons + neutrinos after nu_dec
      !see above eq 12 in PRD 99, 103504 (2019)
      grc_param = 2._mcp + 7._mcp/4._mcp*default_nnu*(4._mcp/11._mcp)**(4._mcp/3._mcp)

      !erick30.11.19. Initialization of the BDE parameters
      !wDGc, Ncolors & Nflavors are read in driver.F90.
      !At this point, wDGc is not expressed in multiples of 1/3 anymore!
      CMB%wDGc = wDGc
      wDGc_param = CMB%wDGc

      !see eq 18 in PRD 99, 103504 (2019)
      CMB%alphaBDE = 2._mcp*(1._mcp+2._mcp/(Ncolors-Nflavors))
      alphaBDE_param = CMB%alphaBDE

      CMB%cpot = Params(10)
      cpot_param = CMB%cpot

      CMB%Lambdac = Params(9)
      Lambdac_param = CMB%Lambdac

      !degrees of freedom of the Dark Group (eq 9 in PRD 99, 103504 (2019))
      gDGc_param = (1._mcp+7._mcp/8._mcp)*(2._mcp*(Ncolors**2-1._mcp)+2._mcp*Ncolors*Nflavors)
      !ac*Lambdac (eq 15 in PRD 99, 103504 (2019))
      acLc_param = ((1._mcp-CMB%wDGc)*rhor0_param*gDGc_param/(2._mcp*CMB%cpot*grc_param))**(1._mcp/4._mcp)*(4._mcp*gsmM/(11._mcp*gsmgut))**(1._mcp/3._mcp)
      CMB%ac = acLc_param/CMB%Lambdac
      ac_param = CMB%ac

      !extra contribution to N_eff by the Dark Group in the early universe
      !see section II. THE BOUND DARK ENERGY MODEL in PRD 99, 103504 (2019)
      CMB%Nex = 2._mcp/(1._mcp-CMB%wDGc)*CMB%cpot*acLc_param**4*8._mcp/7._mcp*(11._mcp/4._mcp)**(4._mcp/3._mcp)/rhog0_param

      !OmegaDGc
      rhoDGc_param= 2._mcp/(1._mcp-wDGc_param)*cpot_param*Lambdac_param**4
      rhom0_param = 3._mcp/c8piG_eVm2*1.e10/(c**2)*(Params(1)+Params(2))/eV2toMpcm2
      CMB%OmegaDGc = rhoDGc_param/(rhom0_param/(CMB%ac**3)+rhor0_param/(CMB%ac**4)+rhoDGc_param)
      OmegaDGc_param = CMB%OmegaDGc

        !erick30.11.19 Initialization of the standard cosmological hard parameters
        !see paramnames/params_CMB_BDE.paramnames
        CMB%reserved = 0
        CMB%ombh2 = Params(1)
        CMB%tau = params(3) !tau, set zre later
        CMB%Omk = Params(4)
        CMB%w = Params(7)
        CMB%wa = Params(8)
        CMB%nnu = Params(11) !3.046
        !Params(5) is now mnu, where mnu is physical standard neutrino mass and we assume standard heating
        CMB%sum_mnu_standard = Params(5)
        CMB%omnuh2 = Params(5)/neutrino_mass_fac*(standard_neutrino_neff/3)**0.75_mcp
        !Params(6) is mass_sterile*Neff_sterile
        CMB%omnuh2_sterile = Params(6)/neutrino_mass_fac
        !we are using interpretation where there are degeneracy_factor neutrinos, each exactly thermal
        !So internally 3.046 or 3.046/3 massive neutrnos. But mnu is the physical integer mass sum.
        if (CMB%omnuh2_sterile >0 .and. CMB%nnu < standard_neutrino_neff) then
            if(present(error))then
                error=-1
            else
                call MpiStop('sterile neutrino mass required Neff>3.046')
            end if
        end if

        CMB%omnuh2 = CMB%omnuh2 + CMB%omnuh2_sterile
        CMB%omch2 = Params(2)
        !erick30.11.19 We consider the simplest scenario where neutrinos remain massless.
        !Therefore the total dark matter component is only due to CDM, which implies that CMB%nufrac = 0
        !CMB%omdmh2 = CMB%omch2+ CMB%omnuh2   !original
        CMB%omdmh2 = CMB%omch2
        CMB%nufrac = CMB%omnuh2/CMB%omdmh2

        !erick30.11.19 We compute the fraction of primordial helium using NextBDE
        NextBDE = CMB%Nex
        if (CosmoSettings%bbn_consistency) then
            !CMB%YHe = BBN_YHe%Value(CMB%ombh2,CMB%nnu - standard_neutrino_neff,error) !original
            CMB%YHe = BBN_YHe%Value(CMB%ombh2,NextBDE,error)
        else
            !e.g. set from free parameter..
            CMB%YHe  =Params(12)
        end if

        CMB%iso_cdm_correlated =  Params(13)
        CMB%zre_delta = Params(14)
        CMB%ALens = Params(15)
        CMB%ALensf = Params(16)
        CMB%fdm = Params(17)
        call SetFast(Params,CMB)
    end if

    !erick30.11.19 Here we compute omb, omc and update H0, wx0 & theta 
    CMB%H0 = H0
    CMB%h = CMB%H0/100
    h2 = CMB%h**2
    CMB%wx0 = wx0
    CMB%theta = thetaq
    CMB%omb = CMB%ombh2/h2
    CMB%omc = CMB%omch2/h2
    CMB%omnu = CMB%omnuh2/h2
    CMB%omdm = CMB%omdmh2/h2
    CMB%omv = 1- CMB%omk - CMB%omb - CMB%omdm
    end subroutine SetForH


    !!! Simple parameterization for background data, e.g. Supernovae only (no thermal history)
    subroutine BK_Init(this, Ini, Names, Config)
    class(BackgroundParameterization) :: this
    class(TSettingIni) :: Ini
    class(TParamNames) :: Names
    class(TGeneralConfig), target :: Config

    this%late_time_only = .true.
    call this%Initialize(Ini,Names, 'paramnames/params_background.paramnames', Config)
    call this%SetTheoryParameterNumbers(Names%num_MCMC,0)

    end subroutine BK_Init

    subroutine BK_ParamArrayToTheoryParams(this, Params, CMB)
    class(BackgroundParameterization) :: this
    real(mcp) Params(:)
    class(TTheoryParams), target :: CMB
    real(mcp) omegam, h2

    select type (CMB)
    class is (CMBParams)
        omegam = Params(1)
        CMB%H0 = Params(2)
        CMB%omk = Params(3)
        CMB%omnuh2=Params(4)/neutrino_mass_fac*(standard_neutrino_neff/3)**0.75_mcp
        CMB%w =    Params(5)
        CMB%wa =    Params(6)
        CMB%nnu =    Params(7)

        CMB%h=CMB%H0/100
        h2 = CMB%h**2
        CMB%Yhe=0.24
        CMB%omnu = CMB%omnuh2/h2
        CMB%omb= omegam - CMB%omnu
        CMB%ombh2 = CMB%omb*h2
        CMB%omc=0
        CMB%omch2 = CMB%omc*h2
        CMB%zre=0
        CMB%tau=0
        CMB%omdmh2 = CMB%omch2+ CMB%omnuh2
        CMB%omdm = CMB%omdmh2/h2
        CMB%omv = 1- CMB%omk - CMB%omb - CMB%omdm
        CMB%nufrac=CMB%omnuh2/CMB%omdmh2
        CMB%reserved=0
        CMB%fdm=0
        CMB%iso_cdm_correlated=0
        CMB%Alens=1
    end select
    end subroutine BK_ParamArrayToTheoryParams


    subroutine BK_CalcDerivedParams(this, P, Theory, derived)
    class(BackgroundParameterization) :: this
    real(mcp), allocatable :: derived(:)
    class(TTheoryPredictions), allocatable :: Theory
    real(mcp) :: P(:)
    Type(CMBParams) CMB

    allocate(Derived(1))

    call this%ParamArrayToTheoryParams(P,CMB)

    derived(1) = CMB%omv

    end subroutine BK_CalcDerivedParams


    end module CosmologyParameterizations
