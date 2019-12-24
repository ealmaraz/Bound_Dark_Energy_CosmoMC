    !     Code for Anisotropies in the Microwave Background
    !     by Antony Lewis (http://cosmologist.info/) and Anthony Challinor
    !     See readme.html for documentation. This is a sample driver routine that reads
    !     in one set of parameters and produdes the corresponding output.

    program driver
    use IniFile
    use CAMB
    use LambdaGeneral
    use Lensing
    use AMLUtils
    use Transfer
    use constants
    use Bispectrum
    use CAMBmain
    use NonLinear
    !erick11.04.19 To connect with <params>_quint in equations_bde_v3 (see below)
    use Quint
#ifdef NAGF95
    use F90_UNIX
#endif
    implicit none

    Type(CAMBparams) P

    character(LEN=Ini_max_string_len) numstr, VectorFileName, &
        InputFile, ScalarFileName, TensorFileName, TotalFileName, LensedFileName,&
        LensedTotFileName, LensPotentialFileName,ScalarCovFileName

    !erick11.04.19. Name of the output files for SneIa, BAO and Press-Schechter
    character(len=:), allocatable :: dL_filename, rBAO_filename, PressSchechter_filename
    integer len_dL_filename, len_rBAO_filename, len_PressSchechter_filename

    !erick16.08.17. SNeIa, BAO & Press-Schechter info required for producing the output files
    real(dl) zSNeIa, zSNeIa_min, zSNeIa_max, dL_SNeIa, mu_SNeIa
    integer i_zSNeIa, n_zSNeIa

    real(dl) zBAO, zBAO_min, zBAO_max, DV_BAO, r_BAO
    integer i_zBAO, n_zBAO

    integer i_PS, n_PS
    real(dl) kPS_min, kPS_max, kPS, RadiusPS
    real(dl) NPS_min, NPS_max, NPS
    real(dl) Msolar, MassPS_Mpc, MassPS_Solar
    real(dl) sigmaPS(1)

    !erick11.04.19. rhos in eV**4 today
    real(dl) rhog0, rhon0, rhor0, rhom0, rhol, rhocrit0

    !erick11.04.19. grc : degrees of freedom of photons + neutrinos after nu_dec
    !               gDGc: degrees of freedom of the Dark Group
    !               omega_radiation & omega_de are used to compute Lambdac for use_physical = F
    real(dl) grc, gDGc, omega_radiation, omega_de

    integer i
    character(LEN=Ini_max_string_len) TransferFileNames(max_transfer_redshifts), &
        MatterPowerFileNames(max_transfer_redshifts), outroot, version_check
    real(dl) output_factor, nmassive

#ifdef WRITE_FITS
    character(LEN=Ini_max_string_len) FITSfilename
#endif

    logical bad

    InputFile = ''
    if (GetParamCount() /= 0)  InputFile = GetParam(1)
    if (InputFile == '') error stop 'No parameter input file'

    call Ini_Open(InputFile, 1, bad, .false.)
    if (bad) error stop 'Error opening parameter file'

    Ini_fail_on_not_found = .false.

    outroot = Ini_Read_String('output_root')
    if (outroot /= '') outroot = trim(outroot) // '_'

    highL_unlensed_cl_template = Ini_Read_String_Default('highL_unlensed_cl_template',highL_unlensed_cl_template)

    call CAMB_SetDefParams(P)

    P%WantScalars = Ini_Read_Logical('get_scalar_cls')
    P%WantVectors = Ini_Read_Logical('get_vector_cls',.false.)
    P%WantTensors = Ini_Read_Logical('get_tensor_cls',.false.)

    P%OutputNormalization=outNone
    output_factor = Ini_Read_Double('CMB_outputscale',1.d0)

    P%WantCls= P%WantScalars .or. P%WantTensors .or. P%WantVectors

    P%PK_WantTransfer=Ini_Read_Logical('get_transfer')

    AccuracyBoost  = Ini_Read_Double('accuracy_boost',AccuracyBoost)
    lAccuracyBoost = Ini_Read_Real('l_accuracy_boost',lAccuracyBoost)
    HighAccuracyDefault = Ini_Read_Logical('high_accuracy_default',HighAccuracyDefault)

    P%NonLinear = Ini_Read_Int('do_nonlinear',NonLinear_none)

    P%DoLensing = .false.
    if (P%WantCls) then
        if (P%WantScalars  .or. P%WantVectors) then
            P%Max_l = Ini_Read_Int('l_max_scalar')
            P%Max_eta_k = Ini_Read_Double('k_eta_max_scalar',P%Max_l*2._dl)
            if (P%WantScalars) then
                P%DoLensing = Ini_Read_Logical('do_lensing',.false.)
                if (P%DoLensing) lensing_method = Ini_Read_Int('lensing_method',1)
            end if
            if (P%WantVectors) then
                if (P%WantScalars .or. P%WantTensors) error stop 'Must generate vector modes on their own'
                i = Ini_Read_Int('vector_mode')
                if (i==0) then
                    vec_sig0 = 1
                    Magnetic = 0
                else if (i==1) then
                    Magnetic = -1
                    vec_sig0 = 0
                else
                    error stop 'vector_mode must be 0 (regular) or 1 (magnetic)'
                end if
            end if
        end if

        if (P%WantTensors) then
            P%Max_l_tensor = Ini_Read_Int('l_max_tensor')
            P%Max_eta_k_tensor =  Ini_Read_Double('k_eta_max_tensor',Max(500._dl,P%Max_l_tensor*2._dl))
        end if
    endif

    !  Read initial parameters.

    !erick15.01.18 We import the mode to be evolved. Zero means all of the modes to be evolved
    !erick30.11.19 CosmoMC modification. In CosmoMC we always evolve all of the modes, so
    !              we don't need to read this variable from params.ini
    !fixq = Ini_Read_Double('fixq',0._dl)

    !erick16.09.19. Initialization of 'quint' parameters of module Quint. These parameters are used
    !               here and in equations_bde_v3.f90
    cpot_quint = Ini_Read_Double('cpot',1._dl)
    P%cpot = cpot_quint

    wDGc_quint = Ini_Read_Double('wDGc',1._dl)/3._dl
    P%wDGc = wDGc_quint

    alphaBDE_quint = 2._dl*(1._dl+2._dl/(Ini_Read_Double('Ncolors',3._dl)-Ini_Read_Double('Nflavors',6._dl)))
    P%alphaBDE = alphaBDE_quint

    !erick09.04.19. rho of radiation (in eV**4) and acLc
    rhog0_quint = const_pi**2/15._dl*COBE_CMBTemp**4*(k_B/eV)**4
    rhon0_quint = 7._dl/8._dl*(4._dl/11._dl)**(4._dl/3._dl)*default_nnu*rhog0_quint
    rhor0_quint = rhog0_quint+rhon0_quint
    P%rhor0 = rhor0_quint
    rhor0 = rhor0_quint

    !degrees of freedom of photons + neutrinos after nu_dec (see above eq 12 in PRD 99, 103504 (2019))
    grc = 2._dl + 7._dl/4._dl*default_nnu*(4._dl/11._dl)**(4._dl/3._dl)
    !degrees of freedom of the Dark Group (eq 9 in PRD 99, 103504 (2019))
    gDGc = (1._dl+7._dl/8._dl)*(2*(Ini_Read_Double('Ncolors')**2-1._dl)+2._dl*Ini_Read_Double('Ncolors')*Ini_Read_Double('Nflavors'))
    !ac*Lambdac (eq 15 in PRD 99, 103504 (2019))
    acLc_quint = ((1._dl-wDGc_quint)*rhor0_quint*gDGc/(2._dl*cpot_quint*grc))**(1._dl/4._dl)*(4._dl*gsmM/(11._dl*gsmgut))**(1._dl/3._dl)

    !erick09.04.19 use_physical flag imported from params.ini
    uphys = Ini_Read_Logical('use_physical',.true.)

    if (Ini_Read_Logical('use_physical',.false.)) then
        rhom0_quint = 3._dl/c8piG_eVm2*1.e10/(c**2)*(Ini_Read_Double('ombh2')+Ini_Read_Double('omch2'))/eV2toMpcm2
	      !rhom0 is required in Press-Schechter
	      rhom0 = rhom0_quint
        Lambdac_quint = Ini_Read_Double('Lambdac')
        P%Lambdac = Lambdac_quint
        ac_quint = acLc_quint/Lambdac_quint
        P%ac = ac_quint
        !erick09.04.19 Here we determine H0 from Lambdac, omegabh2 & omegach2
	      call bde_h0(P,P%H0)
        P%omegab = Ini_Read_Double('ombh2')/(P%H0/100)**2
        P%omegac = Ini_Read_Double('omch2')/(P%H0/100)**2
        !erick02.01.16. To neglect massive neutrinos check the modifications in params.ini
        P%omegan = Ini_Read_Double('omnuh2')/(P%H0/100)**2
        P%omegav = 1-Ini_Read_Double('omk')-P%omegab-P%omegac
    else
        !erick11.04.19. We now support use_physical = F. In this case H0 is simply read from params.ini
	      P%H0 = Ini_Read_Double('hubble')
        P%omegab = Ini_Read_Double('omega_baryon')
        P%omegac = Ini_Read_Double('omega_cdm')
        P%omegan = Ini_Read_Double('omega_neutrino')
        P%omegav = 1-Ini_Read_Double('omk')-P%omegab-P%omegac
        !critical density today (in eV**4)
        rhocrit0_quint = (3.*(P%H0/100.)**2/c8piG_eVm2)*(1.e10/(c**2*eV2toMpcm2))
        rhom0_quint = 3./c8piG_eVm2*1.e10/(c**2)*(P%omegab*(P%H0/100.)**2+P%omegac*(P%H0/100.)**2)/eV2toMpcm2
        rhom0 = rhom0_quint
        rhol_quint = rhocrit0_quint-(rhor0_quint+rhom0_quint)

        !erick11.04.19. In this case Lambdac is determined by a binary search in equations_bde_v3
        omega_radiation = rhor0_quint/rhocrit0_quint
	      omega_de = 1._dl-P%omegab-P%omegac-omega_radiation
        Lambdac_quint = bde_search_lambdac(omega_de)
        P%Lambdac = Lambdac_quint
        ac_quint = acLc_quint/Lambdac_quint
        P%ac = ac_quint
	  end if

    P%tcmb   = Ini_Read_Double('temp_cmb',COBE_CMBTemp)
    P%yhe    = Ini_Read_Double('helium_fraction',0.24_dl)

    !erick - comment. Check the modificacions in params.ini, where we set the things properly to consider
    !                 massless neutrinos
    P%Num_Nu_massless  = Ini_Read_Double('massless_neutrinos')
    P%Nu_mass_eigenstates = Ini_Read_Int('nu_mass_eigenstates',1)
    if (P%Nu_mass_eigenstates > max_nu) error stop 'too many mass eigenstates'

    numstr = Ini_Read_String('massive_neutrinos')
    read(numstr, *) nmassive
    if (abs(nmassive-nint(nmassive))>1e-6) error stop 'massive_neutrinos should now be integer (or integer array)'
    read(numstr,*, end=100, err=100) P%Nu_Mass_numbers(1:P%Nu_mass_eigenstates)
    P%Num_Nu_massive = sum(P%Nu_Mass_numbers(1:P%Nu_mass_eigenstates))

    if (P%Num_Nu_massive>0) then
        P%share_delta_neff = Ini_Read_Logical('share_delta_neff', .true.)
        numstr = Ini_Read_String('nu_mass_degeneracies')
        if (P%share_delta_neff) then
            if (numstr/='') write (*,*) 'WARNING: nu_mass_degeneracies ignored when share_delta_neff'
        else
            if (numstr=='') error stop 'must give degeneracies for each eigenstate if share_delta_neff=F'
            read(numstr,*) P%Nu_mass_degeneracies(1:P%Nu_mass_eigenstates)
        end if
        numstr = Ini_Read_String('nu_mass_fractions')
        if (numstr=='') then
            if (P%Nu_mass_eigenstates >1) error stop 'must give nu_mass_fractions for the eigenstates'
            P%Nu_mass_fractions(1)=1
        else
            read(numstr,*) P%Nu_mass_fractions(1:P%Nu_mass_eigenstates)
        end if
    end if

    !JD 08/13 begin changes for nonlinear lensing of CMB + LSS compatibility
    !P%Transfer%redshifts -> P%Transfer%PK_redshifts and P%Transfer%num_redshifts -> P%Transfer%PK_num_redshifts
    !in the P%WantTransfer loop.
    if (((P%NonLinear==NonLinear_lens .or. P%NonLinear==NonLinear_both) .and. P%DoLensing) &
        .or. P%PK_WantTransfer) then
    P%Transfer%high_precision=  Ini_Read_Logical('transfer_high_precision',.false.)
    else
        P%transfer%high_precision = .false.
    endif
    if (P%NonLinear/=NonLinear_none) call NonLinear_ReadParams(DefIni)

    if (P%PK_WantTransfer)  then
        P%WantTransfer  = .true.
        P%transfer%kmax          =  Ini_Read_Double('transfer_kmax')
        P%transfer%k_per_logint  =  Ini_Read_Int('transfer_k_per_logint')
        P%transfer%PK_num_redshifts =  Ini_Read_Int('transfer_num_redshifts')

        transfer_interp_matterpower = Ini_Read_Logical('transfer_interp_matterpower ', transfer_interp_matterpower)
        transfer_power_var = Ini_read_int('transfer_power_var',transfer_power_var)
        if (P%transfer%PK_num_redshifts > max_transfer_redshifts) error stop 'Too many redshifts'
        do i=1, P%transfer%PK_num_redshifts
            P%transfer%PK_redshifts(i)  = Ini_Read_Double_Array('transfer_redshift',i,0._dl)
            transferFileNames(i)     = Ini_Read_String_Array('transfer_filename',i)
            MatterPowerFilenames(i)  = Ini_Read_String_Array('transfer_matterpower',i)
            if (TransferFileNames(i) == '') then
                TransferFileNames(i) =  trim(numcat('transfer_',i))//'.dat'
            end if
            if (MatterPowerFilenames(i) == '') then
                MatterPowerFilenames(i) =  trim(numcat('matterpower_',i))//'.dat'
            end if
            if (TransferFileNames(i)/= '') &
                TransferFileNames(i) = trim(outroot)//TransferFileNames(i)
            if (MatterPowerFilenames(i) /= '') &
                MatterPowerFilenames(i)=trim(outroot)//MatterPowerFilenames(i)
        end do
    else
        P%Transfer%PK_num_redshifts = 1
        P%Transfer%PK_redshifts = 0
    end if

    if ((P%NonLinear==NonLinear_lens .or. P%NonLinear==NonLinear_both) .and. P%DoLensing) then
        P%WantTransfer  = .true.
        call Transfer_SetForNonlinearLensing(P%Transfer)
    end if

    call Transfer_SortAndIndexRedshifts(P%Transfer)
    !JD 08/13 end changes

    P%transfer%kmax=P%transfer%kmax*(P%h0/100._dl)

    Ini_fail_on_not_found = .false.

    DebugParam = Ini_Read_Double('DebugParam',DebugParam)
    ALens = Ini_Read_Double('Alens',Alens)

    call Reionization_ReadParams(P%Reion, DefIni)
    call InitialPower_ReadParams(P%InitPower, DefIni, P%WantTensors)
    call Recombination_ReadParams(P%Recomb, DefIni)
    if (Ini_HasKey('recombination')) then
        i = Ini_Read_Int('recombination',1)
        if (i/=1) error stop 'recombination option deprecated'
    end if

    call Bispectrum_ReadParams(BispectrumParams, DefIni, outroot)

    if (P%WantScalars .or. P%WantTransfer) then
        P%Scalar_initial_condition = Ini_Read_Int('initial_condition',initial_adiabatic)
        if (P%Scalar_initial_condition == initial_vector) then
            P%InitialConditionVector=0
            numstr = Ini_Read_String('initial_vector',.true.)
            read (numstr,*) P%InitialConditionVector(1:initial_iso_neutrino_vel)
        end if
        if (P%Scalar_initial_condition/= initial_adiabatic) use_spline_template = .false.
    end if

    if (P%WantScalars) then
        ScalarFileName = trim(outroot)//Ini_Read_String('scalar_output_file')
        LensedFileName =  trim(outroot) //Ini_Read_String('lensed_output_file')
        LensPotentialFileName =  Ini_Read_String('lens_potential_output_file')
        if (LensPotentialFileName/='') LensPotentialFileName = concat(outroot,LensPotentialFileName)
        ScalarCovFileName =  Ini_Read_String_Default('scalar_covariance_output_file','scalCovCls.dat',.false.)
        if (ScalarCovFileName/='') then
            has_cl_2D_array = .true.
            ScalarCovFileName = concat(outroot,ScalarCovFileName)
        end if
    end if

    if (P%WantTensors) then
        TensorFileName =  trim(outroot) //Ini_Read_String('tensor_output_file')
        if (P%WantScalars)  then
            TotalFileName =  trim(outroot) //Ini_Read_String('total_output_file')
            LensedTotFileName = Ini_Read_String('lensed_total_output_file')
            if (LensedTotFileName/='') LensedTotFileName= trim(outroot) //trim(LensedTotFileName)
        end if
    end if

    if (P%WantVectors) then
        VectorFileName =  trim(outroot) //Ini_Read_String('vector_output_file')
    end if

#ifdef WRITE_FITS
    if (P%WantCls) then
        FITSfilename =  trim(outroot) //Ini_Read_String('FITS_filename',.true.)
        if (FITSfilename /='') then
            inquire(file=FITSfilename, exist=bad)
            if (bad) then
                open(unit=18,file=FITSfilename,status='old')
                close(18,status='delete')
            end if
        end if
    end if
#endif


    Ini_fail_on_not_found = .false.

    !optional parameters controlling the computation
    P%AccuratePolarization = Ini_Read_Logical('accurate_polarization',.true.)
    P%AccurateReionization = Ini_Read_Logical('accurate_reionization',.false.)
    P%AccurateBB = Ini_Read_Logical('accurate_BB',.false.)
    P%DerivedParameters = Ini_Read_Logical('derived_parameters',.true.)

    version_check = Ini_Read_String('version_check')
    if (version_check == '') then
        !tag the output used parameters .ini file with the version of CAMB being used now
        call TNameValueList_Add(DefIni%ReadValues, 'version_check', version)
    else if (version_check /= version) then
        write(*,*) 'WARNING: version_check does not match this CAMB version'
    end if
    !Mess here to fix typo with backwards compatibility
    if (Ini_HasKey('do_late_rad_trunction')) then
        DoLateRadTruncation = Ini_Read_Logical('do_late_rad_trunction',.true.)
        if (Ini_HasKey('do_late_rad_truncation')) error stop 'check do_late_rad_xxxx'
    else
        DoLateRadTruncation = Ini_Read_Logical('do_late_rad_truncation',.true.)
    end if

    if (HighAccuracyDefault) then
        DoTensorNeutrinos = .true.
    else
        DoTensorNeutrinos = Ini_Read_Logical('do_tensor_neutrinos',DoTensorNeutrinos )
    end if

    FeedbackLevel = Ini_Read_Int('feedback_level',FeedbackLevel)

    output_file_headers = Ini_Read_Logical('output_file_headers',output_file_headers)

    P%MassiveNuMethod  = Ini_Read_Int('massive_nu_approx',Nu_best)

    ThreadNum      = Ini_Read_Int('number_of_threads',ThreadNum)
    use_spline_template = Ini_Read_Logical('use_spline_template',use_spline_template)

    if (do_bispectrum) then
        lSampleBoost   = 50
    else
        lSampleBoost   = Ini_Read_Double('l_sample_boost',lSampleBoost)
    end if
    if (outroot /= '') then
        if (InputFile /= trim(outroot) //'params.ini') then
            call Ini_SaveReadValues(trim(outroot) //'params.ini',1)
        else
            write(*,*) 'Output _params.ini not created as would overwrite input'
        end if
    end if

    !erick11.04.19. We comment this line to close the inifile after creating the sigma8 & fsigma8 files
    !call Ini_Close

    if (.not. CAMB_ValidateParams(P)) error stop 'Stopped due to parameter error'

#ifdef RUNIDLE
    call SetIdle
#endif

    if (global_error_flag==0) call CAMB_GetResults(P)

    !erick29.01.16. Here we create the file containing the luminosity distance for the BDE model.
    !               We can get dL by two ways:
    !               a) Using directly modules/function LuminosityDistance(z)
    !               b) Using modules/function ComovingRadialDistance(z) and then (1+z)*ComovingRadialDistance(z)
    !               Both methods are consistent -the difference of the values they provide is lower than 1e-7.
    !               The columns are:
    !                         a --- z --- dL (in Mpc) --- dL*H0 (dimensionless) --- mu (distance moduli)
    zSNeIa_max = Ini_Read_Double('zSNeIa_max',1.5_dl)
    n_zSNeIa = 1000   !number of points in the output file

    len_dL_filename = len(trim(outroot)) + 6
    allocate(character(len_dL_filename) :: dL_filename)
    dL_filename = trim(outroot)//'dL.dat'
    call CreateTxtFile(dL_filename,100)

    zSNeIa_min = 1.e-5_dl
    do i_zSNeIa=1, n_zSNeIa
      zSNeIa = zSNeIa_min + (i_zSNeIa-1)*(zSNeIa_max-zSNeIa_min)/(n_zSNeIa-1)
      dL_SNeIa = LuminosityDistance(zSNeIa)
      mu_SNeIa = 5._dl*log10(dL_SNeIa*1.e6_dl/10._dl)
      write(100,'(5E13.5)') 1/(1._dl+zSNeIa), zSNeIa, dL_SNeIa, 1.e5_dl*(P%H0/100)*dL_SNeIa/c, mu_SNeIa
    end do

    !erick16.02.16. Here we create the file containing the r_BAO estimations for the BDE model.
    !               We get r_BAO through the subroutine modules/BAO_D_v and the value of the comoving
    !               sound horizon at the drag epoch computed in modules/inithermo. The columns are:
    !                       a --- z --- rs(z_drag) --- DV(z) --- r_BAO(z)
    !               Since r_BAO is divergent at z=0, we start from z_min=zBAO_max/10.
    zBAO_max = Ini_Read_Double('zBAO_max',1.5_dl)
    n_zBAO = 1000  !number of points in the output file

    len_rBAO_filename = len(trim(outroot)) + 8
    allocate(character(len_rBAO_filename) :: rBAO_filename)
    rBAO_filename = trim(outroot)//'rBAO.dat'
    call CreateTxtFile(rBAO_filename,101)

    zBAO_min = 1.e-5_dl
    do i_zBAO=1, n_zBAO
        zBAO = zBAO_min + (i_zBAO-1)*(zBAO_max-zBAO_min)/(n_zBAO-1)
        DV_BAO = BAO_D_v(zBAO)
        r_BAO = P%rsdrag/DV_BAO
        write(101,'(5E13.5)') 1/(1._dl+zBAO), zBAO, P%rsdrag, DV_BAO, r_BAO
    end do

    if (global_error_flag/=0) then
        write(*,*) 'Error result '//trim(global_error_message)
        error stop
    endif


    if (P%PK_WantTransfer) then
        call Transfer_SaveToFiles(MT,TransferFileNames)
        call Transfer_SaveMatterPower(MT,MatterPowerFileNames)
        call Transfer_output_sig8(MT)
        !erick11.04.19. Once we create the sigma8 & fsigma8 files, we close params.ini
	      call Ini_Close

        !erick16.08.17. Here we create the file containing the information related to LSS according to the
        !               Press-Schechter theory.
        !               The columns are:
        ! log(M_pert [h^-1 M_solar]) --- M_pert [h^-1 M_solar] ---  k_pert [h Mpc^-1] --- R_pert [h^-1 Mpc] --- sigma_R

        !Solar mass in kg
        Msolar = 1.98855E+30_dl
        n_PS = 1000   !number of points in the output file

        !min & max of k (in Mpc**-1)
        kPS_min = 0.001_dl
        kPS_max = 12._dl

        !min & max of N=log(M), where M is in Mpc**-1. In these expressions, rhom0 is in eV**4 so we need to convert it
        !into Mpc**{-4}. The radius scale of the perturbation is given by: R_pert = pi/k_pert. Here
        !                                M_pert = \frac{4}{3}\pi*R_pert^3*\rho_{m0}
        NPS_min = log(4._dl/3._dl*const_pi**4/kPS_max**3*rhom0*eV2toMpcm2**2)
        NPS_max = log(4._dl/3._dl*const_pi**4/kPS_min**3*rhom0*eV2toMpcm2**2)

        len_PressSchechter_filename = len(trim(outroot)) + 18
        allocate(character(len_PressSchechter_filename) :: PressSchechter_filename)
        PressSchechter_filename = trim(outroot)//'PressSchechter.dat'
        call CreateTxtFile(PressSchechter_filename,104)
        do i_PS=1, n_PS
           NPS = NPS_min+(i_PS-1)*(NPS_max-NPS_min)/(n_PS-1)
           MassPS_Mpc = exp(NPS)
           kPS = (4._dl/3._dl*const_pi**4*rhom0/MassPS_Mpc*eV2toMpcm2**2)**(1._dl/3._dl)   !kPS is in Mpc^-1
           RadiusPS = const_pi/kPS            !so far R_pert is in Mpc
           RadiusPS = RadiusPS*(P%H0/100)     !here R_pert is in h^-1 Mpc
           call Transfer_Get_SigmaR(MT,RadiusPS,sigmaPS,8,8,1)
           MassPS_Solar = MassPS_Mpc*(hbar/(c*Mpc*Msolar))  !this is M_pert in M_solar units:
           MassPS_Solar = MassPS_Solar*(P%H0/100)           !this is M_pert in h^-1 M_solar
           write(104,'(5E20.10E2)') log(MassPS_Solar), MassPS_Solar, kPS/(P%H0/100), RadiusPS, sigmaPS(1)
           sigmaPS(1) = 0
        end do
    end if


    if (P%WantCls) then
        call output_cl_files(ScalarFileName, ScalarCovFileName, TensorFileName, TotalFileName, &
            LensedFileName, LensedTotFilename, output_factor)

        call output_lens_pot_files(LensPotentialFileName, output_factor)

        if (P%WantVectors) then
            call output_veccl_files(VectorFileName, output_factor)
        end if

#ifdef WRITE_FITS
        if (FITSfilename /= '') call WriteFitsCls(FITSfilename, CP%Max_l)
#endif
    end if

    call CAMB_cleanup
    stop

100 stop 'Must give num_massive number of integer physical neutrinos for each eigenstate'
    end program driver


#ifdef RUNIDLE
    !If in Windows and want to run with low priorty so can multitask
    subroutine SetIdle
    USE DFWIN
    Integer dwPriority
    Integer CheckPriority

    dwPriority = 64 ! idle priority
    CheckPriority = SetPriorityClass(GetCurrentProcess(), dwPriority)

    end subroutine SetIdle
#endif
