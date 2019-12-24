====================
CosmoMC implementation of the Bound Dark Energy model
        V(phi) = cpot * Lambdac^{4+alpha}/phi^alpha
Original CosmoMC distribution   : November 2016
CAMB distribution: v3.1.4 version (September 2019)
Modifications : indicated by 'erickdd.mm.yy'
Comments      : indicated by 'erick - comment'
====================


Change log 
* The hard parameters of BDE are: 
        omegabh2, omegach2, Lambdac, cpot, ns, As & tau
  H0 and 100*theta_{MC} are derived. 
* The input parameters of the Dark Group are wDGc (the EoS at ac in units of 1/3), and 
  the colors and flavors of the Dark Group
* use_physical = F is not supported yet


Modified files (main modifications)
*source: Minor modifications in driver.F90, settings.f90 & Makefile.
         In CosmologyTypes.f90 we declare BDE and other extra parameters as well output
         z's for LSS and BAO information.
         In Calculator_CAMB.f90 we connect CAMB to CosmoMC variables.
         In CosmologyParameterizations.f90 we solve the background evolution and compute
         H0, wBDE0 & 100*theta_{MC}. The hard work is done here!

* camb: Minor modifications in cmbmain.f90, inidriver.F90 & modules.f90.
        Slight modifications in equations_bde_v3.f90: these modifications are aimed to plug
        our CAMB distribution into CosmoMC without introducing too much extra code.
        We note that the BDE subroutines to compute the background evolution now run in 
        source/CosmologyParameterizations.f90
        
* odepack: The odepack folder with the DLSODA subroutines is placed in the main directory
           of CosmoMC. See source/Makefile

* paramnames: We introduce/set in order basic and derived parameters in 
              params_CMB_BDE.paramnames, derived_theory.paramnames & derived_LSS.paramnames 
              

LAST MODIFICATION: 30th November 2019


