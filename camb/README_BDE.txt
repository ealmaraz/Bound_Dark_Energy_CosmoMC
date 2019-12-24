====================
CAMB implementation of the Bound Dark Energy model
	V(phi) = cpot * Lambdac^{4+alpha}/phi^alpha
Original CAMB distribution   : november 2016 + equations_quint.f90 template
Modifications : indicated by 'erickdd.mm.yy'
Comments      : indicated by 'erick - comment'
====================


Change log 
* If use_physical = True, the basic parameters of the model are: 
	omegabh2, omegach2, Lambdac, ns, As & tau
  and H0 is derived. 
  If use_physical = False, the basic parameters of the model are: 
	omegabh2, omegach2, H0, ns, As & tau
  and Lambdac is derived.
* The input parameters of the Dark Group are cpot (the proportionality constant of the potential),
  wDGc (the EoS at ac in units of 1/3), and the colors and flavors of the Dark Group
* Format for the name of the output files from output_root in params.ini 
* Vector and tensor perturbations are not supported



Modified files (main modifications)
* equations_bde_v3.f90: We introduce/modify the subroutines/functions:
			    bde_h0, bde_omega_today, bde_dlsoda_equations, jac, bde_search_lambdac, 
			    bde_dverk_equations, bde_v, bde_background, dtauda, 
			    GaugeInterface_EvolveScal, output, initial, and derivs
			We use the fortran ODEPACK library (specifically, the DLSODA solver) to deal 
			with the initial stiffness of the system of ODEs describing the background
			evolution of BDE after condensation. Then we use dverk implemented in 
                        $CAMB/subroutines.f90 to solve the system up to present time. On the effect
			of different numerical integrators in cosmological quantities, see
			analysis/integrators_test.ipynb

*inidriver.F90: We introduce code to read BDE parameters and support for use_physical = F. We introduce
		code to create output files for SNeIa, BAO and Press-Schechter information.

* cmbmain: We introduce code to compute the growth factors and compute linear perturbations in
	    the Newtonian gauge (in subroutine CalcScalarSources)

* modules.f90: Some slight modifications

* params.ini: We specify the input parameters of the model

* Makefile: Here we link ODEPACK to CAMB (see the instructions to accomplish this)



Outputs
Besides the usual output files containing the info of the CMB, P(k), and transfers, we also create the 
following files:

1) Background information is saved in output_root_background.dat containing the data in the following format:
     a --- N --- x --- y --- l --- phi_camb --- phidot_camb --- V_camb --- Omegar --- Omegam --- OmegaDE --- wDE 
If we set use_physical = F, we append a suffix 'nophysical' after 'background'

2) SNeIa information is saved in output_root_dL.dat containing the data in the following format:
     a --- z --- dL (in Mpc) --- dL*H0 (dimensionless) --- mu (distance moduli)

3) BAO information is saved in output_root_rBAO.dat containing the data in the following format:
     a --- z --- rs(z_drag) (in Mpc) --- DV(z) (in Mpc) ---  r_BAO(z) (dimensionless)

4) sigma8 & fsigma8 information is saved in output_root_sigma8.dat and output_root_fsigma8.dat containing the 
   data in the following format:
     z --- sigma8 (fsigma8)

5) Press-Schechter information is saved in output_root_PressSchechter.dat containing the data in the 
following format:
     log(M_pert [h^-1 M_solar]) --- M_pert [h^-1 M_solar] ---  k_pert [h Mpc^-1] --- R_pert [h^-1 Mpc] --- sigma_R 
     
6) The information of the perturbations is saved in output_root_k_<mode>.dat containing the data in the following format: 
        tau    ---  a   --- keta  ---      dc     ---    db       ---  vb   --- dg       ---    qg    
        dn    ---   qn    --- y(EV%r_ix+2) ---  y(EV%r_ix+3) --- dphi --- dphi' --- deltam(syn) 
        growth(syn) --- Hconf --- dc'(Syn) --- db'(Syn) --- dc(new) --- db(new) --- delgtam(new) 
        dg(new) --- dn(new)

7) The information of the growth factors is saved in output_root_growth_factors_k_<mode>.dat containing the data 
   in the following format: 
     tau --- a --- z --- Hconf --- D1 --- a*dD1/dtau --- D2 --- a*dD2/dtau


LAST MODIFICATION: 22nd September 2019

