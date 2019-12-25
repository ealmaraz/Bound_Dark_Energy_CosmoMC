# Bound Dark Energy #

This code presents the implementation in CosmoMC of the Bound Dark Energy model studied in:

1. A. de la Macorra and E. Almaraz, <em> Phys. Rev. Lett. </em> **121**, 161303 (2018); https://arxiv.org/abs/1805.01510
2. E. Almaraz and A. de la Macorra, <em> Phys. Rev. D </em> **99**, 103504 (2019); https://arxiv.org/abs/1812.01133,

where Dark Energy is described by a dynamically-generated scalar field with an inverse power-law potential:
<img src="https://latex.codecogs.com/svg.latex?\Large&space;V(\phi)=c\frac{\Lambda_c^{4+\alpha}}{\phi^{\alpha}}" title="\Large V(\phi)=c\frac{\Lambda_c^{4+\alpha}}{\phi^{\alpha}}" />

Here <img src="https://latex.codecogs.com/svg.latex?\Large&space;c" title="\Large c" /> is a constant, <img src="https://latex.codecogs.com/svg.latex?\Large&space;\Lambda_c" title="\Large \Lambda_c" /> is the condensation scale of the scalar field and <img src="https://latex.codecogs.com/svg.latex?\Large&space;\alpha=2/3" title="\Large \alpha=2/3" />.

If you use this code, I would be very grateful if you cite these two papers. 

### Requisites ###
1. ifort (I use version 15.0.7 20160518)
2. odepack
3. Open MPI

### Installation ###
1. Go to `$CosmoMC`
2. Download `odepack.f`, `odepack_sub1.f` & `odepack_sub2.f` from:<br>
  https://people.sc.fsu.edu/~jburkardt/f77_src/odepack/odepack.html <br>
  and place these files in `$CosmoMC/odepack`
3. Compile ODEPACK:
```
   cd $CosmoMC/odepack
   ifort -c odepack.f odepack_sub1.f odepack_sub2.f
   ar qc libodepack.a *.o
```                                
4. Link these libraries to CosmoMC (see Makefile)
5. Compile CosmoMC. Open `$CosmoMC/source` and type 
```
make
```
6. Run the code:
```
mpirun -np <number_of_chains> ./cosmomc test.ini
```

Note that in order to use Planck data you have to install the likelihood code first. See CosmoMC website for more information: https://cosmologist.info/cosmomc/

### Change log ###
Modifications and comments are indicated by `erickdd.mm.yy`. Main modifications done in:
1. `source/CosmologyTypes.f90`
2. `source/Calculator_CAMB.f90`
3. `source/CosmologyParameterizations.f90`
4. `paramnames/params_CMB_BDE.paramnames`

See README_BDE.txt for more details.<br> 
`$CosmoMC/chains` contains the chains leading to the results reported in [1] & [2] <br>
`$CosmoMC/chains/bestfits` contains two best fits: `v3` corresponding to an older version of the code (also reported in [1] & [2]) and `v3.1.4`, which was computed using this version. We have set <img src="https://latex.codecogs.com/svg.latex?\Large&space;c" title="\Large c" /> = 1 in our analysis. <br>

The data used to find the constraints are: <br>

3. Ade, P. A. R. et al, <em> Astron. Astrophys. </em> **594**, A13 (2016)
4. Betoule, M. et al, <em> Astron. Astrophys. </em> **568**, A22 (2014)
5. Ross, A. J. et al, <em> Mon. Not. R. Astron. Soc. </em> **449**, 835-847 (2015) 
6. Beutler, F. et al, <em> Mon. Not. R. Astron. Soc.</em> **416**, 3017-3032 (2011)
7. Gil-Marín, H. et al, <em> Mon. Not. R. Astron. Soc.</em> **460**, 4210-4219 (2016)

### Acknowledgments ###
This code is built upon the CosmoMC public code (November 2016 release; https://cosmologist.info/cosmomc/) and our BDE implementation in CAMB: https://github.com/ealmaraz/Bound_Dark_Energy_CAMB. The ODEPACK Fortran library is also used: https://people.sc.fsu.edu/~jburkardt/f77_src/odepack/odepack.html

### Contact ###
For comments, suggestions, etc, feel free to contact me: erickalmaraz@gmail.com
