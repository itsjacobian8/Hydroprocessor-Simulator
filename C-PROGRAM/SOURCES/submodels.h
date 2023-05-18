#ifndef submodels_h_
#define submodels_h_

#include<stddef.h>

// Gas Phase Properties
typedef struct gas_phase
{
    double  qgtMin;
    double  qgtMax;
    double  qgtIncrement;
    double  rhog;
} gasPhase;

typedef struct liquid_phase
{
    double  qlvr;
    double  rhol;
    double  mul;
} liquidPhase;

typedef struct solid_phase
{
    double  mcat;
    double  Vp;
    double  Ap;
    double  rhop;  
} solidPhase;

typedef struct Interfacial
{
    double  sigma;
} interfacial;

typedef struct Column
{
    double  P;
    double  Dc;
    double  drec;
    double  vsep;
    double  hset;
} column;

typedef struct Grid
{
    double  nr;
    double  nor;
    double  dor;
} grid;

typedef struct bubble_size_distribution
{
    double  classWidth;
    double  minCutOff;
    double  maxCutOff;
} bubbleSizeDistribution;

typedef struct Separator
{
    double  beta1;
    double  beta2;
    double  Rstep;
} separator;

typedef struct Convergence
{
    double  absTol;
    double  relTol;
    size_t  maxIter;
} convergence;


/**
 * @brief Computes volume equivalent particle diameter given particle volume
 * 
 * @param Vp particle volume [m^3]
 * @return double 
 */
double volumeEquivalentDiameter
(
    double Vp
);

/**
 * @brief Computes particle sphericity given particle volume and surface area
 * 
 * @param Vp particle volume [m^3]
 * @param Ap particle surface area [m^2]
 * @return double 
 */
double sphericity
(
    double Vp, 
    double Ap
);

/**
 * @brief Computes particle Archimedes number
 * 
 * @param rhop particle density [kg/m^3]
 * @param rhol fluid (i.e. liquid) density [kg/m^3]
 * @param mul fluid (i.e. liquid) viscosity [Pa.s]
 * @param dv volume equivalent particle diameter [m]
 * @return double 
 */
double ArchimedesNumber
(
    double rhop, 
    double rhol, 
    double mul, 
    double dv
);

/**
 * @brief Computes void fraction at minimum fluidization
 * 
 * @param phi particle sphericity [-]
 * @return double 
 */
double minimumVoidFraction
(
    double phi
);

/**
 * @brief Computes two phase (i.e. liquid-solid) minimum fluidization velocity 
 * 
 * @param rhol fluid (i.e. liquid) density [kg/m^3]
 * @param mul fluid (i.e. liquid) viscosity [Pa.s]
 * @param dv volume equivalent particle diameter [m]
 * @param phi particle sphericity [-]
 * @param Arp particle Archimedes number [-]
 * @param emfo void fraction at minimum fluidization [-]
 * @return double 
 */
double LSMinimumFluidizationVelocity
(
    double rhol, 
    double mul, 
    double dv, 
    double phi, 
    double Arp, 
    double emfo
);

/**
 * @brief Computes three phase (gas-liquid-solid) minimum fluidization velocity
 * 
 * @param rhop particle density [kg/m^3]
 * @param rhol liquid density [kg/m^3]
 * @param rhog gas density [kg/m^3]
 * @param mul liquid viscosity [Pa.s]
 * @param Ugbed bed superficial gas velocity [m/s]
 * @param dv volume equivalent particle diameter [m]
 * @param phi particle sphericity [-]
 * @param emfo void fraction at minimum fluidization (two phase system) [-]
 * @param Umfo two phase minimum fluidization velocity [m/s]
 * @param relTol relative tolerance [-]
 * @return double 
 */
double GLSMinimumFluidizationVelocity
(
    double rhop, 
    double rhol, 
    double rhog, 
    double mul, 
    double Ugbed, 
    double dv, 
    double phi, 
    double emfo, 
    double Umfo, 
    double relTol
);

/**
 * @brief Computes wall effect parameter
 * 
 * @param dv volume equivalent particle diameter [m]
 * @param Dc column's diameter [m]
 * @return double 
 */
double wallEffectParameter
(
    double dv, 
    double Dc
);

/**
 * @brief Computes Richardson-Zaki coefficient
 * 
 * @param Arp particle Archimedes number [-]
 * @param dv volume equivalent particle diameter [m]
 * @param Dc column's diameter [m]
 * @return double 
 */
double RichardsonZakiCoefficient
(
    double Arp, 
    double dv, 
    double Dc
);

/**
 * @brief Computes particle terminal velocity
 * 
 * @param rhol liquid density [kg/m^3]
 * @param mul liquid viscosity [Pa.s]
 * @param Arp particle Archimedes number [-]
 * @param dv volume equivalent particle diameter [m]
 * @param phi particle sphericity [-]
 * @return double 
 */
double terminalVelocity
(
    double rhol, 
    double mul, 
    double Arp, 
    double dv, 
    double phi
);

/**
 * @brief Computes particles holdup
 * 
 * @param Ugbed bed superficial gas velocity [m/s]
 * @param Ulbed bed superficial liquid velocity [m/s]
 * @param k wall effect parameter [-]
 * @param n Richardson-Zaki Coefficient [-]
 * @param upinf particle terminal velocity [m/s]
 * @return double 
 */
double particlesHoldup
(
    double Ugbed, 
    double Ulbed, 
    double k, 
    double n, 
    double upinf
);

/**
 * @brief Computes bed height
 * 
 * @param drec recycle line diameter [m]
 * @param Dc column diameter [m]
 * @param epsilonp particles holdup [-]
 * @param mcat mass of catalyst inventory [kg]
 * @param rhop particle density [kg/m^3]
 * @return double 
 */
double bedHeight
(
    double drec, 
    double Dc, 
    double epsilonp, 
    double mcat, 
    double rhop
);

/**
 * @brief Chebyshev Polynomial Relation 3 - degree 7 for numerator/denominator
 * 
 * @param z z-score
 * @return double 
 */
double R377(double z);

/**
 * @brief Chebyshev Polynomial Relation 4 - degree 7 for numerator/denominator
 * 
 * @param z z-score
 * @return double 
 */
double R477(double z);

/**
 * @brief Chebyshev Polynomial Relation 5 - degree 7 for numerator/denominator
 * 
 * @param z z-score
 * @return double 
 */
double R577(double z);

/**
 * @brief Computes standard normal quantile
 * 
 * @param cProb cumulative probability
 * @return double 
 */
double standardNormalQuantile(double cProb);

/**
 * @brief Computes lognormal quantile
 * 
 * @param cProb cumulative probability
 * @param mu lognormal mean
 * @param s lognormal standard deviation
 * @return double 
 */
double lognormalQuantile(double cProb, double mu, double s);

/**
 * @brief Computes lognormal probability density
 * 
 * @param cbi bubble chord length
 * @param mu lognormal mean
 * @param s lognormal standard deviation
 * @return double 
 */
double lognormalPDF(double cbi, double mu, double s);

/**
 * @brief Computes lognormal distribution parameters
 * 
 * @param rhol liquid density
 * @param rhog gas density
 * @param mul liquid viscosity
 * @param Ugbed bed superficial gas velocity
 * @param Ulbed bed superficial liquid velocity
 * @param dor distributor outlet orifice diameter
 * @param nor number of outlet orifices per riser
 * @param nr number of risers
 * @param Dc column diameter
 * @param drec recycle line diameter
 * @param mu lognormal mean(to be computed)
 * @param s lognormal standard deviation(to be computed)
 */
void lognormalDistributionParameters
(
    double rhol, 
    double rhog, 
    double mul, 
    double Ugbed, 
    double Ulbed, 
    double dor, 
    double nor, 
    double nr, 
    double Dc, 
    double drec, 
    double *mu, 
    double *s
);

/**
 * @brief Converts bubble chord length to volume equivalent diameter
 * 
 * @param cbi bubble chord length
 * @param rhol liquid density
 * @param rhog gas density
 * @param sigma gas-liquid interfacial tension
 * @param relTol relative tolerance
 * @return double 
 */
double chordToDiameter
(
    double cbi, 
    double rhol, 
    double rhog, 
    double sigma, 
    double relTol
);

/**
 * @brief Computes bubble drag coefficient
 * 
 * @param dbi bubble diameter
 * @param egi individual gas holdup
 * @param usi bubble slip velocity
 * @param rhol liquid density
 * @param rhog gas density
 * @param sigma gas-liquid interfacial tension
 * @param mul liquid viscosity
 * @param p dimensionless operating pressure
 * @return double 
 */
double MachDrag
(
    double dbi, 
    double egi, 
    double usi, 
    double rhol, 
    double rhog, 
    double sigma, 
    double mul, 
    double p
);

/**
 * @brief Computes bubble slip velocity
 * 
 * @param dbi bubble diameter
 * @param egi individual gas holdup
 * @param rhol liquid density
 * @param rhog gas density
 * @param sigma gas-liquid interfacial tension
 * @param mul liquid viscosity
 * @param p dimensionless operating pressure
 * @param relTol relative tolerance
 * @return double 
 */
double bubbleSlipVelocity(double dbi, double egi, double rhol, double rhog, double sigma, double mul, double p, double relTol);

/**
 * @brief Objective function for iteratively computing individual gas holdups
 * 
 * @param egi individual gas holup (initial guess)
 * @param dbi bubble diameter
 * @param rhol liquid density
 * @param rhog gas density
 * @param sigma gas-liquid interfacial tension
 * @param mul liquid viscosity
 * @param p dimensionless operating pressure
 * @param Ugbed bed superficial gas velocity
 * @param Ulbed bed superficial liquid velocity
 * @param relTol relative tolerance
 * @return double 
 */
double objectiveFunction
(
    double egi, 
    double dbi, 
    double rhol, 
    double rhog, 
    double sigma, 
    double mul, 
    double p, 
    double Ugbed, 
    double Ulbed, 
    double relTol
);

/**
 * @brief Custom bisection method for computing individual gas holdups
 * 
 * @param dbi bubble diameter
 * @param rhol liquid density
 * @param rhog gas density
 * @param sigma gas-liquid interfacial tension
 * @param mul liquid viscosity
 * @param p dimensionless operating pressure
 * @param Ugbed bed superficial gas velocity
 * @param Ulbed bed superficial liquid velocity
 * @param absTol absolute tolerance
 * @param relTol relative tolerance
 * @return double 
 */
double bisect
(
    double dbi, 
    double rhol, 
    double rhog, 
    double sigma, 
    double mul, 
    double p, 
    double Ugbed, 
    double Ulbed, 
    double absTol, 
    double relTol
);

/**
 * @brief Generates bubble size and slip velocity distribution
 * 
 * @param rhol liquid density
 * @param rhog gas density
 * @param sigma gas-liquid interfacial tension
 * @param mul liquid viscosity
 * @param p dimensionless operating pressure
 * @param Ugbed bed superficial gas velocity
 * @param Ulbed bed superficial liquid velocity
 * @param cbiMin minimum bubble chord length
 * @param classWidth bubble class width
 * @param nClasses number of bubble classes
 * @param absTol absolute tolerance 
 * @param relTol relative tolerance
 * @param mu lognormal mean
 * @param s lognormal standard deviation
 * @param bsvd bubble size and slip velocity distribution (to be filled)
 */
void bubbleSizeAndVelocityDistribution
(
    double rhol, 
    double rhog, 
    double sigma, 
    double mul, 
    double p, 
    double Ugbed, 
    double Ulbed, 
    double cbiMin, 
    double classWidth, 
    size_t nClasses, 
    double absTol, 
    double relTol, 
    double mu, 
    double s, 
    double **bsvd
);

/**
 * @brief Computes gas-liquid separation efficiency of the internal liquid recycle pan
 * 
 * @param rhol liquid density
 * @param rhog gas density
 * @param mul liquid viscosity
 * @param sigma gas-liquid interfacial tension
 * @param Ugbed bed superficial gas velocity
 * @param Ulbed bed superficial liquid velocity
 * @param dor distributor outlet orifice diameter
 * @param nor number of outlet orifices per riser
 * @param nr number of risers
 * @param Dc column diameter
 * @param drec recycle line diameter
 * @param R liquid recycle ratio
 * @param vsep gas-liquid separator volume
 * @param p dimensionless operating volume
 * @param classWidth bubble class width
 * @param beta1 geometry dependent parameter
 * @param beta2 bubble and geometry dependent parameter
 * @param minCutOff minimum cutoff cumulative probability
 * @param maxCutOff maximum cutoff cumulative probability
 * @param absTol absolute tolerance
 * @param relTol relative tolerance
 * @return double 
 */
double GLSeparator
(
    double rhol, 
    double rhog, 
    double mul, 
    double sigma, 
    double Ugbed, 
    double Ulbed, 
    double dor, 
    double nor, 
    double nr, 
    double Dc, 
    double drec, 
    double R, 
    double vsep, 
    double p, 
    double classWidth, 
    double beta1, 
    double beta2, 
    double minCutOff, 
    double maxCutOff,  
    double absTol, 
    double relTol
);

/**
 * @brief Updates recycle gas and liquid flow rates
 * 
 * @param Ugbed bed superficial gas velocity
 * @param Ulbed bed superficial liquid velocity
 * @param Dc column diameter
 * @param drec recycle line diameter
 * @param R liquid recycle ratio
 * @param eta gas-liquid separation efficiency
 * @param qgr gas recycle flow rate
 * @param qlr liquid recycle flow rate
 */
void recycleFlowRates
(
    double Ugbed, 
    double Ulbed, 
    double Dc, 
    double drec, 
    double R, 
    double eta, 
    double *qgr, 
    double *qlr
);

/**
 * @brief Updates bed gas and liquid flow rate
 * 
 * @param qgt fresh treat gas flow rate
 * @param qlvr fresh bitumen flow rate
 * @param Dc column diameter
 * @param drec recycle line diameter
 * @param qgr recycle gas flow rate
 * @param qlr recycle liquid flow rate
 * @param Ugbed bed superficial gas velocity
 * @param Ulbed bed superficial liquid velocity
 */
void bedFlowRates
(
    double qgt, 
    double qlvr, 
    double Dc, 
    double drec, 
    double qgr, 
    double qlr, 
    double *Ugbed, 
    double *Ulbed
);

/**
 * @brief Computes freeboard gas holdup
 * 
 * @param rhol liquid density
 * @param rhog gas density
 * @param sigma gas-liquid interfacial tension
 * @param mul liquid viscosity
 * @param p dimensionless operating pressure
 * @param Ugbed bed superficial gas velocity
 * @param Ulbed bed superficial liquid velocity
 * @param dor distributor outlet orifice diameter
 * @param nor number of outlet orifices per riser
 * @param nr number of risers
 * @param Dc column diameter
 * @param drec recycle line diameter
 * @param absTol absolute tolerance
 * @param relTol relative tolerance
 * @param minCutOff minimum cutoff cumulative probability
 * @param maxCutOff maximum cutoff cumulative probability
 * @param classWidth bubble class width
 * @return double 
 */
double freeboardGasHoldup
(
    double rhol, 
    double rhog, 
    double sigma, 
    double mul, 
    double p, 
    double Ugbed, 
    double Ulbed, 
    double dor, 
    double nor, 
    double nr, 
    double Dc, 
    double drec, 
    double absTol,
    double relTol,
    double minCutOff,
    double maxCutOff, 
    double classWidth
);

void readGasProperties(gasPhase *properties, char *input_dir);

void readLiquidProperties(liquidPhase *properties, char *input_dir);

void readSolidProperties(solidPhase *properties, char *input_dir);

void readInterfacialProperties(interfacial *properties, char *input_dir);

void readColumnParameters(column *parameters, char *input_dir);

void readGridParameters(grid *parameters, char *input_dir);

void readBubbleSizeDistributionParameters(bubbleSizeDistribution *parameters, char *input_dir);

void readSeparatorParameters(separator *parameters, char *input_dir);

void readConvergenceCriteria(convergence *criteria, char *input_dir);

/**
 * @brief Iterative solver using known mass of catalyst inventory for convergence
 * 
 * @param params struct of input parameters
 * @param output_dir output directory where bubble size and velocity distribution data will be stored
 */
void solve
(
    gasPhase gasProperties,
    liquidPhase liquidProperties,
    solidPhase  solidProperties,
    interfacial interfcialProperties,
    column  columnParameters,
    grid    gridParameters,
    bubbleSizeDistribution  bubbleSizeDistributionParameters,
    separator   separatorParameters,
    convergence convergenceCriteria,
    char *output_dir
);

#endif