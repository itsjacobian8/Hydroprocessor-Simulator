#include<math.h>
#include<assert.h>
#include<stdlib.h>
#include<string.h>
#include<stdio.h>
#include<locale.h>
#include<yaml.h>
#include<time.h>

#include "submodels.h"

#define pi  M_PI
#define g   9.81

double volumeEquivalentDiameter
(
    double Vp
)
{
    return cbrt((6.0*Vp)/pi);
}

double sphericity
(
    double Vp, 
    double Ap
)
{
    return cbrt(pi*pow(6.0*Vp, 2.0))/Ap;
}

double ArchimedesNumber
(
    double rhop, 
    double rhol, 
    double mul, 
    double dv
)
{
    return (rhol*(rhop - rhol)*g*pow(dv,3.0))/pow(mul, 2.0);
}

double minimumVoidFraction
(
    double phi
)
{
    return 0.415/cbrt(phi);
}

double LSMinimumFluidizationVelocity
(
    double rhol, 
    double mul, 
    double dv, 
    double phi, 
    double Arp, 
    double emfo
)
{
    // constants
    double C1 = 150.0*(1.0 - emfo)/(3.5*phi);
    double C2 = pow(emfo, 3.0)/1.75;

    // particle Reynolds number at minimum fluidization
    double Remfo = sqrt(pow(C1, 2.0) + C2*Arp) -  C1;

    return (mul*Remfo)/(rhol*dv);
}

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
)
{
    double relError, Umf_prev, emf, alphamf, Argls, C1, C2, Remf;

    // initialize three-phase minimum fluidization velocity using two-phase value
    double Umf = Umfo;

    do
    {
        // assign current minimum fluidization velocity to previous value
        Umf_prev = Umf;

        // void fraction at minimum fluidization
        emf = emfo*(1.0 - 0.34*(1.0 - Umf_prev/Umfo) + 0.22*pow(1.0-Umf_prev/Umfo,2.0));

        // gas holdup at minimum fluidization
        alphamf = (0.16*Ugbed)/(emf*(Ugbed + Umf_prev));

        // Archimedes number
        Argls = (rhol*(rhop - (rhog*alphamf + rhol*(1.0-alphamf)))*g*pow(dv,3.0))/pow(mul,2.0);

        // constants
        C1 = (150.0*(1.0 - emf))/(3.5*phi);
        C2 = (pow(emf,3.0)*pow(1.0 - alphamf,3.0))/1.75;

        // minimum Reynolds number
        Remf = sqrt(pow(C1,2.0) + C2*Argls) - C1;

        // update minimum fluidization velocity
        Umf = (mul*Remf)/(rhol*dv);

        // relative error on minimum fluidization velocity
        relError = fabs((Umf - Umf_prev)/Umf_prev);
    } while (relError > relTol);
    
    return Umf;
}

double wallEffectParameter
(
    double dv, 
    double Dc
)
{
    return 1.0 - 1.33*(dv/Dc);
}

double RichardsonZakiCoefficient
(
    double Arp, 
    double dv, 
    double Dc
)
{
    double A = 0.043*pow(Arp,0.57)*(1.0 - 1.24*pow(dv/Dc,0.27));
    return (4.8 + 2.4*A)/(1.0 + A);
}

double terminalVelocity
(
    double rhol, 
    double mul, 
    double Arp, 
    double dv, 
    double phi
)
{
    // particle Reynolds number at terminal conditions
    double Repinf = cbrt(Arp)/(18.0/pow(Arp, 2.0/3.0) + (2.335 - 1.744*phi)/pow(Arp, 1.0/6.0));
    
    return (mul*Repinf)/(rhol*dv);
}

double particlesHoldup
(
    double Ugbed, 
    double Ulbed, 
    double k, 
    double n, 
    double upinf
)
{
    return 1.0 - (pow(Ulbed/(k*upinf),1.0/n) * (1.0 + 0.22*pow(Ugbed/Ulbed,0.92)));
}

double bedHeight
(
    double drec, 
    double Dc, 
    double epsilonp, 
    double mcat, 
    double rhop
)
{
    // column's cross-sectional area
    double Ac = 0.25*pi*(pow(Dc,2.0) - pow(drec,2.0));

    return mcat/(rhop*Ac*epsilonp);
}

double R377(double z)
{
    double  numerator   =   2.5090809287301226727e3*pow(z,7) + 
                            3.3430575583588128105e4*pow(z,6) + 
                            6.7265770927008700853e4*pow(z,5) + 
                            4.5921953931549871457e4*pow(z,4) +
                            1.3731693765509461125e4*pow(z,3) + 
                            1.9715909503065514427e3*pow(z,2) + 
                            1.3314166789178437745e2*z +
                            3.3871328727963666080;
    
    double  denominator =   5.2264952788528545610e3*pow(z,7) + 
                            2.8729085735721942674e4*pow(z,6) + 
                            3.9307895800092710610e4*pow(z,5) + 
                            2.1213794301586595867e4*pow(z,4) + 
                            5.3941960214247511077e3*pow(z,3) + 
                            6.8718700749205790830e2*pow(z,2) + 
                            4.2313330701600911252*10*z + 
                            1.0000000000000000000;

    return numerator/denominator;
}

double R477(double z)
{
    double  numerator   =   2.01033439929228813265e-7*pow(z,7) + 
                            2.71155556874348757815e-5*pow(z,6) + 
                            1.24266094738807843860e-3*pow(z,5) + 
                            2.65321895265761230930e-2*pow(z,4) +
                            2.96560571828504891230*0.1*pow(z,3) + 
                            1.78482653991729133580*pow(z,2) + 
                            5.46378491116411436990*z +
                            6.65790464350110377720;
    
    double  denominator =   2.04426310338993978564e-15*pow(z,7) + 
                            1.42151175831644588870e-7*pow(z,6) + 
                            1.84631831751005468180e-5*pow(z,5) + 
                            7.86869131145613259100e-4*pow(z,4) + 
                            1.48753612908506148525e-2*pow(z,3) + 
                            1.36929880922735805310*0.1*pow(z,2) + 
                            5.99832206555887937690*0.1*z + 
                            1.00000000000000000000;

    return numerator/denominator;
}

double R577(double z)
{
    double  numerator   =   7.74545014278341407640e-4*pow(z,7) + 
                            2.27238449892691845833e-2*pow(z,6) + 
                            2.41780725177450611770*0.1*pow(z,5) + 
                            1.27045825245236838258*pow(z,4) +
                            3.64784832476320460504*pow(z,3) + 
                            5.76949722146069140550*pow(z,2) + 
                            4.63033784615654529590*z +
                            1.42343711074968357734;
    
    double  denominator =   1.05075007164441684324e-9*pow(z,7) + 
                            5.47593808499534494600e-4*pow(z,6) + 
                            1.51986665636164571966e-2*pow(z,5) + 
                            1.48103976427480074590*0.1*pow(z,4) + 
                            6.89767334985100004550*0.1*pow(z,3) + 
                            1.67638483018380384940*pow(z,2) + 
                            2.05319162663775882187*z + 
                            1.00000000000000000000;

    return numerator/denominator;
}

double standardNormalQuantile(double cProb)
{
    double stdnormQuantile, z;

    double limit1 = 0.425;
    double limit2 = 1e-25;
    double limit3 = 0.075;
    double limit4 = 0.925;
    double limit5 = (1.0 - 1e-25);

    if(fabs(cProb - 0.5) < limit1)
    {
        z = pow(0.425, 2.0) - pow(cProb - 0.5, 2.0);
        stdnormQuantile = (cProb - 0.5)*R377(z);
    }
    else if(cProb <= limit2)
    {
        z = sqrt(-log(cProb)) - 5.0;
        stdnormQuantile = -R477(z);
    }
    else if(cProb > limit2 && cProb < limit3)
    {
        z = sqrt(-log(cProb)) - 1.6;
        stdnormQuantile = -R577(z);
    }
    else if(cProb > limit4 && cProb < limit5)
    {
        z = sqrt(-log(1.0 - cProb)) - 1.6;
        stdnormQuantile = R577(z);
    }
    else if (cProb >= limit5)
    {
        z = sqrt(-log(1.0 - cProb)) - 5.0;
        stdnormQuantile = R477(z);
    }

    return stdnormQuantile;
}

double lognormalQuantile
(
    double cProb, 
    double mu, 
    double s
)
{
    // compute standard normal quantile
    double stdnormQuantile = standardNormalQuantile(cProb);
    
    return exp(mu + s*stdnormQuantile);
}

double lognormalPDF
(
    double cbi, 
    double mu, 
    double s
)
{
    return 1.0/(cbi*s*sqrt(2.0*pi)) * exp(-pow(log(cbi) - mu, 2.0)/(2.0*pow(s, 2.0)));
}

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
)
{
    // cross-sectional area of outlet orifice
    double Aor = 0.25*pi*pow(dor,2.0);

    // cross-sectional area of the column
    double Ac = 0.25*pi*(pow(Dc,2.0) - pow(drec,2.0));

    // orifice-based superficial gas velocity
    double Ugor = Ugbed/(nr*nor)*(Ac/Aor);

    // orifice-based superficial liquid velocity
    double Ulor = Ulbed/(nr*nor)*(Ac/Aor);

    // orifice-based dimensionless groups
    double Re = (rhol*Ulor*dor)/mul;
    double Fr = Ugor/sqrt(g*dor);
    double rhoRatio = rhog/rhol;

    // dimensionless lognormal mean
    double muStar = -0.5305 - 0.1469*log(rhoRatio) - 0.4229*log(Re) + 0.345*log(Fr);

    // lognormal mean
    *mu =  muStar + log(1000.0*dor);

    // lognormal normal standard deviation
    *s = 0.7594 - 0.0426*log(rhoRatio) - 0.0201*log(Re) + 0.0596*log(Fr);
}

double chordToDiameter
(
    double cbi, 
    double rhol, 
    double rhog, 
    double sigma,
    double relTol
)
{
    // Eotvos number
    double Eo = ((rhol - rhog)*g*pow(0.0015*cbi, 2.0))/sigma;

    // bubble aspect ratio
    double E = 1.0/(1.0 + 0.163*pow(Eo, 0.757));

    // initial estimate of bubble diameter
    double dbi = 1.5*cbi*pow(E, -2.0/3.0);

    double dbi_prev, relError;
    do
    {
        // assign current bubble diameter estimate to previous value
        dbi_prev = dbi;

        // update Eotvos number
        Eo = ((rhol - rhog)*g*pow(0.001*dbi_prev, 2.0))/sigma;

        // update bubble aspect ratio
        E = 1.0/(1.0 + 0.163*pow(Eo, 0.757));

        // update bubble diameter 
        dbi = 1.5*cbi*pow(E, -2.0/3.0);

        // compute relative error
        relError = fabs((dbi - dbi_prev)/dbi);
    } while (relError > relTol);

    return dbi;
}

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
)
{
    // bubble Reynolds number
    double Reb = (rhol*usi*(0.001*dbi))/mul;

    // bubble Eotvos number
    double Eo = ((rhol - rhog)*g*pow(0.001*dbi,2.0))/sigma;

    return 24.0/Reb * (1.0 + pow(Reb*Eo, 0.47))*pow(1.0 - egi, 2.46*(1.0 - exp(-0.13*p/Eo)));
}

double bubbleSlipVelocity
(
    double dbi, 
    double egi, 
    double rhol, 
    double rhog, 
    double sigma, 
    double mul, 
    double p,
    double relTol
)
{
    double usi, usi_prev, Cdi, Eo, E, relError;

    // initialize bubble slip velocity
    usi = 0.001;

    do
    {
        // assign current slip velocity to previous value
        usi_prev = usi;

        // bubble drag coefficient
        Cdi = MachDrag(dbi, egi, usi_prev, rhol, rhog, sigma, mul, p);

        // bubble Eotvos number
        Eo = ((rhol - rhog)*g*pow(0.001*dbi, 2.0))/sigma;

        // bubble aspect ratio
        E = 1.0/(1.0 + 0.163*pow(Eo, 0.757));

        // update slip velocity
        usi = sqrt((4.0/3.0)*(g*0.001*dbi)/Cdi*pow(E,2.0/3.0)*(rhol - rhog)/rhol);

        // relative error on slip velocity
        relError = fabs((usi - usi_prev)/usi_prev);
    } while (relError > relTol);

    return usi;
}

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
)
{
    // compute bubble slip velocity 
    double usi = bubbleSlipVelocity(dbi, egi, rhol, rhog, sigma, mul, p, relTol);

    return (Ugbed/fmax(egi, 1e-8) - Ulbed/(1.0 - fmin(egi, 1.0 - 1e-8)) - usi);
}

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
)
{
    double egi_low, egi_high, egi, objFuncLow, objFuncMid;

    // initialize lower end of solution bracket
    egi_low = 0.0;

    // initialize upper end of solution bracket
    egi_high = 1.0;

    while (fabs(egi_high - egi_low) > absTol)
    {
        // mid-point of the solution bracket
        egi = (egi_high + egi_low)/2.0;

        // evaluate objective function at the mid-point of the solution bracket
        objFuncMid = objectiveFunction(egi, dbi, rhol, rhog, sigma, mul, p, Ugbed, Ulbed, relTol);

        // check if objective function evaluated at mid-point value returned zero
        if (objFuncMid == 0.0)
        {
            // solution has converged -- break from the loop
            break;
        }
        
        // evaluate objective function at the lower end of the solution bracket
        objFuncLow = objectiveFunction(egi_low, dbi, rhol, rhog, sigma, mul, p, Ugbed, Ulbed, relTol);

        
        // update lower and upper ends of the solution bracket
        if (objFuncLow*objFuncMid < 0.0)
        {
            // lower and mid-point evaluations have opposite signs
            // solution must be between lower end and mid-point
            // assign current mid-point value to upper limit of the bracket
            egi_high = egi;
        }
        else if (objFuncLow*objFuncMid > 0.0)
        {
            // lower and mid-point evaluations have the same sign
            // solution must be between mid-point and upper limit of the bracket
            // assign current mid-point value to lower limit of the bracket
            egi_low = egi;
        }
    }
    return egi;
}

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
)
{
    // populate bubble size and velocity distribution array
    double cbi, dbi, egi, usi, pdf;
    for (size_t i = 0; i < nClasses; i++)
    {
        // bubble chord length
        cbi = cbiMin + (double)(i*classWidth);

        // convert bubble chord length to volume equivalent diameter
        dbi = chordToDiameter(cbi, rhol, rhog, sigma, relTol);
        bsvd[i][0] = dbi;

        // calculate class holdup
        egi = bisect(dbi, rhol, rhog, sigma, mul, p, Ugbed, Ulbed, absTol, relTol);

        // calculate bubble slip velocity
        usi = bubbleSlipVelocity(dbi, egi, rhol, rhog, sigma, mul, p, relTol);
        bsvd[i][1] = usi;

        // calculate class probability density
        pdf = lognormalPDF(cbi, mu, s);
        bsvd[i][2] = pdf;
    }
}

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
)
{
    // compute lognormal distribution parameters
    double mu, s;
    lognormalDistributionParameters(rhol, rhog, mul, Ugbed, Ulbed, dor, nor, nr, Dc, drec, &mu, &s);

    // compute lower and upper bounds of bubble size distribution
    double cbiMin = lognormalQuantile(minCutOff, mu, s);
    double cbiMax = lognormalQuantile(maxCutOff, mu, s);

    // number of bubble classes
    size_t nClasses = (size_t)(round((cbiMax - cbiMin)/classWidth + 1));

    // allocate memory for bubble size and velocity distribution array
    double **bsvd = malloc(sizeof(double *)*nClasses);
    assert(bsvd != NULL);

    for (size_t i = 0; i < nClasses; i++)
    {
        bsvd[i] = malloc(3*sizeof(double));
        assert(bsvd[i] != NULL);
    }
    
    // Generate bubble size and velocity distribution
    bubbleSizeAndVelocityDistribution(rhol, rhog, sigma, mul, p, Ugbed, Ulbed, cbiMin, classWidth, nClasses, absTol, relTol, mu, s, bsvd);

    // compute gas-liquid separation efficiency
    double Reb, Eo, Bi, kprime, etai, pdfr, sumpdfr, sumpdfb;
    sumpdfr = 0.0;
    sumpdfb = 0.0;

    double Ac = 0.25*pi*(pow(Dc,2.0) - pow(drec, 2.0));
    for (size_t i = 0; i < nClasses; i++)
    {
        // bubble Reynolds number
        Reb = (rhol*bsvd[i][1]*(0.001*bsvd[i][0]))/mul;

        // bubble Eotvos number
        Eo = ((rhol - rhog)*g*pow(0.001*bsvd[i][0], 2.0))/sigma;

        // constant
        Bi = sqrt(Eo)/(beta2*Reb*tanh(bsvd[i][0]));

        // liquid residence time through the recycle pan
        kprime = vsep/(R*Ulbed*Ac);

        // individual gas-liquid separation efficiency
        etai = beta1*log(kprime) + Bi;

        // probability of entraining a bubble class in the liquid recycle
        pdfr = (1.0 - fmax(fmin(etai, 1.0), 0.0))*R*bsvd[i][2];

        // accumulate bed and recycle probabilities
        sumpdfr = sumpdfr + pdfr;
        sumpdfb = sumpdfb + bsvd[i][2];
        
        // free memory that is no longer required
        free(bsvd[i]);
    }
    
    // free remaining non-required memory
    free(bsvd);

    return 1.0 - sumpdfr/(R*sumpdfb);
}

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
)
{
    // column's cross-sectional area
    double Ac = 0.25*pi*(pow(Dc, 2.0) - pow(drec, 2.0));

    // gas recycle flow rate
    *qgr = R*(1.0 - eta)*Ugbed*Ac;

    // liquid recycle flow rate
    *qlr = R*Ulbed*Ac;
}

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
)
{
    // column's cross-sectional area
    double Ac = 0.25*pi*(pow(Dc, 2.0) - pow(drec, 2.0));

    // bed superficial gas velocity
    *Ugbed = (qgt + qgr)/Ac;

    // bed superficial liquid velocity
    *Ulbed = (qlvr + qlr)/Ac;
}

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
)
{
    // compute lognormal distribution parameters
    double mu, s;
    lognormalDistributionParameters(rhol, rhog, mul, Ugbed, Ulbed, dor, nor, nr, Dc, drec, &mu, &s);

    // compute lower and upper bounds of bubble size distribution
    double cbiMin = lognormalQuantile(minCutOff, mu, s);
    double cbiMax = lognormalQuantile(maxCutOff, mu, s);

    // number of bubble classes
    size_t nClasses = (size_t)(round((cbiMax - cbiMin)/classWidth + 1));

    // allocate memory for bubble size and velocity distribution array
    double **bsvd = malloc(sizeof(double *)*nClasses);
    assert(bsvd != NULL);

    for (size_t i = 0; i < nClasses; i++)
    {
        bsvd[i] = malloc(3*sizeof(double));
        assert(bsvd[i] != NULL);
    }
    
    // Generate bubble size and velocity distribution
    bubbleSizeAndVelocityDistribution(rhol, rhog, sigma, mul, p, Ugbed, Ulbed, cbiMin, classWidth, nClasses, absTol, relTol, mu, s, bsvd);

    // compute freeboard gas holdup
    double egi, egnum, egden;
    egnum = 0.0;
    egden = 0.0;
    for (size_t i = 0; i < nClasses; i++)
    {
        // calculate class holdup
        egi = bisect(bsvd[i][0], rhol, rhog, sigma, mul, p, Ugbed, Ulbed, absTol, relTol);

        // acummulate numerator and denominator
        egnum = egnum + egi*bsvd[i][2];
        egden = egden + bsvd[i][2];

        // free memory that is no longer required
        free(bsvd[i]);
    }
    // free remaining non-required memory
    free(bsvd);

    return egnum/egden;
}

void readGasProperties(gasPhase *properties, char *input_dir)
{
    clock_t start = clock();

    // input file
    char fname[100];
    sprintf(fname, "%s/gas.yaml", input_dir);

    FILE* fh = fopen(fname, "r");

    yaml_parser_t   parser;
    yaml_token_t    token;

    char *stringToDPointer;

    if(!yaml_parser_initialize(&parser))
    {
        fputs("Failed to initialize parser!\n",stderr);
    }

    if(fh == NULL)
    {
        fputs("Failed to open file!\n",stderr);
    }

    // set input file
    yaml_parser_set_input_file(&parser, fh);

    printf("GAS PHASE PROPERTIES\n");
    printf("--------------------\n");
    do
    {
        yaml_parser_scan(&parser, &token);

        if (token.type == YAML_SCALAR_TOKEN)
        {
            const char * key = (const char *) token.data.scalar.value;
            
            if (!strcmp(key, (const char *) "qgtMin"))
            {
                // read token state
                yaml_parser_scan(&parser, &token);

                // read token value
                yaml_parser_scan(&parser, &token);

                // store value
                properties->qgtMin = strtod((const char *) token.data.scalar.value, &stringToDPointer);

                printf("minimum gas flow rate [ m^3/s ]: %lf\n",properties->qgtMin);
            }
            else if (!strcmp(key, (const char *) "qgtMax"))
            {
                yaml_parser_scan(&parser, &token);
                yaml_parser_scan(&parser, &token);

                properties->qgtMax = strtod((const char *) token.data.scalar.value, &stringToDPointer);
                printf("maximum gas flow rate [ m^3/s ]: %lf\n", properties->qgtMax);
            }
            else if (!strcmp(key, (const char *) "qgtIncrement"))
            {
                yaml_parser_scan(&parser, &token);
                yaml_parser_scan(&parser, &token);

                properties->qgtIncrement = strtod((const char *) token.data.scalar.value, &stringToDPointer);
                printf(" gas flow rate  increment [ m^3/s ]: %lf\n",properties->qgtIncrement);
            }
            else if (!strcmp(key, (const char *) "rhog"))
            {
                yaml_parser_scan(&parser, &token);
                yaml_parser_scan(&parser, &token);

                properties->rhog = strtod((const char *) token.data.scalar.value, &stringToDPointer);
                printf("gas density [ kg/m^3 ]: %lf\n",properties->rhog);
            }
        }
        
    } while (token.type != YAML_STREAM_END_TOKEN);
    
    clock_t end = clock();

    double duration = (double) (end - start);
    printf("execution time [ seconds ]: %lf\n\n", duration);

    yaml_token_delete(&token);
    yaml_parser_delete(&parser);
    fclose(fh);
}

void readLiquidProperties(liquidPhase *properties, char *input_dir)
{

    clock_t start = clock();

    char fname[100];
    sprintf(fname, "%s/liquid.yaml", input_dir);

    FILE* fh = fopen(fname, "r");

    yaml_parser_t   parser;
    yaml_token_t    token;

    char *stringToDPointer;

    if(!yaml_parser_initialize(&parser))
    {
        fputs("failed to initialize parser!\n",stderr);
    }

    if(fh == NULL)
    {
        fputs("failed to open file!\n",stderr);
    }

    // set input file
    yaml_parser_set_input_file(&parser, fh);

    printf("LIQUID PHASE PROPERTIES\n");
    printf("-----------------------\n");
    do
    {
        yaml_parser_scan(&parser, &token);

        if (token.type == YAML_SCALAR_TOKEN)
        {
            const char * key = (const char *) token.data.scalar.value;
            
            if (!strcmp(key, (const char *) "qlvr"))
            {
                yaml_parser_scan(&parser, &token);
                yaml_parser_scan(&parser, &token);

                properties->qlvr = strtod((const char *) token.data.scalar.value, &stringToDPointer);

                printf("liquid flow rate [ m^3/s ]: %lf\n",properties->qlvr);
            }
            else if (!strcmp(key, (const char *) "rhol"))
            {
                yaml_parser_scan(&parser, &token);
                yaml_parser_scan(&parser, &token);

                properties->rhol = strtod((const char *) token.data.scalar.value, &stringToDPointer);
                printf("liquid density [ kg/m^3 ]: %lf\n",properties->rhol);
            }
            else if (!strcmp(key, (const char *) "mul"))
            {
                yaml_parser_scan(&parser, &token);
                yaml_parser_scan(&parser, &token);

                properties->mul = strtod((const char *) token.data.scalar.value, &stringToDPointer);
                printf("liquid viscosity [ Pa.s ]: %lf\n",properties->mul);
            }
        }
        
    } while (token.type != YAML_STREAM_END_TOKEN);
    
    clock_t end = clock();

    double duration = (double) (end - start);
    printf("execution time [ seconds ]: %lf\n\n", duration);
    yaml_token_delete(&token);
    yaml_parser_delete(&parser);
    fclose(fh);
}

void readSolidProperties(solidPhase *properties, char *input_dir)
{
    
    clock_t start = clock();

    char fname[100];
    sprintf(fname, "%s/solids.yaml", input_dir);

    FILE* fh = fopen(fname, "r");

    yaml_parser_t   parser;
    yaml_token_t    token;

    char *stringToDPointer;

    if(!yaml_parser_initialize(&parser))
    {
        fputs("Failed to initialize parser!\n",stderr);
    }

    if(fh == NULL)
    {
        fputs("Failed to open file!\n",stderr);
    }

    // set input file
    yaml_parser_set_input_file(&parser, fh);

    printf("SOLID PHASE PROPERTIES\n");
    printf("-----------------------\n");
    do
    {
        yaml_parser_scan(&parser, &token);

        if (token.type == YAML_SCALAR_TOKEN)
        {
            const char * key = (const char *) token.data.scalar.value;
            
            if (!strcmp(key, (const char *) "mcat"))
            {
                yaml_parser_scan(&parser, &token);
                yaml_parser_scan(&parser, &token);

                properties->mcat = strtod((const char *) token.data.scalar.value, &stringToDPointer);

                printf("mass of catalyst inventory [ kg ]: %lf\n",properties->mcat);
            }
            else if (!strcmp(key, (const char *) "Vp"))
            {
                yaml_parser_scan(&parser, &token);
                yaml_parser_scan(&parser, &token);

                properties->Vp = strtod((const char *) token.data.scalar.value, &stringToDPointer);
                printf("particle volume [ m^3 ]: %lf\n",properties->Vp);
            }
            else if (!strcmp(key, (const char *) "Ap"))
            {
                yaml_parser_scan(&parser, &token);
                yaml_parser_scan(&parser, &token);

                properties->Ap = strtod((const char *) token.data.scalar.value, &stringToDPointer);
                printf("particle surface area [ m^2 ]: %lf\n",properties->Ap);
            }
            else if (!strcmp(key, (const char *) "rhop"))
            {
                yaml_parser_scan(&parser, &token);
                yaml_parser_scan(&parser, &token);

                properties->rhop = strtod((const char *) token.data.scalar.value, &stringToDPointer);
                printf("particle density [ kg/m^3 ]: %lf\n",properties->rhop);
            }
        }
        
    } while (token.type != YAML_STREAM_END_TOKEN);
    
    clock_t end = clock();

    double duration = (double) (end - start);
    printf("execution time [ seconds ]: %lf\n\n", duration);
    yaml_token_delete(&token);
    yaml_parser_delete(&parser);
    fclose(fh);
}

void readInterfacialProperties(interfacial *properties, char *input_dir)
{
    
    clock_t start = clock();

    char fname[100];
    sprintf(fname, "%s/interfacial.yaml", input_dir);

    FILE* fh = fopen(fname, "r");

    yaml_parser_t   parser;
    yaml_token_t    token;

    char *stringToDPointer;

    if(!yaml_parser_initialize(&parser))
    {
        fputs("Failed to initialize parser!\n",stderr);
    }

    if(fh == NULL)
    {
        fputs("Failed to open file!\n",stderr);
    }

    // set input file
    yaml_parser_set_input_file(&parser, fh);

    printf("INTERFACIAL PROPERTIES\n");
    printf("----------------------\n");
    do
    {
        yaml_parser_scan(&parser, &token);

        if (token.type == YAML_SCALAR_TOKEN)
        {
            const char * key = (const char *) token.data.scalar.value;
            
            if (!strcmp(key, (const char *) "sigma"))
            {
                yaml_parser_scan(&parser, &token);
                yaml_parser_scan(&parser, &token);

                properties->sigma = strtod((const char *) token.data.scalar.value, &stringToDPointer);

                printf("gas-liquid interfacial tension [ N/m ]: %lf\n",properties->sigma);
            }
            
        }
        
    } while (token.type != YAML_STREAM_END_TOKEN);
    
    clock_t end = clock();

    double duration = (double) (end - start);
    printf("execution time [ seconds ]: %lf\n\n", duration);
    yaml_token_delete(&token);
    yaml_parser_delete(&parser);
    fclose(fh);
}

void readColumnParameters(column *parameters, char *input_dir)
{
    
    clock_t start = clock();

    char fname[100];
    sprintf(fname, "%s/column.yaml", input_dir);

    FILE* fh = fopen(fname, "r");

    yaml_parser_t   parser;
    yaml_token_t    token;

    char *stringToDPointer;

    if(!yaml_parser_initialize(&parser))
    {
        fputs("Failed to initialize parser!\n",stderr);
    }

    if(fh == NULL)
    {
        fputs("Failed to open file!\n",stderr);
    }

    // set input file
    yaml_parser_set_input_file(&parser, fh);

    printf("COLUMN PARAMETERS\n");
    printf("-----------------\n");
    do
    {
        yaml_parser_scan(&parser, &token);

        if (token.type == YAML_SCALAR_TOKEN)
        {
            const char * key = (const char *) token.data.scalar.value;
            
            if (!strcmp(key, (const char *) "P"))
            {
                yaml_parser_scan(&parser, &token);
                yaml_parser_scan(&parser, &token);

                parameters->P = strtod((const char *) token.data.scalar.value, &stringToDPointer);

                printf("operating pressure [ MPa ]: %lf\n",parameters->P);
            }
            else if (!strcmp(key, (const char *) "Dc"))
            {
                yaml_parser_scan(&parser, &token);
                yaml_parser_scan(&parser, &token);

                parameters->Dc = strtod((const char *) token.data.scalar.value, &stringToDPointer);
                printf("column diameter [ m ]: %lf\n",parameters->Dc);
            }
            else if (!strcmp(key, (const char *) "drec"))
            {
                yaml_parser_scan(&parser, &token);
                yaml_parser_scan(&parser, &token);

                parameters->drec = strtod((const char *) token.data.scalar.value, &stringToDPointer);
                printf("recycle line diameter [ m^2 ]: %lf\n",parameters->drec);
            }
            else if (!strcmp(key, (const char *) "vsep"))
            {
                yaml_parser_scan(&parser, &token);
                yaml_parser_scan(&parser, &token);

                parameters->vsep = strtod((const char *) token.data.scalar.value, &stringToDPointer);
                printf("separator volume [ m^3 ]: %lf\n",parameters->vsep);
            }
            else if (!strcmp(key, (const char *) "hset"))
            {
                yaml_parser_scan(&parser, &token);
                yaml_parser_scan(&parser, &token);

                parameters->hset = strtod((const char *) token.data.scalar.value, &stringToDPointer);
                printf("desired bed height [ m ]: %lf\n",parameters->hset);
            }
        }
        
    } while (token.type != YAML_STREAM_END_TOKEN);
    
    clock_t end = clock();

    double duration = (double) (end - start);
    printf("execution time [ seconds ]: %lf\n\n", duration);
    yaml_token_delete(&token);
    yaml_parser_delete(&parser);
    fclose(fh);
}

void readGridParameters(grid *parameters, char *input_dir)
{
    
    clock_t start = clock();

    char fname[100];
    sprintf(fname, "%s/grid.yaml", input_dir);

    FILE* fh = fopen(fname, "r");

    yaml_parser_t   parser;
    yaml_token_t    token;

    char *stringToDPointer;

    if(!yaml_parser_initialize(&parser))
    {
        fputs("Failed to initialize parser!\n",stderr);
    }

    if(fh == NULL)
    {
        fputs("Failed to open file!\n",stderr);
    }

    // set input file
    yaml_parser_set_input_file(&parser, fh);

    printf("GRID PARAMETERS\n");
    printf("---------------\n");
    do
    {
        yaml_parser_scan(&parser, &token);

        if (token.type == YAML_SCALAR_TOKEN)
        {
            const char * key = (const char *) token.data.scalar.value;
            
            if (!strcmp(key, (const char *) "nr"))
            {
                yaml_parser_scan(&parser, &token);
                yaml_parser_scan(&parser, &token);

                parameters->nr = strtod((const char *) token.data.scalar.value, &stringToDPointer);

                printf("number of risers [ - ]: %lf\n",parameters->nr);
            }
            else if (!strcmp(key, (const char *) "nor"))
            {
                yaml_parser_scan(&parser, &token);
                yaml_parser_scan(&parser, &token);

                parameters->nor = strtod((const char *) token.data.scalar.value, &stringToDPointer);
                printf("number of outlet orifices per riser [ - ]: %lf\n",parameters->nor);
            }
            else if (!strcmp(key, (const char *) "dor"))
            {
                yaml_parser_scan(&parser, &token);
                yaml_parser_scan(&parser, &token);

                parameters->dor = strtod((const char *) token.data.scalar.value, &stringToDPointer);
                printf("outlet orifice diameter [ m ]: %lf\n",parameters->dor);
            }
        }
        
    } while (token.type != YAML_STREAM_END_TOKEN);
    
    clock_t end = clock();

    double duration = (double) (end - start);
    printf("execution time [ seconds ]: %lf\n\n", duration);
    yaml_token_delete(&token);
    yaml_parser_delete(&parser);
    fclose(fh);
}

void readBubbleSizeDistributionParameters(bubbleSizeDistribution *parameters, char *input_dir)
{
    
    clock_t start = clock();

    char fname[100];
    sprintf(fname, "%s/bubbleSizeDistribution.yaml", input_dir);

    FILE* fh = fopen(fname, "r");

    yaml_parser_t   parser;
    yaml_token_t    token;

    char *stringToDPointer;

    if(!yaml_parser_initialize(&parser))
    {
        fputs("Failed to initialize parser!\n",stderr);
    }

    if(fh == NULL)
    {
        fputs("Failed to open file!\n",stderr);
    }

    // set input file
    yaml_parser_set_input_file(&parser, fh);

    printf("BUBBLE SIZE DISTRIBUTION PARAMETERS\n");
    printf("-----------------------------------\n");
    do
    {
        yaml_parser_scan(&parser, &token);

        if (token.type == YAML_SCALAR_TOKEN)
        {
            const char * key = (const char *) token.data.scalar.value;
            
            if (!strcmp(key, (const char *) "classWidth"))
            {
                yaml_parser_scan(&parser, &token);
                yaml_parser_scan(&parser, &token);

                parameters->classWidth = strtod((const char *) token.data.scalar.value, &stringToDPointer);

                printf("bubble class width [ mm ]: %lf\n",parameters->classWidth);
            }
            else if (!strcmp(key, (const char *) "minCutOff"))
            {
                yaml_parser_scan(&parser, &token);
                yaml_parser_scan(&parser, &token);

                parameters->minCutOff = strtod((const char *) token.data.scalar.value, &stringToDPointer);
                printf("cumulative probability for determining minimum bubble size [ - ]: %lf\n",parameters->minCutOff);
            }
            else if (!strcmp(key, (const char *) "maxCutOff"))
            {
                yaml_parser_scan(&parser, &token);
                yaml_parser_scan(&parser, &token);

                parameters->maxCutOff = strtod((const char *) token.data.scalar.value, &stringToDPointer);
                printf("cumulative probability for determining maximum bubble size [ - ]: %lf\n",parameters->maxCutOff);
            }
        }
        
    } while (token.type != YAML_STREAM_END_TOKEN);
    
    clock_t end = clock();

    double duration = (double) (end - start);
    printf("execution time [ seconds ]: %lf\n\n", duration);
    yaml_token_delete(&token);
    yaml_parser_delete(&parser);
    fclose(fh);
}

void readSeparatorParameters(separator *parameters, char *input_dir)
{
    
    clock_t start = clock();

    char fname[100];
    sprintf(fname, "%s/separator.yaml", input_dir);

    FILE* fh = fopen(fname, "r");

    yaml_parser_t   parser;
    yaml_token_t    token;

    char *stringToDPointer;

    if(!yaml_parser_initialize(&parser))
    {
        fputs("Failed to initialize parser!\n",stderr);
    }

    if(fh == NULL)
    {
        fputs("Failed to open file!\n",stderr);
    }

    // set input file
    yaml_parser_set_input_file(&parser, fh);

    printf("SEPARATOR PARAMETERS\n");
    printf("--------------------\n");
    do
    {
        yaml_parser_scan(&parser, &token);

        if (token.type == YAML_SCALAR_TOKEN)
        {
            const char * key = (const char *) token.data.scalar.value;
            
            if (!strcmp(key, (const char *) "beta1"))
            {
                yaml_parser_scan(&parser, &token);
                yaml_parser_scan(&parser, &token);

                parameters->beta1 = strtod((const char *) token.data.scalar.value, &stringToDPointer);

                printf("beta1 [ 1/s ]: %lf\n",parameters->beta1);
            }
            else if (!strcmp(key, (const char *) "beta2"))
            {
                yaml_parser_scan(&parser, &token);
                yaml_parser_scan(&parser, &token);

                parameters->beta2 = strtod((const char *) token.data.scalar.value, &stringToDPointer);
                printf("beta2 [ 1/mm ]: %lf\n",parameters->beta2);
            }
            else if (!strcmp(key, (const char *) "Rstep"))
            {
                yaml_parser_scan(&parser, &token);
                yaml_parser_scan(&parser, &token);

                parameters->Rstep = strtod((const char *) token.data.scalar.value, &stringToDPointer);
                printf("liquid recycle ratio step size [ - ]: %lf\n",parameters->Rstep);
            }
        }
        
    } while (token.type != YAML_STREAM_END_TOKEN);
    
    clock_t end = clock();

    double duration = (double) (end - start);
    printf("execution time [ seconds ]: %lf\n\n", duration);
    yaml_token_delete(&token);
    yaml_parser_delete(&parser);
    fclose(fh);
}

void readConvergenceCriteria(convergence *criteria, char *input_dir)
{
    
    clock_t start = clock();

    char fname[100];
    sprintf(fname, "%s/convergence.yaml", input_dir);

    FILE* fh = fopen(fname, "r");

    yaml_parser_t   parser;
    yaml_token_t    token;

    char *stringToDPointer;

    if(!yaml_parser_initialize(&parser))
    {
        fputs("Failed to initialize parser!\n",stderr);
    }

    if(fh == NULL)
    {
        fputs("Failed to open file!\n",stderr);
    }

    // set input file
    yaml_parser_set_input_file(&parser, fh);

    printf("CONVERGENCE CRITERIA\n");
    printf("--------------------\n");
    do
    {
        yaml_parser_scan(&parser, &token);

        if (token.type == YAML_SCALAR_TOKEN)
        {
            const char * key = (const char *) token.data.scalar.value;
            
            if (!strcmp(key, (const char *) "absTol"))
            {
                yaml_parser_scan(&parser, &token);
                yaml_parser_scan(&parser, &token);

                criteria->absTol = strtod((const char *) token.data.scalar.value, &stringToDPointer);

                printf("absolute tolerance [ dependent ]: %lf\n",criteria->absTol);
            }
            else if (!strcmp(key, (const char *) "relTol"))
            {
                yaml_parser_scan(&parser, &token);
                yaml_parser_scan(&parser, &token);

                criteria->relTol = strtod((const char *) token.data.scalar.value, &stringToDPointer);
                printf("relative tolerance [ - ]: %lf\n",criteria->relTol);
            }
            else if (!strcmp(key, (const char *) "maxIter"))
            {
                yaml_parser_scan(&parser, &token);
                yaml_parser_scan(&parser, &token);

                criteria->maxIter = strtod((const char *) token.data.scalar.value, &stringToDPointer);
                printf("maximum iterations [ - ]: %ld\n",criteria->maxIter);
            }
        }
        
    } while (token.type != YAML_STREAM_END_TOKEN);
    
    clock_t end = clock();

    double duration = (double) (end - start);
    printf("execution time [ seconds ]: %lf\n\n", duration);
    yaml_token_delete(&token);
    yaml_parser_delete(&parser);
    fclose(fh);
}

void solve
(
    gasPhase gasProperties,
    liquidPhase liquidProperties,
    solidPhase  solidProperties,
    interfacial interfacialProperties,
    column  columnParameters,
    grid    gridParameters,
    bubbleSizeDistribution  bubbleSizeDistributionParameters,
    separator   separatorParameters,
    convergence convergenceCriteria,
    char *output_dir
)
{
    // create a file for storing summary results and write column headings to it.
    char summary_filename[100];
    sprintf(summary_filename, "%s/summary_results.dat", output_dir);
        
    FILE* summary_file;
    summary_file = fopen(summary_filename, "w");

    fprintf(summary_file, "Ugt[m/s]\tUlvr[m/s]\tUgbed[m/s]\tUlbed[m/s]\tUmf[m/s]\tegbed\tepsilonl\tegfb\teta\tR\thbed[m]\tRgas\tmu\ts\n");

    // simulation and output parameters
    double qgt, Ugt, Ulvr, Ugbed, Ulbed, Umf, eta, qgr, qlr, epsilonp, hbed, egfb, egbed, epsilonl, R, Rgas, mu, s, error;
    double dv, phi, Arp, emfo, Umfo, k, n, upinf;

    // number of inlet gas flow rates to be simulated
    size_t nGasFlows = (size_t) round((gasProperties.qgtMax - gasProperties.qgtMin)/gasProperties.qgtIncrement + 1.0);

    // dimensionless operating pressure 
    double p = columnParameters.P/0.10;

    for (size_t i = 0; i < nGasFlows; i++)
    {
        clock_t start = clock();

        // initialize/reset recycle gas and liquid flow rates
        qgr = 0.0;
        qlr = 0.0;

        // initialize/reset liquid recycle ratio
        R = 0.0;

        // while-loop counter
        size_t counter = 0;

        // gas flow rate to be simulated
        qgt = gasProperties.qgtMin + (double)(i)* (gasProperties.qgtIncrement);

        printf("Inlet gas flow rate: %.2lf m\u00B3/s\n",qgt);
        printf("-------------------------------------------\n");
        do
        {
            // update bed flow rates
            bedFlowRates
            (
                qgt, 
                liquidProperties.qlvr, 
                columnParameters.Dc, 
                columnParameters.drec, 
                qgr, 
                qlr, 
                &Ugbed, 
                &Ulbed
            );

            // calculate minimum fluidization velocity for a two phase system
            // these could be taken out of the do-while loop, but remains here to make it easy to follow calculation steps
            dv = volumeEquivalentDiameter(solidProperties.Vp); // volume equivalent particle diameter
            phi = sphericity(solidProperties.Vp, solidProperties.Ap); // particle sphericity
            Arp = ArchimedesNumber(solidProperties.rhop, liquidProperties.rhol, liquidProperties.mul, dv); // particle-liquid Archimedes number
            emfo = minimumVoidFraction(phi);    // void fraction at minimum fluidization (two phase system)
            Umfo = LSMinimumFluidizationVelocity(liquidProperties.rhol, liquidProperties.mul, dv, phi, Arp, emfo); // minimum fluidization velocity (two phase system)

            // calculate minimum fluidization velocity for a three phase system
            Umf = GLSMinimumFluidizationVelocity(solidProperties.rhop, liquidProperties.rhol, gasProperties.rhog, liquidProperties.mul, Ugbed, dv, phi, emfo, Umfo, convergenceCriteria.relTol);

            // ensure bed liquid velocity is greater than minimum fluidization velocity for the three phase system
            while (Ulbed < Umf)
            {
                // increase liquid recycle ratio since bed liquid velocity is still less than minimum fluidization velocity
                R = fmax(fmin(R + separatorParameters.Rstep, 1.0), 0.0);

                // calculate/update gas-liquid separation efficiency of the recycle pan
                eta = GLSeparator
                (
                    liquidProperties.rhol, 
                    gasProperties.rhog, 
                    liquidProperties.mul, 
                    interfacialProperties.sigma, 
                    Ugbed, 
                    Ulbed, 
                    gridParameters.dor, 
                    gridParameters.nor, 
                    gridParameters.nr, 
                    columnParameters.Dc, 
                    columnParameters.drec, 
                    R, 
                    columnParameters.vsep, 
                    p, 
                    bubbleSizeDistributionParameters.classWidth, 
                    separatorParameters.beta1, 
                    separatorParameters.beta2, 
                    bubbleSizeDistributionParameters.minCutOff, 
                    bubbleSizeDistributionParameters.maxCutOff, 
                    convergenceCriteria.absTol, 
                    convergenceCriteria.relTol
                );
                
                // update recycle flow rates
                recycleFlowRates(Ugbed, Ulbed, columnParameters.Dc, columnParameters.drec, R, eta, &qgr, &qlr);

                // update bed flow rate
                bedFlowRates(qgt, liquidProperties.qlvr, columnParameters.Dc, columnParameters.drec, qgr, qlr, &Ugbed, &Ulbed);

                // update minimum fluidization velocity for the three phase system 
                // no need to update two phase value since none of its input parameters has changed
                Umf = GLSMinimumFluidizationVelocity(solidProperties.rhop, liquidProperties.rhol, gasProperties.rhog, liquidProperties.mul, Ugbed, dv, phi, emfo, Umfo, convergenceCriteria.relTol);
            }

            // calculate particles holdup
            // k, n & upinf could also be taken out of the do-while loop
            k = wallEffectParameter(dv, columnParameters.Dc); // wall effect parameter
            n = RichardsonZakiCoefficient(Arp, dv, columnParameters.Dc);  // Richardson-Zaki coefficient
            upinf = terminalVelocity(liquidProperties.rhol, liquidProperties.mul, Arp, dv, phi);    // particle terminal velocity
            epsilonp = particlesHoldup(Ugbed, Ulbed, k, n, upinf);  // particles holdup
            
            // calculate bed height
            hbed = bedHeight(columnParameters.drec, columnParameters.Dc, epsilonp, solidProperties.mcat, solidProperties.rhop);

            // deviation of calculated bed height from desired value
            error = (hbed - columnParameters.hset);

            // break from loop if absolute error on bed height is smaller than absolute tolerance
            if(fabs(error) < convergenceCriteria.absTol)
            {
                printf("SOLUTION CONVERGED.\n");
                break;
            }

            // otherwise update liquid recycle ratio appropriately
            R = fmax(fmin(R - separatorParameters.Rstep*error, 1.0), 0.0);
            
            // update gas-liquid separation efficiency of the recycle pan
            eta = GLSeparator
            (
                liquidProperties.rhol, 
                gasProperties.rhog, 
                liquidProperties.mul, 
                interfacialProperties.sigma, 
                Ugbed, 
                Ulbed, 
                gridParameters.dor, 
                gridParameters.nor, 
                gridParameters.nr, 
                columnParameters.Dc, 
                columnParameters.drec, 
                R, 
                columnParameters.vsep, 
                p, 
                bubbleSizeDistributionParameters.classWidth, 
                separatorParameters.beta1, 
                separatorParameters.beta2, 
                bubbleSizeDistributionParameters.minCutOff, 
                bubbleSizeDistributionParameters.maxCutOff, 
                convergenceCriteria.absTol, 
                convergenceCriteria.relTol
            );
                
            // update recycle flow rates
            recycleFlowRates(Ugbed, Ulbed, columnParameters.Dc, columnParameters.drec, R, eta, &qgr, &qlr);

            // display simulation progress
            printf("iteration: %ld, current bed height [m]: %.4lf, deviation from desired bed height [m]: %.4lf\n", counter, hbed, error);
            
            // update do-while loop counter
            counter++;
        } while (counter < convergenceCriteria.maxIter);

        // calculate freeboard gas holdup
        egfb = freeboardGasHoldup(liquidProperties.rhol, gasProperties.rhog, interfacialProperties.sigma, liquidProperties.mul,p, Ugbed, Ulbed, gridParameters.dor, 
        gridParameters.nor, gridParameters.nr, columnParameters.Dc, columnParameters.drec, convergenceCriteria.absTol, convergenceCriteria.relTol, bubbleSizeDistributionParameters.minCutOff, bubbleSizeDistributionParameters.maxCutOff, bubbleSizeDistributionParameters.classWidth);
        
        // calculate bed gas holdup
        egbed = egfb/1.3;

        // calculate bed liquid holdup
        epsilonl = 1.0 - egbed - epsilonp;

        // calculate recycle gas fraction
        Rgas = R*(1.0 - eta);

        // calculate inlet gas velocity
        Ugt = qgt/(0.25*pi*(pow(columnParameters.Dc, 2.0) - pow(columnParameters.drec, 2.0)));

        // calculate inlet gas velocity
        Ulvr = liquidProperties.qlvr/(0.25*pi*(pow(columnParameters.Dc, 2.0) - pow(columnParameters.drec, 2.0)));

        // write summary results to file
        fprintf(summary_file,"%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n",
        Ugt, Ulvr,Ugbed, Ulbed, Umf, egbed, epsilonl, egfb, eta, R, hbed, Rgas, mu, s);


        // Generate bubble size and slip velocity distribution for storage
        lognormalDistributionParameters(liquidProperties.rhol, gasProperties.rhog, liquidProperties.mul, 
        Ugbed, Ulbed, gridParameters.dor, gridParameters.nor, gridParameters.nr, columnParameters.Dc, columnParameters.drec, &mu, &s);

        double cbiMin = lognormalQuantile(bubbleSizeDistributionParameters.minCutOff, mu, s);
        double cbiMax = lognormalQuantile(bubbleSizeDistributionParameters.maxCutOff, mu, s);

        size_t nClasses = (size_t)(round((cbiMax - cbiMin)/bubbleSizeDistributionParameters.classWidth + 1.0));

        double **bsvd = malloc(sizeof(double *)*nClasses);
        assert(bsvd != NULL);

        for (size_t j = 0; j < nClasses; j++)
        {
            bsvd[j] = malloc(3*sizeof(double));
            assert(bsvd[j] != NULL);
        }
    
        bubbleSizeAndVelocityDistribution(liquidProperties.rhol, gasProperties.rhog, interfacialProperties.sigma, liquidProperties.mul, p, 
        Ugbed, Ulbed, cbiMin, bubbleSizeDistributionParameters.classWidth, nClasses, convergenceCriteria.absTol, convergenceCriteria.relTol, mu, s, bsvd);

        // write bubble size and slip velocity distribution to file
        char dist_filename[100];
        sprintf(dist_filename, "%s/QGT-%ld.dat", output_dir, (size_t)(100*qgt));
        
        FILE* dist_file;
        dist_file = fopen(dist_filename, "w");

        fprintf(dist_file, "%s\t%s\t%s\n","dbi[mm]", "usi[m/s]", "pdf");

        for (size_t k = 0; k < nClasses; k++)
        {
            fprintf(dist_file, "%lf\t%lf\t%lf\n",bsvd[k][0], bsvd[k][1], bsvd[k][2]);
        
            free(bsvd[k]);
        }
        free(bsvd);
        fclose(dist_file);

        clock_t end = clock();
        double  duration = (double) (end - start);
        printf("execution time [ seconds ]: %lf\n\n\n", duration);
    }

    // close summary results file
    fclose(summary_file);
}
