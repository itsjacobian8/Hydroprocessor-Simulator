#include<stdio.h>
#include<stdlib.h>
#include "submodels.h"

int main(int argc, char *argv[])
{

    if (argc !=3)
    {
        printf("Usage: ./simulator <input directory> <output directory>\n");
        abort();
    }
    
    // input directory
    char *inputDir = argv[1];

    // output directory
    char *outputDir = argv[2];
    
    // read in gas phase properties
    gasPhase gasProperties;
    readGasProperties(&gasProperties, inputDir);
    
    // read in liquid phase properties
    liquidPhase liquidProperties;
    readLiquidProperties(&liquidProperties, inputDir);

    // read in solid phase properties
    solidPhase  solidProperties;
    readSolidProperties(&solidProperties, inputDir);

    // read in interfacial properties
    interfacial interfacialProperties;
    readInterfacialProperties(&interfacialProperties, inputDir);

    // read in column properties
    column  columnParameters;
    readColumnParameters(&columnParameters, inputDir);

    // read in grid parameters
    grid   gridParameters;
    readGridParameters(&gridParameters, inputDir);

    // read in bubble size distribution parameters
    bubbleSizeDistribution  bubbleSizeDistributionParameters;
    readBubbleSizeDistributionParameters(&bubbleSizeDistributionParameters, inputDir);

    // read in gas-liquid separation submodel parameters
    separator   separatorParameters;
    readSeparatorParameters(&separatorParameters, inputDir);

    // read convergence criteria
    convergence convergenceCriteria;
    readConvergenceCriteria(&convergenceCriteria, inputDir);

    // call to solver
    solve
    (
        gasProperties, 
        liquidProperties, 
        solidProperties, 
        interfacialProperties, 
        columnParameters,
        gridParameters,
        bubbleSizeDistributionParameters,
        separatorParameters,
        convergenceCriteria,
        outputDir
    );

    return 0;
}

