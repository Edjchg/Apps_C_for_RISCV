// Copyright (c) 2007 Intel Corp.

// Black-Scholes
// Analytical method for calculating European Options
//
//
// Reference Source: Options, Futures, and Other Derivatives, 3rd Edition, Prentice
// Hall, John C. Hull,

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

//double max_otype, min_otype ;
//double max_sptprice, min_sptprice;
//double max_strike, min_strike;
//double max_rate, min_rate ;
//double max_volatility, min_volatility;
//double max_otime, min_otime ;
//double max_out_price, min_out_price;

#define DIVIDE 120.0

//Precision to use for calculations
#define fptype float

#define NUM_RUNS 1

typedef struct OptionData_
{
    fptype s;        // spot price
    fptype strike;   // strike price
    fptype r;        // risk-free interest rate
    fptype divq;     // dividend rate
    fptype v;        // volatility
    fptype t;        // time to maturity or option expiration in years
                     //     (1yr = 1.0, 6mos = 0.5, 3mos = 0.25, ..., etc)
    char OptionType; // Option type.  "P"=PUT, "C"=CALL
    fptype divs;     // dividend vals (not used in this test)
    fptype DGrefval; // DerivaGem Reference Value
} OptionData;

OptionData *data;
fptype *prices;
int numOptions;

int *otype;
fptype *sptprice;
fptype *strike;
fptype *rate;
fptype *volatility;
fptype *otime;
int numError = 0;

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
// Cumulative Normal Distribution Function
// See Hull, Section 11.8, P.243-244
#define inv_sqrt_2xPI 0.39894228040143270286

fptype CNDF(fptype InputX)
{
    int sign;

    fptype OutputX;
    fptype xInput;
    fptype xNPrimeofX;
    fptype expValues;
    fptype xK2;
    fptype xK2_2, xK2_3;
    fptype xK2_4, xK2_5;
    fptype xLocal, xLocal_1;
    fptype xLocal_2, xLocal_3;

    // Check for negative value of InputX
    if (InputX < 0.0)
    {
        InputX = -InputX;
        sign = 1;
    }
    else
        sign = 0;

    xInput = InputX;

    // Compute NPrimeX term common to both four & six decimal accuracy calcs
    expValues = exp(-0.5f * InputX * InputX);
    xNPrimeofX = expValues;
    xNPrimeofX = xNPrimeofX * inv_sqrt_2xPI;

    xK2 = 0.2316419 * xInput;
    xK2 = 1.0 + xK2;
    xK2 = 1.0 / xK2;
    xK2_2 = xK2 * xK2;
    xK2_3 = xK2_2 * xK2;
    xK2_4 = xK2_3 * xK2;
    xK2_5 = xK2_4 * xK2;

    xLocal_1 = xK2 * 0.319381530;
    xLocal_2 = xK2_2 * (-0.356563782);
    xLocal_3 = xK2_3 * 1.781477937;
    xLocal_2 = xLocal_2 + xLocal_3;
    xLocal_3 = xK2_4 * (-1.821255978);
    xLocal_2 = xLocal_2 + xLocal_3;
    xLocal_3 = xK2_5 * 1.330274429;
    xLocal_2 = xLocal_2 + xLocal_3;

    xLocal_1 = xLocal_2 + xLocal_1;
    xLocal = xLocal_1 * xNPrimeofX;

    //printf("# xLocal: %10.10f\n", xLocal);

    xLocal = 1.0 - xLocal;

    OutputX = xLocal;

    //printf("# Output: %10.10f\n", OutputX);

    if (sign)
    {
        OutputX = 1.0 - OutputX;
    }

    return OutputX;
}

//////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////
fptype BlkSchlsEqEuroNoDiv(fptype sptprice,
                           fptype strike, fptype rate, fptype volatility,
                           fptype time, int otype, float timet, fptype *N1, fptype *N2)
{
    fptype OptionPrice;

    // local private working variables for the calculation
    //fptype xStockPrice;
    //fptype xStrikePrice;
    fptype xRiskFreeRate;
    fptype xVolatility;
    fptype xTime;
    fptype xSqrtTime;

    fptype logValues;
    fptype xLogTerm;
    fptype xD1;
    fptype xD2;
    fptype xPowerTerm;
    fptype xDen;
    fptype d1;
    fptype d2;
    fptype FutureValueX;
    fptype NofXd1;
    fptype NofXd2;
    fptype NegNofXd1;
    fptype NegNofXd2;

    //xStockPrice = sptprice;
    //xStrikePrice = strike;
    xRiskFreeRate = rate;
    xVolatility = volatility;
    xTime = time;

    xSqrtTime = sqrt(xTime);

    logValues = log(sptprice / strike);

    xLogTerm = logValues;

    xPowerTerm = xVolatility * xVolatility;
    xPowerTerm = xPowerTerm * 0.5;

    xD1 = xRiskFreeRate + xPowerTerm;
    xD1 = xD1 * xTime;
    xD1 = xD1 + xLogTerm;

    xDen = xVolatility * xSqrtTime;
    xD1 = xD1 / xDen;
    xD2 = xD1 - xDen;

    d1 = xD1;
    d2 = xD2;

    NofXd1 = CNDF(d1);

    if (NofXd1 > 1.0)
        printf("Greater than one!");
    //printf("# d1: %10.10f\n", NofXd1);

    NofXd2 = CNDF(d2);
    if (NofXd2 > 1.0)
        printf("Greater than one!");
    //printf("# d2: %10.10f\n", NofXd2);

    *N1 = NofXd1;
    *N2 = NofXd2;

    FutureValueX = strike * (exp(-(rate) * (time)));
    if (otype == 0)
    {
        OptionPrice = (sptprice * NofXd1) - (FutureValueX * NofXd2);
    }
    else
    {
        NegNofXd1 = (1.0 - NofXd1);
        NegNofXd2 = (1.0 - NofXd2);
        OptionPrice = (FutureValueX * NegNofXd2) - (sptprice * NegNofXd1);
    }

    return OptionPrice;
}

double normalize(double in, double min, double max, double min_new, double max_new)
{
    return (((in - min) / (max - min)) * (max_new - min_new)) + min_new;
}

int bs_thread(void *tid_ptr)
{
    int i, j;

    int tid = *(int *)tid_ptr;
    int start = tid * (numOptions);
    int end = start + (numOptions);
    fptype price_orig;

    for (j = 0; j < NUM_RUNS; j++)
    {
        for (i = start; i < end; i++)
        {
            /* Calling main function to calculate option value based on 
             * Black & Scholes's equation.
             */
            fptype price;
            fptype N1, N2;

            double dataIn[6];
            double dataOut[1];

            dataIn[0] = sptprice[i];
            dataIn[1] = strike[i];
            dataIn[2] = rate[i];
            dataIn[3] = volatility[i];
            dataIn[4] = otime[i];
            dataIn[5] = otype[i];

            //#pragma parrot(input, "blackscholes", [6]dataIn)

            price_orig = BlkSchlsEqEuroNoDiv(sptprice[i], strike[i],
                                             rate[i], volatility[i], otime[i],
                                             otype[i], 0, &N1, &N2);
            dataOut[0] = price_orig;

            //#pragma parrot(output, "blackscholes", [1]<0.1; 0.9>dataOut)

            price_orig = dataOut[0];
            prices[i] = price_orig;
        }
    }
    return 0;
}

int main(int argc, char **argv)
{
    FILE *file;
    int i;
    int loopnum;
    fptype *buffer;
    int *buffer2;
    int rv;

    char *inputFile = argv[1];
    char *outputFile = argv[2];

    rv = 1;
    numOptions = 40;

    OptionData inputData[40] = {
        {40},
        {42.00, 40.00, 0.1000, 0.00, 0.20, 0.50, 'C', 0.00, 4.759423036851750055},
        {42.00, 40.00, 0.1000, 0.00, 0.20, 0.50, 'P', 0.00, 0.808600016880314021},
        {100.00, 100.00, 0.0500, 0.00, 0.15, 1.00, 'P', 0.00, 3.714602051381290071},
        {60.00, 65.00, 0.0800, 0.00, 0.30, 0.25, 'C', 0.00, 2.133371966735750025},
        {100.00, 90.00, 0.1000, 0.00, 0.10, 0.10, 'C', 0.00, 10.895610714793999563},
        {100.00, 90.00, 0.1000, 0.00, 0.10, 0.50, 'C', 0.00, 14.421570828843300660},
        {100.00, 90.00, 0.1000, 0.00, 0.10, 1.00, 'C', 0.00, 18.630859120667498274},
        {100.00, 100.00, 0.1000, 0.00, 0.10, 0.10, 'C', 0.00, 1.814984118378420108},
        {100.00, 100.00, 0.1000, 0.00, 0.10, 0.50, 'C', 0.00, 5.850273604284979889},
        {100.00, 100.00, 0.1000, 0.00, 0.10, 1.00, 'C', 0.00, 10.308147243666800463},
        {100.00, 110.00, 0.1000, 0.00, 0.10, 0.10, 'C', 0.00, 0.003523074865584340},
        {100.00, 110.00, 0.1000, 0.00, 0.10, 0.50, 'C', 0.00, 1.140722843827409960},
        {100.00, 110.00, 0.1000, 0.00, 0.10, 1.00, 'C', 0.00, 4.216747020308850402},
        {100.00, 90.00, 0.1000, 0.00, 0.25, 0.10, 'C', 0.00, 11.135244618346700207},
        {100.00, 90.00, 0.1000, 0.00, 0.25, 0.50, 'C', 0.00, 16.092638844092299166},
        {100.00, 90.00, 0.1000, 0.00, 0.25, 1.00, 'C', 0.00, 21.163454658480098658},
        {100.00, 100.00, 0.1000, 0.00, 0.25, 0.10, 'C', 0.00, 3.659962660310000171},
        {100.00, 100.00, 0.1000, 0.00, 0.25, 0.50, 'C', 0.00, 9.582231441086729973},
        {100.00, 100.00, 0.1000, 0.00, 0.25, 1.00, 'C', 0.00, 14.975798441718900733},
        {100.00, 110.00, 0.1000, 0.00, 0.25, 0.10, 'C', 0.00, 0.589613262035412977},
        {100.00, 110.00, 0.1000, 0.00, 0.25, 0.50, 'C', 0.00, 5.123575416865319809},
        {100.00, 110.00, 0.1000, 0.00, 0.25, 1.00, 'C', 0.00, 10.160055944516599880},
        {100.00, 90.00, 0.1000, 0.00, 0.50, 0.10, 'C', 0.00, 12.919509619564699676},
        {100.00, 90.00, 0.1000, 0.00, 0.50, 0.50, 'C', 0.00, 21.438280655724401669},
        {100.00, 90.00, 0.1000, 0.00, 0.50, 1.00, 'C', 0.00, 28.643647264488201643},
        {100.00, 100.00, 0.1000, 0.00, 0.50, 0.10, 'C', 0.00, 6.779936664291260406},
        {100.00, 100.00, 0.1000, 0.00, 0.50, 0.50, 'C', 0.00, 16.263193147074300526},
        {100.00, 100.00, 0.1000, 0.00, 0.50, 1.00, 'C', 0.00, 23.926745214162700393},
        {100.00, 110.00, 0.1000, 0.00, 0.50, 0.10, 'C', 0.00, 3.061909566605539812},
        {100.00, 110.00, 0.1000, 0.00, 0.50, 0.50, 'C', 0.00, 12.155688815568700178},
        {100.00, 110.00, 0.1000, 0.00, 0.50, 1.00, 'C', 0.00, 19.929858372066298955},
        {100.00, 90.00, 0.1000, 0.00, 0.10, 0.10, 'P', 0.00, 0.000095752219082624},
        {100.00, 90.00, 0.1000, 0.00, 0.10, 0.50, 'P', 0.00, 0.032219033907538303},
        {100.00, 90.00, 0.1000, 0.00, 0.10, 1.00, 'P', 0.00, 0.066226743903883001},
        {100.00, 100.00, 0.1000, 0.00, 0.10, 0.10, 'P', 0.00, 0.819967493295226002},
        {100.00, 100.00, 0.1000, 0.00, 0.10, 0.50, 'P', 0.00, 0.973216054356396021},
        {100.00, 100.00, 0.1000, 0.00, 0.10, 1.00, 'P', 0.00, 0.791889047262799961},
        {100.00, 110.00, 0.1000, 0.00, 0.10, 0.10, 'P', 0.00, 8.909004787274069415},
        {100.00, 110.00, 0.1000, 0.00, 0.10, 0.50, 'P', 0.00, 5.775959538905960144}
        };

    // alloc spaces for the option data
    data = (OptionData *)malloc(numOptions * sizeof(OptionData));
    prices = (fptype *)malloc(numOptions * sizeof(fptype));
    for (loopnum = 0; loopnum < numOptions; ++loopnum)
    {
        data[loopnum].s = inputData[loopnum].s;
        data[loopnum].strike = inputData[loopnum].strike;
        data[loopnum].r = inputData[loopnum].r;
        data[loopnum].divq = inputData[loopnum].divq;
        data[loopnum].v = inputData[loopnum].v;
        data[loopnum].t = inputData[loopnum].t;
        data[loopnum].OptionType = inputData[loopnum].OptionType;
        data[loopnum].divs = inputData[loopnum].divs;
        data[loopnum].DGrefval = inputData[loopnum].DGrefval;

        rv = 9;

    }
#define PAD 256
#define LINESIZE 64

    buffer = (fptype *)malloc(5 * numOptions * sizeof(fptype) + PAD);
    sptprice = (fptype *)(((unsigned long long)buffer + PAD) & ~(LINESIZE - 1));
    strike = sptprice + numOptions;
    rate = strike + numOptions;
    volatility = rate + numOptions;
    otime = volatility + numOptions;

    buffer2 = (int *)malloc(numOptions * sizeof(fptype) + PAD);
    otype = (int *)(((unsigned long long)buffer2 + PAD) & ~(LINESIZE - 1));

    for (i = 0; i < numOptions; i++)
    {
        otype[i] = (data[i].OptionType == 'P') ? 1 : 0;
        sptprice[i] = data[i].s / DIVIDE;
        strike[i] = data[i].strike / DIVIDE;
        rate[i] = data[i].r;
        volatility[i] = data[i].v;
        otime[i] = data[i].t;
    }

    //serial version
    int tid = 0;
    bs_thread(&tid);

    //Write prices to output file
    //prices is the final output
    for (i = 0; i < numOptions; i++)
    {
        printf("%.18f\n", prices[i]);
    }

    free(data);
    free(prices);

    return 0;
}