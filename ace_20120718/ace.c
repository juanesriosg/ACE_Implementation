/**
 * @file ace.c
 * @brief ACE automatic color enhancement
 * @author Pascal Getreuer <getreuer@gmail.com>
 * 
 * Copyright (c) 2012, Pascal Getreuer
 * All rights reserved.
 * 
 * This program is free software: you can redistribute it and/or modify it
 * under, at your option, the terms of the GNU General Public License as 
 * published by the Free Software Foundation, either version 3 of the 
 * License, or (at your option) any later version, or the terms of the 
 * simplified BSD license.
 *
 * You should have received a copy of these licenses along with this program.
 * If not, see <http://www.gnu.org/licenses/> and
 * <http://www.opensource.org/licenses/bsd-license.html>.
 */

#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <fftw3.h>
#ifdef _OPENMP
#include <omp.h>
#endif
#include "ace.h"

/** @brief Degree of the polynomial approximation, must be odd */
#define DEGREE          9


/** @brief Compute small factorials */
static double Factorial(int n)
{
    static const double Table[15] = {1, 1, 2, 6, 24, 120, 720, 5040, 40320,
        362880, 3628800, 39916800, 479001600, 6227020800, 87178291200};
    
    if(n < 0)
        return 0;
    else if(n < 15)
        return Table[n];
    else
        return n * Factorial(n - 1);
}


/** @brief Compute binomial coefficients */
static double BinomCoeff(int m, int n)
{
    return Factorial(m) / ( Factorial(n) * Factorial(m - n));
}


/** @brief FFT-based convolution */
static void Convolve(float *BlurredTrans, const float *OmegaTrans, 
    long NumPixels, fftwf_plan ForwardPlan, fftwf_plan InversePlan)
{
    long i;
    
    fftwf_execute(ForwardPlan);
    
    for(i = 0; i < NumPixels; i++)
        BlurredTrans[i] *= OmegaTrans[i];
    
    fftwf_execute(InversePlan);    
}


/** @brief Evaluate integer power, hardcoded for degrees 1 to 9 */
static void IntPow(float *Dest, const float *Src, size_t NumSamples, int m)
{
    size_t i;
    
    switch(m)
    {
    case 1:
        memcpy(Dest, Src, sizeof(float)*NumSamples);
        break;
    case 2:
        for(i = 0; i < NumSamples; i++)
        {
            float x = Src[i];
            Dest[i] = x * x;
        }
        break;
    case 3:
        for(i = 0; i < NumSamples; i++)
        {
            float x = Src[i];
            Dest[i] = x * x * x;
        }
        break;
    case 4:
        for(i = 0; i < NumSamples; i++)
        {
            float x = Src[i];
            float x2 = x * x;
            Dest[i] = x2 * x2;
        }
        break;
    case 5:
        for(i = 0; i < NumSamples; i++)
        {
            float x = Src[i];
            float x2 = x * x;
            Dest[i] = x2 * x2 * x;
        }
        break;
    case 6:
        for(i = 0; i < NumSamples; i++)
        {
            float x = Src[i];
            float x2 = x * x;
            Dest[i] = x2 * x2 * x2;
        }
        break;
    case 7:
        for(i = 0; i < NumSamples; i++)
        {
            float x = Src[i];
            float x2 = x * x;
            Dest[i] = x2 * x2 * x2 * x;
        }
        break;
    case 8:
        for(i = 0; i < NumSamples; i++)
        {
            float x = Src[i];
            float x2 = x * x;
            float x4 = x2 * x2;
            Dest[i] = x4 * x4;
        }
        break;
    case 9:
        for(i = 0; i < NumSamples; i++)
        {
            float x = Src[i];
            float x2 = x * x;
            float x4 = x2 * x2;
            Dest[i] = x4 * x4 * x;
        }
        break;
    default:
        for(i = 0; i < NumSamples; i++)
            Dest[i] = (float)pow(Src[i], m);
        break;
    }
}


/** @brief Compute the FFT of omega(x,y) = 1/sqrt(x^2 + y^2) */
static int ComputeOmegaTrans(float *OmegaTrans, 
    float *Omega, int Width, int Height)
{
    fftwf_plan Plan = NULL;
    long PadNumPixels = ((long)Width + 1) * ((long)Height + 1);
    double Sum = 0;
    long i, x, y;
    
    for(y = 0, i = 0; y <= Height; y++)
        for(x = 0; x <= Width; x++, i++)
        {
            Omega[i] = (x == 0 && y == 0) ? 0
                : 1.0f/sqrt(x*x + y*y);
            Sum += ((x == 0 || x == Width) ? 1 : 2) 
                * ((y == 0 || y == Height) ? 1 : 2) 
                * Omega[i];
        }
    
    Sum *= 4*PadNumPixels;
    
    for(i = 0; i < PadNumPixels; i++)
        Omega[i] /= Sum;
    
    if(!(Plan = fftwf_plan_r2r_2d(Height + 1, Width + 1, Omega, OmegaTrans,
        FFTW_REDFT00, FFTW_REDFT00, FFTW_ESTIMATE | FFTW_DESTROY_INPUT)))
        return 0;
    
    fftwf_execute(Plan);
    fftwf_destroy_plan(Plan);
    
    /* Cut last row and column from KernelTrans */
    for(y = 1, i = Width; y < Height; y++, i += Width)
        memmove(OmegaTrans + i, OmegaTrans + i + y, sizeof(float)*Width); 
    
    return 1;    
}


/**
 * @brief ACE automatic color enhancement
 * @param u the enhanced output image
 * @param f the input image in planar row-major order
 * @param Width, Height image dimensions
 * @param Alpha the slope parameter (>=1), larger implies stronger enhancement
 * 
 * This routine perform ACE enhancement using the fast O(N log N) algorithm of
 * Bertalmio et al.  The slope parameter must be an integer or half-integer 
 * between 1 and 8 (1, 1.5, 2, 2.5, ..., 7.5, 8).  The slope function is 
 * approximated by an \f$L^infinity\f$-optimal degree 9 polynomial.
 * 
 * If OpenMP is enabled, the main computation loop is parallelized.  Note: in
 * the parallel computation, values are summed in a nondeterministic order 
 * (i.e., in whichever order threads finish the computations).  Due to 
 * rounding effects, addition is not exactly associative and the output varies
 * slightly between runs (+/- 1 intensity level).
 */
int AceEnhanceImage(float *u, const float *f, 
    int Width, int Height, float Alpha)
{
    /* L^infinity-optimal degree 9 polynomial coefficients.  Since polynomials
       are odd, only the odd term coefficients are stored. */
    static const float SlopeCoeff[15][(DEGREE + 1)/2] = {
    /* Alpha = 1 */
    {1.0f, 0, 0, 0, 0},    
    /* Alpha = 1.5 */    
    {1.33743875f, 1.55213754f, -3.02825657f, -0.12350511f, 1.28325061f},
    /* Alpha = 2 */
    {1.85623249f, 3.82397125f, -19.70879455f, 26.15510902f, -11.15375327f},
    /* Alpha = 2.5 */
    {2.79126397f, -1.30687551f, -10.57298680f, 20.02623286f, -9.98284231f},
    /* Alpha = 3 */
    {3.51036396f, -6.31644952f, 0.92439798f, 9.32834829f, -6.50264005f},
    /* Alpha = 3.5 */
    {4.15462973f, -11.85851451f, 16.03418150f, -7.07985902f, -0.31040920f},
    /* Alpha = 4 */
    {4.76270090f, -18.23743983f, 36.10529118f, -31.35677926f, 9.66532431f},
    /* Alpha = 4.5 */
    {5.34087782f, -25.67018163f, 63.87617747f, -70.15437134f, 27.66951403f},
    /* Alpha = 5 */
    {5.64305564f, -28.94026159f, 74.52401661f, -83.54012582f, 33.39343065f},
    /* Alpha = 5.5 */
    {5.92841230f, -32.11619291f, 85.01764165f, -96.84966316f, 39.11863693f},
    /* Alpha = 6 */
    {6.19837979f, -35.18789052f, 95.28157108f, -109.95601312f, 44.78177264f},
    /* Alpha = 6.5 */
    {6.45529995f, -38.16327397f, 105.31193936f, -122.83169063f, 50.36462504f},
    /* Alpha = 7 */
    {6.69888108f, -41.02503190f, 115.02784036f, -135.35603880f, 55.81014424f},
    /* Alpha = 7.5 */
    {6.92966632f, -43.76867314f, 124.39645141f, -147.47363378f, 61.09053024f},
    /* Alpha = 8 */    
    {7.15179080f, -46.43557440f, 133.54648929f, -159.34156394f, 66.27157886f}
    };
#ifdef _OPENMP
    const int NumThreads = omp_get_max_threads();
#else
    const int NumThreads = 1;
#endif
    const long NumPixels = ((long)Width) * ((long)Height);
    const long PadNumPixels = ((long)Width + 1) * ((long)Height + 1);
    float *(*Arrays)[2] = NULL;
    float *OmegaTrans = NULL, *PolyCoeffs = NULL;
    fftwf_plan (*Plans)[2] = NULL;    
    long i;
    int Iter, n, Channel, Success = 0;
    
    /* Allocate memory */
    if(!(OmegaTrans = (float *)fftwf_malloc(sizeof(float)*PadNumPixels))
        || !ComputeOmegaTrans(OmegaTrans, u, Width, Height)
        || !(PolyCoeffs = (float *)malloc(sizeof(float)*10*10)))
        goto Catch;
    
    if(!(Arrays = (float *(*)[2])malloc(sizeof(float *)*2*NumThreads)))
        goto Catch;
    
    for(i = 0; i < NumThreads; i++)
        Arrays[i][0] = Arrays[i][1] = NULL;
    
    if(!(Plans = (fftwf_plan (*)[2])malloc(sizeof(fftwf_plan)*2*NumThreads)))
        goto Catch;
    
    for(i = 0; i < NumThreads; i++)
        Plans[i][0] = Plans[i][1] = NULL;
    
    /* For each thread, allocate a workspace array and create FFTW transform
       plans.  These will be used for computing convolutions in parallel. */
    for(i = 0; i < NumThreads; i++)
        if(!(Arrays[i][0] = (float *)
            fftwf_malloc(sizeof(float)*NumPixels))
            || !(Arrays[i][1] = (float *)
            fftwf_malloc(sizeof(float)*NumPixels))
            /* DCT-II with source Arrays[i][0] and dest Arrays[i][1] */
            || !(Plans[i][0] = fftwf_plan_r2r_2d(Height, Width, 
            Arrays[i][0], Arrays[i][1], FFTW_REDFT10, FFTW_REDFT10, 
            FFTW_ESTIMATE | FFTW_DESTROY_INPUT))
            /* DCT-III with source Arrays[i][1] and dest Arrays[i][0] */
            || !(Plans[i][1] = fftwf_plan_r2r_2d(Height, Width, 
            Arrays[i][1], Arrays[i][0], FFTW_REDFT01, FFTW_REDFT01, 
            FFTW_ESTIMATE | FFTW_DESTROY_INPUT)))
        goto Catch;
    
    /* Select polynomial from SlopeCoeff table */
    i = (int)(2*Alpha + 0.5f) - 2;
    
    if(i < 0)
        i = 0;
    else if(i > 14)
        i = 14;
    
    /* Precompute coefficients */
    for(n = 0; n <= DEGREE; n++)
    {
        int m;
        
        for(m = n + ((n % 2 == 0) ? 1 : 0); m <= DEGREE; m += 2)
            PolyCoeffs[10*n + m] = pow(-1, m - n + 1)
                * SlopeCoeff[i][(m - 1)/2] * BinomCoeff(m, n);
    }

    /* Special case for n = zero term */
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic) private(i)
#endif
    for(Channel = 0; Channel < 3; Channel++)
    {
        const float *Src = f + NumPixels * Channel;
        float *Dest = u + NumPixels * Channel;
        
        for(i = 0; i < NumPixels; i++)
        {
            float Temp = Src[i];
            float TempSqr = Temp*Temp;                
            float a = PolyCoeffs[DEGREE];
            int m = DEGREE;
            
            while(m >= 2)
                a = a*TempSqr + PolyCoeffs[m -= 2];
            
            Dest[i] = a * Temp;
        }
    }
    
    /* Most of the computation time is spent in this loop. */
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic) private(i)
#endif
    for(Iter = 0; Iter < 3*DEGREE; Iter++)
    {
#ifdef _OPENMP
        const int ThreadId = omp_get_thread_num();
#else
        const int ThreadId = 0;
#endif
        int Channel = Iter / DEGREE;    /* = current image channel  */
        int n = 1 + (Iter % DEGREE);    /* = current summation term */
        const float *Src = f + NumPixels * Channel;
        float *Dest = u + NumPixels * Channel;
        float *Blurred = Arrays[ThreadId][0];

        /* Compute Blurred = Src to the nth power */
        IntPow(Blurred, Src, NumPixels, n);
        /* Convolve Blurred with Omega */
        Convolve(Arrays[ThreadId][1], OmegaTrans, NumPixels,
            Plans[ThreadId][0], Plans[ThreadId][1]);
        
        for(i = 0; i < NumPixels; i++)
        {
            float Temp = Src[i];
            float TempSqr = Temp*Temp;
            float a = PolyCoeffs[(DEGREE + 1)*n + DEGREE];
            int m = DEGREE;
            
            while(m - n >= 2)
            {
                m -= 2;
                a = a*TempSqr + PolyCoeffs[(DEGREE + 1)*n + m];
            }
            
            if(n % 2 == 0)
                a *= Temp;
            
#ifdef _OPENMP            
            Blurred[i] *= a;
#else
            Dest[i] += a * Blurred[i];
#endif
        }
        
#ifdef _OPENMP
#pragma omp critical
        for(i = 0; i < NumPixels; i++)
            Dest[i] += Blurred[i];
#endif        
    }

#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic) private(i)
#endif
    for(Channel = 0; Channel < 3; Channel++)
    {
        float *Dest = u + NumPixels * Channel;
        float MaxValue = -1e30f;
        
        for(i = 0; i < NumPixels; i++)
            if(Dest[i] > MaxValue)
                MaxValue = Dest[i];
        
        MaxValue *= 2;
        
        for(i = 0; i < NumPixels; i++)
            Dest[i] = 0.5f + Dest[i] / MaxValue;
    }
    
    Success = 1;
Catch:
    /* Free memory */
    for(i = 0; i < NumThreads; i++)
        for(n = 1; n >= 0; n--)
        {
            if(Plans[i][n])
                fftwf_destroy_plan(Plans[i][n]);
            if(Arrays[i][n])
                fftwf_free(Arrays[i][n]);
        }
        
    if(Plans)
        free(Plans);
    if(Arrays)
        free(Arrays);
    if(PolyCoeffs)
        free(PolyCoeffs);
    if(OmegaTrans)
        fftwf_free(OmegaTrans);
    fftwf_cleanup();
    return Success;
}
