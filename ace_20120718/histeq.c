/**
 * @file histeq.c
 * @brief Histogram equalization
 * @author Pascal Getreuer <getreuer@gmail.com>
 *
 * 
 * Copyright (c) 2011-2012, Pascal Getreuer
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

#include <ctype.h>
#include <stdio.h>
#include <string.h>
#include "imageio.h"

#define DEFAULT_NUMBINS           256

/** @brief struct of program parameters */
typedef struct
{
    /** @brief Input file name */
    char *InputFile;
    /** @brief Output file name */
    char *OutputFile;
    /** @brief Quality for saving JPEG images (0 to 100) */
    int JpegQuality;
    
    /** @brief Number of histogram bins */
    long NumBins;   
} programparams;


static int ParseParams(programparams *Param, int argc, char *argv[]);

static void PrintHelpMessage()
{
    puts("Histogram equalization, P. Getreuer 2011\n");
    puts("Usage: histeq [options] <input file> <output file>\n"
        "Only " READIMAGE_FORMATS_SUPPORTED " images are supported.\n");
    puts("Options:\n");
    puts("   -b <number>     number of histogram bins (default 256)");
#ifdef USE_LIBJPEG
    puts("   -q <number>     quality for saving JPEG images (0 to 100)\n");
#endif
}


int EqualizeImage(float *Image, int Width, int Height, int NumBins)
{       
    const long NumPixels = ((long)Width) * ((long)Height);
    const long NumEl = 3 * NumPixels;
    const int NumBinsMinusOne = NumBins - 1;
    float *Histogram1D[3] = {NULL, NULL, NULL};
    double CumHist;
    long i;
    int n, Channel, Success = 0;
    
    if(!(Histogram1D[0] = (float *)Malloc(sizeof(float)*3*NumBins)))
        goto Catch;
    
    Histogram1D[1] = Histogram1D[0] + NumBins;
    Histogram1D[2] = Histogram1D[0] + 2*NumBins;
    
    for(Channel = 0; Channel < 3; Channel++)
        for(i = 0; i < NumBins; i++)
            Histogram1D[Channel][i] = 0;
    
    /* Accumate channel histograms */
    for(i = 0; i < NumEl; i += 3)
    {
        Histogram1D[0][(int)(Image[i + 0]*NumBinsMinusOne + 0.5f)]++;
        Histogram1D[1][(int)(Image[i + 1]*NumBinsMinusOne + 0.5f)]++;
        Histogram1D[2][(int)(Image[i + 2]*NumBinsMinusOne + 0.5f)]++;
    }
    
    for(Channel = 0; Channel < 3; Channel++)
        for(i = 0; i < NumBins; i++)
            Histogram1D[Channel][i] /= NumPixels;
    
    /* Convert histograms to equalization maps */
    for(Channel = 0; Channel < 3; Channel++)
    {
        for(n = 0, CumHist = 0; n < NumBins; n++)
        {
            CumHist += Histogram1D[Channel][n];
            Histogram1D[Channel][n] = CumHist - Histogram1D[Channel][n]/2;
        }
        
        Histogram1D[Channel][0] = 0;
        Histogram1D[Channel][NumBins - 1] = 1;
    }

    /* Equalize the image */
    for(i = 0; i < NumEl; i += 3)
    {
        Image[i + 0] = Histogram1D[0][
            (int)(Image[i + 0]*NumBinsMinusOne + 0.5f)];
        Image[i + 1] = Histogram1D[1][
            (int)(Image[i + 1]*NumBinsMinusOne + 0.5f)];
        Image[i + 2] = Histogram1D[2][
            (int)(Image[i + 2]*NumBinsMinusOne + 0.5f)];
    }
    
    Success = 1;
Catch:    
    Free(Histogram1D[0]);
    return Success;
}


int main(int argc, char **argv)
{
    programparams Param;
    float *Image = NULL;
    unsigned long TimeStart;
    int Width, Height, Status = 0;
    
    if(!ParseParams(&Param, argc, argv))
        return 0;
    
    if(!(Image = (float *)ReadImage(&Width, &Height, Param.InputFile,
        IMAGEIO_RGB | IMAGEIO_FLOAT)))
        goto Catch;
    
    TimeStart = Clock();
    
    if(!EqualizeImage(Image, Width, Height, Param.NumBins))
        goto Catch;    
    
    printf("CPU Time: %.3f s\n", 0.001f*(Clock() - TimeStart));
    
    if(!WriteImage(Image, Width, Height, Param.OutputFile,
        IMAGEIO_RGB | IMAGEIO_FLOAT, Param.JpegQuality))
        goto Catch;
    
    Status = 0;
Catch:
    if(Image)
        free(Image);
    return Status;
}


static int ParseParams(programparams *Param, int argc, char *argv[])
{
    static char *DefaultOutputFile = (char *)"out.png";
    char *OptionString;
    char OptionChar;
    int i;


    if(argc < 2)
    {
        PrintHelpMessage();
        return 0;
    }

    /* Set parameter defaults */
    Param->InputFile = 0;
    Param->OutputFile = DefaultOutputFile;
    Param->JpegQuality = 85;
    Param->NumBins = DEFAULT_NUMBINS;
    
    for(i = 1; i < argc;)
    {
        if(argv[i] && argv[i][0] == '-')
        {
            if((OptionChar = argv[i][1]) == 0)
            {
                ErrorMessage("Invalid parameter format.\n");
                return 0;
            }

            if(argv[i][2])
                OptionString = &argv[i][2];
            else if(++i < argc)
                OptionString = argv[i];
            else
            {
                ErrorMessage("Invalid parameter format.\n");
                return 0;
            }

            switch(OptionChar)
            {                 
            case 'b':
                Param->NumBins = atoi(OptionString);
                
                if(Param->NumBins <= 1 || Param->NumBins > 10000)
                {
                    ErrorMessage("Number of bins must be between 2 and 10000.\n");
                    return 0;
                }
                break;            
#ifdef LIBJPEG_SUPPORT
            case 'q':
                Param->JpegQuality = atoi(OptionString);

                if(Param->JpegQuality <= 0 || Param->JpegQuality > 100)
                {
                    fprintf(stderr, "JPEG quality must be between 0 and 100.\n");
                    return 0;
                }
                break;
#endif
            case '-':
                PrintHelpMessage();
                return 0;
            default:
                if(isprint(OptionChar))
                    ErrorMessage("Unknown option \"-%c\".\n", OptionChar);
                else
                    ErrorMessage("Unknown option.\n");

                return 0;
            }

            i++;
        }
        else
        {
            if(!Param->InputFile)
                Param->InputFile = argv[i];
            else
                Param->OutputFile = argv[i];

            i++;
        }
    }

    if(!Param->InputFile)
    {
        PrintHelpMessage();
        return 0;
    }
    return 1;
}
