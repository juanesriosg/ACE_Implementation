/**
 * @file acecli.c
 * @brief ACE automatic color enhancement command line program
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

/**
 * @mainpage
 * @verbinclude readme.txt
 */

#include <math.h>
#include <string.h>
#include <ctype.h>
#include "ace.h"
#include "imageio.h"

#define VERBOSE 0


/** @brief struct of program parameters */
typedef struct
{    
    /** @brief Input file name */
    char *InputFile;
    /** @brief Output file name */
    char *OutputFile;
    /** @brief Quality for saving JPEG images (0 to 100) */
    int JpegQuality;   
    /** @brief Slope parameter, larger implies stronger in enhancement */
    float Alpha;
} programparams;


int ParseParams(programparams *Param, int argc, char *argv[]);

/** @brief Print program usage help message */
void PrintHelpMessage()
{
    puts("ACE automatic color enhancement , P. Getreuer 2012\n");
    puts("Usage: ace [options] <input file> <output file>\n\n"
        "Only " READIMAGE_FORMATS_SUPPORTED " images are supported.\n");
    puts("Options:");
    puts("  -a <number>  alpha, stronger implies stronger enhancement\n");
#ifdef USE_LIBJPEG
    puts("  -q <number>  quality for saving JPEG images (0 to 100)\n");
#endif
    puts("Example: ");
    puts("   ace input.bmp output.bmp");
}


int main(int argc, char *argv[])
{
    programparams Param;
    float *f = NULL, *u = NULL;
    unsigned long TimeStart;
    int Width, Height;
    int Status = 1;
    
    
    if(!ParseParams(&Param, argc, argv))
        return 0;
            
    /* Read the input image */
    if(!(f = (float *)ReadImage(&Width, &Height, Param.InputFile,
        IMAGEIO_FLOAT | IMAGEIO_PLANAR | IMAGEIO_RGB)))
        goto Catch;

    /* Allocate the output image */
    if(!(u = (float *)Malloc(sizeof(float)*3*
        ((long int)Width)*((long int)Height))))
        goto Catch;
    
    printf("Enhancing %dx%d image, alpha = %.4f\n", 
        Width, Height, Param.Alpha);
    TimeStart = Clock();
       
    /* ACE enhancement */
    if(!AceEnhanceImage(u, f, Width, Height, Param.Alpha))
    {
        ErrorMessage("Error in computation.\n");
        goto Catch;  
    }
    
    printf("CPU Time: %.3f s\n", 0.001f*(Clock() - TimeStart));
    
    /* Write the output image */
    if(!WriteImage(u, Width, Height, Param.OutputFile, 
        IMAGEIO_FLOAT | IMAGEIO_PLANAR | IMAGEIO_RGB, Param.JpegQuality))
        goto Catch;
#if VERBOSE > 0
    else
        printf("Output written to \"%s\".\n", Param.OutputFile);
#endif
    
    Status = 0;	/* Finished successfully, set exit status to zero. */
    
Catch:
    Free(u);
    Free(f);
    return Status;
}


int ParseParams(programparams *Param, int argc, char *argv[])
{
    static char *DefaultOutputFile = (char *)"out.bmp";
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
    Param->Alpha = 5;
    
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
            case 'a':   /* Read slope parameter alpha */
                Param->Alpha = atof(OptionString);
                
                if(Param->Alpha < 1 || Param->Alpha > 8)
                {
                    ErrorMessage("Alpha must be between 1 and 8.\n");
                    return 0;
                }
                break;
#ifdef USE_LIBJPEG
            case 'q':
                Param->JpegQuality = atoi(OptionString);

                if(Param->JpegQuality <= 0 || Param->JpegQuality > 100)
                {
                    ErrorMessage("JPEG quality must be between 0 and 100.\n");
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
