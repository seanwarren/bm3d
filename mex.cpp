#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <math.h>
#include <string.h>

#include "bm3d.h"
#include "utilities.h"

#include <mex.h>

#define YUV       0
#define YCBCR     1
#define OPP       2
#define RGB       3
#define DCT       4
#define BIOR      5
#define HADAMARD  6
#define NONE      7

using namespace std;

// c: pointer to original argc
// v: pointer to original argv
// o: option name after hyphen
// d: default value (if NULL, the option takes no argument)
const char *pick_option(int nrhs, const mxArray *prhs[], const char *o, const char *d) {
  int id = d ? 1 : 0;

  const int buflen = 1024;
  char v[buflen];
  char vi[buflen];

  for (int i = 0; i < nrhs - id; i++)
  {
     if (mxIsChar(prhs[i]) && mxIsChar(prhs[i + 1]))
     {
        mxGetString(prhs[i], v, buflen);
        mxGetString(prhs[i+id], vi, buflen);
        if (0 == strcmp(v, o))
           return (d == NULL) ? v : vi;
     }
  }
  return d;
}

/**
 * @file   main.cpp
 * @brief  Main executable file. Do not use lib_fftw to
 *         process DCT.
 *
 * @author MARC LEBRUN  <marc.lebrun@cmla.ens-cachan.fr>
 * @author SEAN WARREN  <s.warren@garvan.org.au>
 */
void mexFunction(int nlhs, mxArray *plhs[],
   int nrhs, const mxArray *prhs[])
{
   //! Variables initialization
   const char *_tau_2D_hard = pick_option(nrhs, prhs, "tau_2d_hard", "bior");
   const char *_tau_2D_wien = pick_option(nrhs, prhs, "tau_2d_wien", "dct");
   const char *_color_space = pick_option(nrhs, prhs, "color_space", "opp");
   const char *_patch_size = pick_option(nrhs, prhs, "patch_size", "0"); // >0: overrides default
   const char *_nb_threads = pick_option(nrhs, prhs, "nb_threads", "0");
   const bool useSD_1 = pick_option(nrhs, prhs, "useSD_hard", NULL) != NULL;
   const bool useSD_2 = pick_option(nrhs, prhs, "useSD_wien", NULL) != NULL;
   const bool verbose = pick_option(nrhs, prhs, "verbose", NULL) != NULL;

   //! Check parameters
   const unsigned tau_2D_hard = (strcmp(_tau_2D_hard, "dct") == 0 ? DCT :
      (strcmp(_tau_2D_hard, "bior") == 0 ? BIOR : NONE));
   if (tau_2D_hard == NONE)
   {
      mexWarnMsgIdAndTxt("bm3d", "tau_2d_hard is not known."); return;
   }

   const unsigned tau_2D_wien = (strcmp(_tau_2D_wien, "dct") == 0 ? DCT :
      (strcmp(_tau_2D_wien, "bior") == 0 ? BIOR : NONE));
   if (tau_2D_wien == NONE)
   {
      mexWarnMsgIdAndTxt("bm3d", "tau_2d_wien is not known"); return;
   }

   const unsigned color_space = (strcmp(_color_space, "rgb") == 0 ? RGB :
      (strcmp(_color_space, "yuv") == 0 ? YUV :
      (strcmp(_color_space, "ycbcr") == 0 ? YCBCR :
         (strcmp(_color_space, "opp") == 0 ? OPP : NONE))));
   if (color_space == NONE)
   {
      mexWarnMsgIdAndTxt("bm3d", "color_space is not known."); return;
   }

   const int patch_size = atoi(_patch_size);
   if (patch_size < 0)
   {
      mexWarnMsgIdAndTxt("bm3d", "The patch_size parameter must not be negative."); return;
   }
   
   const int nb_threads = atoi(_nb_threads);
   if (nb_threads < 0)
   {
      mexWarnMsgIdAndTxt("bm3d", "The nb_threads parameter must not be negative."); return;
   }


   //! Declarations
   vector<float> img_noisy, img_basic, img_denoised;
   unsigned width, height, chnls;

   mxAssert(nrhs >= 2, "usage: input sigma [basic]\n\
      [-tau_2d_hard{ dct,bior } (default: bior)]\n\
      [-useSD_hard]\n\
      [-tau_2d_wien{ dct,bior } (default: dct)]\n\
      [-useSD_wien]\n\
      [-color_space{ rgb,yuv,opp,ycbcr } (default: opp)]\n\
      [-patch_size{ 0,8,... } (default: 0, auto size, 8 or 12 depending on sigma)]\n\
      [-nb_threads(default: 0, auto number)]\n\
      [-verbose]");

   mxAssert(mxIsSingle(prhs[0]), "Input image must be of type single");
   int ndim = mxGetNumberOfDimensions(prhs[0]);
   const mwSize* dims = mxGetDimensions(prhs[0]);
   width = dims[0];
   height = dims[1];
   
   if (ndim == 3)
      chnls = dims[2];
   else if (ndim > 3)
   {
      mexWarnMsgIdAndTxt("bm3d", "Expected data to be 2 or 3 dimensional"); return;
   }

   img_noisy.resize(width * height * chnls);
   std::copy_n((float*) mxGetData(prhs[0]), width * height * chnls, img_noisy.begin());

   mxAssert(mxIsNumeric(prhs[1]), "Sigma must be num");
   float fSigma = mxGetScalar(prhs[1]); 

   //! Denoising
   if (run_bm3d(fSigma, img_noisy, img_basic, img_denoised, width, height, chnls,
      useSD_1, useSD_2, tau_2D_hard, tau_2D_wien, color_space, patch_size,
      nb_threads, verbose)
      != EXIT_SUCCESS)
   {
      mexWarnMsgIdAndTxt("bm3d", "BM3D failed"); return;
   }

   if (nlhs >= 0)
   {
      mwSize dims[3] = { width, height, chnls };
      plhs[0] = mxCreateNumericArray(3, dims, mxSINGLE_CLASS, mxREAL);
      std::copy(img_basic.begin(), img_basic.end(), (float*) mxGetData(plhs[0]));
   }
}
