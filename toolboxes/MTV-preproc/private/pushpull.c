#include <math.h>
#include "mex.h"
#include "shoot_boundary.h"
#include "fast.h" // fast approximations of exp/log/erf/...

//-------------------------------------------------------------------------
// Constants
//-------------------------------------------------------------------------

static const double one_div_root_two_pi = 3.989422804014326779399460599343818684e-01;
static const double one_div_root_two    = 7.071067811865475244008443621048490392e-01;
static const double ln_two              = 6.931471805599453094172321214581765680e-01;
static const double one_div_ln_two      = 1.442695040888963387004650940070860088;
static const double one_div_root_ln_two = 1.201122408786449824447117862291634083;
static const double two_div_root_pi     = 1.1283791670955125738961589031215452;

//-------------------------------------------------------------------------
// Window functions
//-------------------------------------------------------------------------

static float wingaussfast(float x2)
{
    return(fasterexp(-0.5*x2));
}
// DIRAC (actually, DIRAC * GAUSSIAN)
static float windirac(float x, float scl /* unused */, float sigbasis)
{
    x /= sigbasis;
    x *= x;
    return(wingaussfast(x));
}
static float windiracnorm(float scl /* unused */, float sigbasis)
{
    return(one_div_root_two_pi/sigbasis);
}
static float windiraclim(float scl /* unused */, float sigbasis)
{
    return(3*sigbasis);
}
static float windiracderiv(float x, float scl /* unused */, float sigbasis)
{
    sigbasis *= sigbasis;
    return(-(x/sigbasis)*wingaussfast(x*x/sigbasis));
}
// GAUSSIAN (actually, GAUSSIAN * GAUSSIAN)
static float wingauss(float x, float scl, float sigbasis)
{
    static const float default_sig2 = 0.1250*one_div_ln_two;
    
    x = (x*x)/(scl*scl*default_sig2+sigbasis*sigbasis);
    return(wingaussfast(x));
}
static float wingaussnorm(float scl, float sigbasis)
{
    static const float default_sig2 = 0.1250*one_div_ln_two;
    
    return(one_div_root_two_pi/sqrt(scl*scl*default_sig2+sigbasis*sigbasis));
}
static float wingausslim(float scl, float sigbasis)
{
    static const float default_sig2 = 0.1250*one_div_ln_two;
    
    return(3*sqrt((scl*scl*default_sig2+sigbasis*sigbasis)));
}
static float wingaussderiv(float x, float scl, float sigbasis)
{
    static const float default_sig2 = 0.1250*one_div_ln_two;
    
    sigbasis = (scl*scl*default_sig2+sigbasis*sigbasis);
    return(-(x/sigbasis)*wingaussfast((x*x)/sigbasis));
}
// RECT (actually, RECT * GAUSSIAN)
static float winrect(float x, float scl, float sigbasis)
{
    float w1 = one_div_root_two/sigbasis;
    scl = 0.5*scl;
    return(fastererfc(w1*(x-scl))-fastererfc(w1*(x+scl)));   
}
static float winrectnorm(float scl, float sigbasis)
{
    return(0.5/scl);
}
static float winrectlim(float scl, float sigbasis)
{
    return(0.5*scl+3*sigbasis);
}
static float winrectderiv(float x, float scl, float sigbasis)
{
    float w1 = one_div_root_two/sigbasis;
    scl = 0.5*scl;
    float xm = w1*(x-scl);
    float xp = w1*(x+scl);
    return(two_div_root_pi*w1*(fasterexp(-xp*xp)-fasterexp(-xm*xm)));   
}

//-------------------------------------------------------------------------
// Main function
//-------------------------------------------------------------------------

/** Pull or push an image according to a deformation and a window
 *
 * @param dm0   [mwSize 3]          Reference dim    (pull:in / push:out)
 * @param  m1   [mwSize]            Nb pulled voxels (pull:out / push:in)
 * @param  n    [mwSize]            Nb features      (4th dimension)
 * @param Psi   [float m1*3]        Deformation
 * @param  F0   [float prod(dm0)*n] Reference volume (pull:in / push:out)
 * @param  S0   [float]             Count volume     (push:out)
 * @param  F1   [float m1*n]        Pulled volume    (pull:out / push:in)
 * @param code  [uint]              (0=push|1=pull|2=pushc|3=pullc)
 * @param Jac   [float mJ*3*3]      Jacobian of the deformation
 * @param  mJ   [mwSize]            Number of Jacobians (1 or m1)
 * @param func  [double 3]          Selection window
 *                                  (0=dirac|1=gauss|2=rect)
 * @param  D1x  [float m1*n]        Pulled derivatives / x
 * @param  D1y  [float m1*n]        Pulled derivatives / y
 * @param  D1z  [float m1*n]        Pulled derivatives / z
 *
 * The Jacobian must map pulled voxels to reference voxels and it should 
 * be possible to write it as J=R*S, where R is a rotation matrix (R*R'=I) 
 * and S is a diagonal scale matrix.
 */
#define TINY 5e-2f
static void jpushpull(mwSize dm0[], mwSize m1, mwSize n, 
                      float Psi[], float F0[], float S0[], float F1[], 
                      unsigned int code, float Jac[], mwSize mJ, 
                      double win[], float D1x[], float D1y[], float D1z[])
{
    /* SD = 1/(2*sqrt(2*ln(2))) corresponding to FWHM = 1 */
    static const float default_sig = 0.5*one_div_root_two*one_div_root_ln_two;
    
    /* Which derivatives should be computed */
    const mwSize doD1x = (mwSize)D1x;
    const mwSize doD1y = (mwSize)D1y;
    const mwSize doD1z = (mwSize)D1z;
    
    /* Pointers into input/output arrays */
    const float   NaN = mxGetNaN();
    const mwSize  m0  = dm0[0]*dm0[1]*dm0[2], /* Number of reference voxels */
                  oy  = dm0[0],               /* Offsets between lines */
                  oz  = dm0[0]*dm0[1];        /* Offsets between slices */
    float  *px  = Psi;                  /* Deformation in x */
    float  *py  = Psi+m1;               /*                y */
    float  *pz  = Psi+m1*2;             /*                z */
    float  *jxx = Jac;                  /* Jacobian in xx */
    float  *jyx = Jac+mJ;               /*             yx */
    float  *jzx = Jac+mJ*2;             /*             zx */
    float  *jxy = Jac+mJ*3;             /*             xy */
    float  *jyy = Jac+mJ*4;             /*             yy */
    float  *jzy = Jac+mJ*5;             /*             zy */
    float  *jxz = Jac+mJ*6;             /*             xz */
    float  *jyz = Jac+mJ*7;             /*             yz */
    float  *jzz = Jac+mJ*8;             /*             zz */
    
    /* Stuff needed inside the loop */
    float  Jxx, Jyy, Jzz,               /* Jacobian, then Rotation */
           Jxy, Jxz, Jyz,
           Jyx, Jzx, Jzy;
    float  Sx, Sy, Sz;                  /* Scale */
    float  limx, limy, limz;            /* Bounding box in ref space */
    float  klimx, klimy, klimz;         /* Bounding box in kernel space */
    float  norm, normx, normy, normz;   /* Kernel normalisation */
    mwSize i, k, iz, iy, ix;            /* Loop iteration variables */
    
    /* Get slice selection function pointers */
    typedef float (*wintype3)(float,float,float);
    typedef float (*wintype2)(float,float);
    wintype3 *fwin   = (wintype3*)malloc(3*sizeof(wintype3)); /* Window */
    wintype2 *fnorm  = (wintype2*)malloc(3*sizeof(wintype2)); /* Normalisation */
    wintype2 *flim   = (wintype2*)malloc(3*sizeof(wintype2)); /* Support */
    wintype3 *fderiv = (wintype3*)malloc(3*sizeof(wintype3)); /* Derivative */
    for(i=0; i<3; ++i)
    {
        switch((unsigned int)(win[i]))
        {
            case 0:
                fwin[i]   = &windirac;
                fnorm[i]  = &windiracnorm;
                flim[i]   = &windiraclim;
                fderiv[i] = &windiracderiv;
                break;
            case 1:
                fwin[i]   = &wingauss;
                fnorm[i]  = &wingaussnorm;
                flim[i]   = &wingausslim;
                fderiv[i] = &wingaussderiv;
                break;
            case 2:
                fwin[i]   = &winrect;
                fnorm[i]  = &winrectnorm;
                flim[i]   = &winrectlim;
                fderiv[i] = &winrectderiv;
                break;
            default:
                free(fwin);
                free(fnorm);
                free(flim);
                free(fderiv);
                mexErrMsgTxt("Window type must be 0 (dirac), 1 (gauss) or 2 (rect).");
                break;
        }
    }
    
    const int voxJ = (mJ != 1);
    if (!voxJ)
    /* Invert Jacobian if same for all voxels */
    {
        /* Jacobian at point i */
        Jxx = *(jxx);
        Jyy = *(jyy);
        Jzz = *(jzz);
        Jxy = *(jxy);
        Jxz = *(jxz);
        Jyz = *(jyz);
        Jyx = *(jyx);
        Jzx = *(jzx);
        Jzy = *(jzy);
        
        /* Compute scale */
        Sx = sqrt(Jxx*Jxx + Jyx*Jyx + Jzx*Jzx);
        Sy = sqrt(Jxy*Jxy + Jyy*Jyy + Jzy*Jzy);
        Sz = sqrt(Jxz*Jxz + Jyz*Jyz + Jzz*Jzz);
        
        /* Compute rotation */
        Jxx /= Sx;
        Jyx /= Sx;
        Jzx /= Sx;
        Jxy /= Sy;
        Jyy /= Sy;
        Jzy /= Sy;
        Jxz /= Sz;
        Jyz /= Sz;
        Jzz /= Sz;
        
        /* Bounding box of contributing points in ref space */
        klimx = flim[0](Sx, default_sig),
        klimy = flim[1](Sy, default_sig),
        klimz = flim[2](Sz, default_sig);
        limx = fmax(fmax(fabs(Jxx*klimx),fabs(Jxy*klimy)),fabs(Jxz*klimz));
        limy = fmax(fmax(fabs(Jyx*klimx),fabs(Jyy*klimy)),fabs(Jyz*klimz));
        limz = fmax(fmax(fabs(Jzx*klimx),fabs(Jzy*klimy)),fabs(Jzz*klimz));
//         mexPrintf("Limits: %f %f %f\n", limx, limy, limz);
        
        normx = fnorm[0](Sx, default_sig);
        normy = fnorm[1](Sy, default_sig);
        normz = fnorm[2](Sz, default_sig);
        norm  = normx*normy*normz;
    }
    
    float *pf1  = F1;   /* Pointer to pulled volume */
    float *pd1x = D1x;  /* Pointers to pulled derivative */
    float *pd1y = D1y;
    float *pd1z = D1z;
    for (i=0; i<m1; ++i, ++pf1, ++pd1x, ++pd1y, ++pd1z)
    /* loop over pulled voxels */
    {
        float  x, y, z;     /* Coordinates in reference volume */
        x = *(px++)-1.0f;   /* Subtract 1 because of MATLAB indexing */
        y = *(py++)-1.0f;
        z = *(pz++)-1.0f;
        
        if (((code & 2)==2 && mxIsFinite(x) && mxIsFinite(y) && mxIsFinite(z))
            || ((x>=-TINY) && (x<=(float)(dm0[0])-1.0f+TINY)
            &&  (y>=-TINY) && (y<=(float)(dm0[1])-1.0f+TINY)
            &&  (z>=-TINY) && (z<=(float)(dm0[2])-1.0f+TINY)))
        /* If (pullc/pushc) OR (pull/push AND inside bounds) */
        {
            if (voxJ)
            /* Invert Jacobian if voxel-specific */
            {
                free(fwin);
                free(fnorm);
                free(flim);
                free(fderiv);
                mexErrMsgTxt("Non-stationary Jacobian not implemented yet");
            }
            
            /* Voxels inside the bounding box */
            mwSignedIndex  ixmax = floor(x+limx),
                           iymax = floor(y+limy),
                           izmax = floor(z+limz),
                           ixmin = ceil(x-limx),
                           iymin = ceil(y-limy),
                           izmin = ceil(z-limz);
            mwSize         lx    = ixmax - ixmin + 1,
                           ly    = iymax - iymin + 1,
                           lz    = izmax - izmin + 1;
            
            /* Create lookups of voxel locations - for coping with edges */
            mwSize  oox[128], ooy[128], ooz[128];
            for(k=0; k<lx; ++k) oox[k] = bound(ixmin+k, dm0[0]);
            for(k=0; k<ly; ++k) ooy[k] = bound(iymin+k, dm0[1])*oy;
            for(k=0; k<lz; ++k) ooz[k] = bound(izmin+k, dm0[2])*oz;
            float  ddx[128], ddy[128], ddz[128];
            for(k=0; k<lx; ++k) ddx[k] = x-(float)(k+ixmin);
            for(k=0; k<ly; ++k) ddy[k] = y-(float)(k+iymin);
            for(k=0; k<lz; ++k) ddz[k] = z-(float)(k+izmin);
            
            if ((code&1)==1)
            /* Pull */
            {
                            
                /* Loop over contributing voxels in ref space */
                for (iz=0; iz<lz;)
                {
                    float dz    = ddz[iz];
                    float *pf0z = F0 + ooz[iz++];
                    for (iy=0; iy<ly;)
                    {
                        float dy    = ddy[iy];
                        float *pf0y = pf0z + ooy[iy++];
                        for (ix=0; ix<lx;)
                        {
                            float dx    = ddx[ix];
                            float *pf0x = pf0y + oox[ix++];
                            
                            /* Rotate dx towards subject space: Td=R'*d */
                            float Tdx  = Jxx * dx + Jyx * dy + Jzx * dz;
                            float Tdy  = Jxy * dx + Jyy * dy + Jzy * dz;
                            float Tdz  = Jxz * dx + Jyz * dy + Jzz * dz;
                            
                            if (fabs(Tdx)<klimx && fabs(Tdy)<klimy && fabs(Tdz)<klimz)
                            /* If x is near enough to the basis centre */
                            {
                                /* Compute weight */
                                float w = norm * fwin[0](Tdx, Sx, default_sig) 
                                               * fwin[1](Tdy, Sy, default_sig) 
                                               * fwin[2](Tdz, Sz, default_sig);
                            
                                float wdx, wdy, wdz;
                                if (doD1x || doD1y || doD1z)
                                {
                                    /* Compute derivative weights in subject space */
                                    float Twdx = normx * fderiv[0](Tdx, Sx, default_sig);
                                    float Twdy = normy * fderiv[1](Tdy, Sy, default_sig);
                                    float Twdz = normz * fderiv[2](Tdz, Sz, default_sig);
                                    /* Rotate back towards reference space: wd=R*Twd */
                                    wdx  = Jxx * Twdx + Jxy * Twdy + Jxz * Twdz;
                                    wdy  = Jyx * Twdx + Jyy * Twdy + Jyy * Twdz;
                                    wdz  = Jzx * Twdx + Jzy * Twdy + Jzz * Twdz;
                                }
                                
                                /* Write result in output array */
                                float *pf1k  = pf1;
                                float *pd1xk = pd1x;
                                float *pd1yk = pd1y;
                                float *pd1zk = pd1z;
                                for (k=0; k<n; ++k, pf0x += m0, pf1k += m1, 
                                     pd1xk += m1, pd1yk += m1, pd1zk += m1)
                                {
//                                     if (pf0x-F0 <0 || pf0x-F0 >= m0*n)
//                                         mexErrMsgTxt("Out of bound! (pull acc)");
                                    (*pf1k) += w * (*pf0x);
                                    if (doD1x)
                                        (*pd1xk) += wdx * (*pf0x);
                                    if (doD1y)
                                        (*pd1yk) += wdy * (*pf0x);
                                    if (doD1z)
                                        (*pd1zk) += wdz * (*pf0x);
                                }
                            }
                        }
                    }
                }
            }
            else
            /* Push */
            {
                /* Loop over contributing voxels in ref space */
                for (iz=0; iz<lz;)
                {
                    float dz    = ddz[iz];
                    float *pf0z = F0 + ooz[iz];
                    float *ps0z = S0 + ooz[iz++];
                    for (iy=0; iy<ly;)
                    {
                        float dy    = ddy[iy];
                        float *pf0y = pf0z + ooy[iy];
                        float *ps0y = ps0z + ooy[iy++];
                        for (ix=0; ix<lx;)
                        {
                            float dx    = ddx[ix];
                            float *pf0x = pf0y + oox[ix];
                            float *ps0x = ps0y + oox[ix++];
                            
                            /* Rotate dx towards subject space: Td=R'*d */
                            float Tdx  = Jxx * dx + Jyx * dy + Jzx * dz;
                            float Tdy  = Jxy * dx + Jyy * dy + Jzy * dz;
                            float Tdz  = Jxz * dx + Jyz * dy + Jzz * dz;
                            
                            if (fabs(Tdx)<klimx && fabs(Tdy)<klimy && fabs(Tdz)<klimz) 
                            /* If x is near enough to the basis centre */
                            {
                                /* Compute weight */
                                float w = norm * fwin[0](Tdx, Sx, default_sig) 
                                               * fwin[1](Tdy, Sy, default_sig) 
                                               * fwin[2](Tdz, Sz, default_sig);

                                /* Write result in output array */
                                float *pf1k = pf1;
                                for (k=0; k<n; ++k, pf0x += m0, pf1k += m1)
                                {
//                                     if (pf0x-F0 <0 || pf0x-F0 >= m0*n)
//                                         mexErrMsgTxt("Out of bound! (push F0).");
//                                     if (pf1k-F1 <0 || pf1k-F1 >= m1*n)
//                                         mexErrMsgTxt("Out of bound! (push F1).");
                                    (*pf0x) += w * (*pf1k);
                                }
                                /* Count image */
                                if (S0!=0)
                                {
//                                     if (ps0x-S0 <0 || ps0x-S0 >= m0)
//                                         mexErrMsgTxt("Out of bound! (push S0).");
                                    (*ps0x) += w;
                                }
                            }
                        }
                    }
                }
            }
        }
        else if ((code&1)==1)
        /* If (pull/push AND out of bounds) */
        {
            float *pf1k = pf1;
            for(k=0; k<n; ++k, pf1k += m1)
            {
//                 if (pf1k-F1 <0 || pf1k-F1 >= m1*n)
//                     mexErrMsgTxt("Out of bound! (NaN).");
                (*pf1k) = NaN;
            }
        }
    }
    
    /* Free arrays of function pointers */
    free(fwin);
    free(fnorm);
    free(flim);
    free(fderiv);
}

//-------------------------------------------------------------------------
// Wrappers
//-------------------------------------------------------------------------

static void jpullc(mwSize dm0[], mwSize m1, mwSize n, float Psi[],  float F0[], /*@out@*/ float F1[], float Jac[], mwSize mJ, double win[], float D1x[], float D1y[], float D1z[])
{
    jpushpull(dm0, m1, n, Psi, F0, 0, F1, 3, Jac, mJ, win, D1x, D1y, D1z);
}

static void  jpull(mwSize dm0[], mwSize m1, mwSize n, float Psi[],  float F0[], /*@out@*/ float F1[], float Jac[], mwSize mJ, double win[], float D1x[], float D1y[], float D1z[])
{
    jpushpull(dm0, m1, n, Psi, F0, 0, F1, 1, Jac, mJ, win, D1x, D1y, D1z);
}

static void jpushc(mwSize dm0[], mwSize m1, mwSize n, float Psi[], float F1[], /*@out@*/ float F0[], /*@null@@out@*/ float S0[], float Jac[], mwSize mJ, double win[])
{
    jpushpull(dm0, m1, n, Psi, F0, S0, F1, 2, Jac, mJ, win, (float*)0, (float*)0, (float*)0);
}

static void  jpush(mwSize dm0[], mwSize m1, mwSize n, float Psi[], float F1[], /*@out@*/ float F0[], /*@null@@out@*/ float S0[], float Jac[], mwSize mJ, double win[])
{
    jpushpull(dm0, m1, n, Psi, F0, S0, F1, 0, Jac, mJ, win, (float*)0, (float*)0, (float*)0);
}

//-------------------------------------------------------------------------
// Parse MATLAB arguments
//-------------------------------------------------------------------------

static void jpull_mexFunction(int flag, int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    mwSize nd, i, dm0[4], dm1[4];
    const mwSize *dmpsi, *dmjac;

    if (nrhs == 0) mexErrMsgTxt("Incorrect usage");
    if (nrhs == 4)
    {
        if (nlhs > 4) mexErrMsgTxt("Max 4 output argument required");
    }
    else
        mexErrMsgTxt("Four input arguments required");

    for(i=0; i<3; i++)
        if (!mxIsNumeric(prhs[i]) || mxIsComplex(prhs[i]) || mxIsSparse(prhs[i]) || !mxIsSingle(prhs[i]))
            mexErrMsgTxt("Data must be numeric, real, full and single");
    if (!mxIsNumeric(prhs[3]) || mxIsComplex(prhs[3]) ||
          mxIsSparse(prhs[3]) || !mxIsDouble(prhs[3]))
        mexErrMsgTxt("Window must be numeric, real, full and double");

    /* Dimensions of image to resample */
    nd = mxGetNumberOfDimensions(prhs[0]);
    if (nd>4) mexErrMsgTxt("Wrong number of dimensions.");
    dm0[0] = dm0[1] = dm0[2] = dm0[3] = 1;
    for(i=0; i<nd; i++)
        dm0[i] = mxGetDimensions(prhs[0])[i];

    /* Dimensions of deformation field and resulting output volume */
    nd = mxGetNumberOfDimensions(prhs[1]);
    if (nd!=4) mexErrMsgTxt("Wrong number of dimensions.");
    dmpsi = mxGetDimensions(prhs[1]);
    if (dmpsi[3]!=3)
        mexErrMsgTxt("Incompatible dimensions.");
    
    /* Dimensions of Jacobian field */
    nd = mxGetNumberOfDimensions(prhs[2]);
    if (nd!=5) mexErrMsgTxt("Wrong number of dimensions.");
    dmjac = mxGetDimensions(prhs[2]);
    if ((dmjac[0]!=dmpsi[0] || dmjac[1]!=dmpsi[1] || dmjac[2]!=dmpsi[2]) &&
        (dmjac[0]*dmjac[1]*dmjac[2]>1))
        mexErrMsgTxt("Deformation and Jacobian fields should have the same size.");
    if (dmjac[3]!=3 || dmjac[4]!=3)
        mexErrMsgTxt("Incompatible dimensions.");

    /* Dimensions of output volumes */
    dm1[0] = dmpsi[0];
    dm1[1] = dmpsi[1];
    dm1[2] = dmpsi[2];
    dm1[3] = dm0[3];  /* Number of output volumes */
    plhs[0] = mxCreateNumericArray(4,dm1, mxSINGLE_CLASS, mxREAL);
    float * D1x = (float *)0;
    float * D1y = (float *)0;
    float * D1z = (float *)0;
    if (nlhs >= 2)
    {
        plhs[1] = mxCreateNumericArray(4,dm1, mxSINGLE_CLASS, mxREAL);
        D1x     = (float *)mxGetPr(plhs[1]);
    }
    if (nlhs >= 3)
    {
        plhs[2] = mxCreateNumericArray(4,dm1, mxSINGLE_CLASS, mxREAL);
        D1y     = (float *)mxGetPr(plhs[2]);
    }
    if (nlhs >= 4)
    {
        plhs[3] = mxCreateNumericArray(4,dm1, mxSINGLE_CLASS, mxREAL);
        D1z     = (float *)mxGetPr(plhs[3]);
    }
    
    if ((flag&2)==2)
        jpull (dm0, dm1[0]*dm1[1]*dm1[2], dm0[3], (float *)mxGetPr(prhs[1]), (float *)mxGetPr(prhs[0]), (float *)mxGetPr(plhs[0]), (float *)mxGetPr(prhs[2]), dmjac[0]*dmjac[1]*dmjac[2], (double *)mxGetPr(prhs[3]), D1x, D1y, D1z);
    else
        jpullc(dm0, dm1[0]*dm1[1]*dm1[2], dm0[3], (float *)mxGetPr(prhs[1]), (float *)mxGetPr(prhs[0]), (float *)mxGetPr(plhs[0]), (float *)mxGetPr(prhs[2]), dmjac[0]*dmjac[1]*dmjac[2], (double *)mxGetPr(prhs[3]), D1x, D1y, D1z);
}

static void jpush_mexFunction(int flag, int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    float *S0;
    mwSize nd, i, dm0[4], dm1[4];
    const mwSize *dmpsi, *dmjac;

    if ((nrhs != 4) && (nrhs != 5))
        mexErrMsgTxt("Four or five input arguments required");
    if (nlhs  > 2) mexErrMsgTxt("Up to two output arguments required");
    
    for(i=0; i<3; i++)
        if (!mxIsNumeric(prhs[i]) || mxIsComplex(prhs[i]) ||
              mxIsSparse(prhs[i]) || !mxIsSingle(prhs[i]))
            mexErrMsgTxt("Data must be numeric, real, full and single");
    if (!mxIsNumeric(prhs[3]) || mxIsComplex(prhs[3]) ||
          mxIsSparse(prhs[3]) || !mxIsDouble(prhs[3]))
        mexErrMsgTxt("Window must be numeric, real, full and double");

    /* Dimensions of input volumes */
    nd = mxGetNumberOfDimensions(prhs[0]);
    if (nd>4) mexErrMsgTxt("Wrong number of dimensions.");
    dm1[0] = dm1[1] = dm1[2] = dm1[3] = 1;
    for(i=0; i<nd; i++)
        dm1[i] = mxGetDimensions(prhs[0])[i];

    /* Dimensions of deformation */
    nd = mxGetNumberOfDimensions(prhs[1]);
    if (nd!=4) mexErrMsgTxt("Wrong number of dimensions.");
    dmpsi = mxGetDimensions(prhs[1]);
    if (dmpsi[0]!=dm1[0] || dmpsi[1]!=dm1[1] || dmpsi[2]!=dm1[2] || dmpsi[3]!=3)
        mexErrMsgTxt("Incompatible dimensions.");
    
    /* Dimensions of Jacobian field */
    nd = mxGetNumberOfDimensions(prhs[2]);
    if (nd!=5) mexErrMsgTxt("Wrong number of dimensions.");
    dmjac = mxGetDimensions(prhs[2]);
    if ((dmjac[0]!=dmpsi[0] || dmjac[1]!=dmpsi[1] || dmjac[2]!=dmpsi[2]) &&
        (dmjac[0]*dmjac[1]*dmjac[2]>1))
        mexErrMsgTxt("Deformation and Jacobian fields should have the same size.");
    if (dmjac[3]!=3 || dmjac[4]!=3)
        mexErrMsgTxt("Incompatible dimensions.");
    
    /* Dimensions of output volumes */
    if (nrhs>=5)
    {
        if (!mxIsNumeric(prhs[4]) || mxIsComplex(prhs[4]) ||
              mxIsSparse(prhs[4]) || !mxIsDouble(prhs[4]))
            mexErrMsgTxt("Data must be numeric, real, full and double");
        if (mxGetNumberOfElements(prhs[4])!= 3)
            mexErrMsgTxt("Output dimensions must have three elements");
        dm0[0] = (mwSize)floor((double)mxGetPr(prhs[4])[0]);
        dm0[1] = (mwSize)floor((double)mxGetPr(prhs[4])[1]);
        dm0[2] = (mwSize)floor((double)mxGetPr(prhs[4])[2]);
    }
    else
    {
        dm0[0] = dm1[0];
        dm0[1] = dm1[1];
        dm0[2] = dm1[2];
    }
    dm0[3]  = dm1[3];
    plhs[0] = mxCreateNumericArray(4,dm0, mxSINGLE_CLASS, mxREAL);

 
    /* Create count volume if required */
    if (nlhs>=2)
    {
        plhs[1] = mxCreateNumericArray(3,dm0, mxSINGLE_CLASS, mxREAL);
        S0      = (float *)mxGetPr(plhs[1]);
    }
    else
        S0      = (float *)0;
 
    if ((flag&2)==2)
        jpush (dm0, dm1[0]*dm1[1]*dm1[2], dm0[3], (float *)mxGetPr(prhs[1]), (float *)mxGetPr(prhs[0]), (float *)mxGetPr(plhs[0]), S0, (float *)mxGetPr(prhs[2]), dmjac[0]*dmjac[1]*dmjac[2], (double *)mxGetPr(prhs[3]));
    else
        jpushc(dm0, dm1[0]*dm1[1]*dm1[2], dm0[3], (float *)mxGetPr(prhs[1]), (float *)mxGetPr(prhs[0]), (float *)mxGetPr(plhs[0]), S0, (float *)mxGetPr(prhs[2]), dmjac[0]*dmjac[1]*dmjac[2], (double *)mxGetPr(prhs[3]));

}

static void boundary_mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    if ((nlhs<=1) && (nrhs==0))
    {
        mwSize nout[] = {1,1,1};
        plhs[0] = mxCreateNumericArray(2,nout, mxDOUBLE_CLASS, mxREAL);
        mxGetPr(plhs[0])[0] = get_bound();
    }
    else if ((nrhs==1) && (nlhs==0))
    {
        if (!mxIsNumeric(prhs[0]) || mxIsComplex(prhs[0]) || mxIsSparse(prhs[0]) || !mxIsDouble(prhs[0]))
            mexErrMsgTxt("Data must be numeric, real, full and double");
        set_bound(mxGetPr(prhs[0])[0]);
    }
}

//-------------------------------------------------------------------------
// Main MEX
//-------------------------------------------------------------------------

#include<string.h>

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    set_bound(get_bound());

    if ((nrhs>=1) && mxIsChar(prhs[0]))
    {
        int buflen;
        char *fnc_str;
        buflen = mxGetNumberOfElements(prhs[0]);
        fnc_str = (char *)mxCalloc(buflen+1,sizeof(mxChar));
        mxGetString(prhs[0],fnc_str,buflen+1);

        if (!strcmp(fnc_str,"pull"))
        {
            mxFree(fnc_str);
            jpull_mexFunction(2, nlhs, plhs, nrhs-1, &prhs[1]);
        }
        else if (!strcmp(fnc_str,"pullc"))
        {
            mxFree(fnc_str);
            jpull_mexFunction(0, nlhs, plhs, nrhs-1, &prhs[1]);
        }
        else if (!strcmp(fnc_str,"push"))
        {
            mxFree(fnc_str);
            jpush_mexFunction(2, nlhs, plhs, nrhs-1, &prhs[1]);
        }
        else if (!strcmp(fnc_str,"pushc"))
        {
            mxFree(fnc_str);
            jpush_mexFunction(0, nlhs, plhs, nrhs-1, &prhs[1]);
        }
        else if (!strcmp(fnc_str,"boundary")  || !strcmp(fnc_str,"bound"))
        {
            mxFree(fnc_str);
            boundary_mexFunction(nlhs, plhs, nrhs-1, &prhs[1]);
        }
        else
        {
            mxFree(fnc_str);
            mexErrMsgTxt("Option not recognised.");
        }
    }
    else
    {
        mexErrMsgTxt("Option not recognised.");
    }
}