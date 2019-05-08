#include "mex.h"
#include <math.h>
#include <string.h>

#define MAXCLASSES 1024

/* mex_mrf("apply",uint8(Z),single(lnG),single(w)) */
void apply(mwSize dm[], unsigned char oZ[], float lnG[], float w[], float nZ[])
{
    mwSize i0, i1, i2, k, m, n;
    float a[MAXCLASSES], sa[MAXCLASSES];
    unsigned char *oz0 = NULL, *oz1 = NULL;
    float *nz0 = NULL, *nz1 = NULL;
    int it;

    m = dm[0]*dm[1]*dm[2];        
    
    for(it=0; it<2; it++) 
    {
        mwSize i2start = it%2;
        for(i2=0; i2<dm[2]; i2++) /* Inferior -> Superior */
        {
            mwSize i1start = (i2start == (i2%2));
            for(i1=0; i1<dm[1]; i1++) /* Posterior -> Anterior */
            {
                mwSize i0start = (i1start == (i1%2));                
                oz1 = oZ + dm[0]*(i1+dm[1]*i2);
                nz1 = nZ  + dm[0]*(i1+dm[1]*i2);
                
                for(i0=i0start; i0<dm[0]; i0+=2) /* Left -> Right */
                {
                    unsigned char *qq = NULL;

                    /* Pointers to current voxel in first volume */                    
                    oz0 = oz1 + i0;
                    nz0 = nz1 + i0;

                    /* Initialise neighbour counts to zero */
                    for(k=0; k<dm[3]; k++) a[k] = 0.0;

                    /* Count neighbours of each class */
                    if(i2>0)       /* Inferior */
                    {

                        qq = oz0 - dm[0]*dm[1];
                        for(k=0; k<dm[3]; k++) a[k] += qq[k*m]*w[2];
                    }

                    if(i2<dm[2]-1) /* Superior */
                    {
                        qq = oz0 + dm[0]*dm[1];
                        for(k=0; k<dm[3]; k++) a[k] += qq[k*m]*w[2];
                    }

                    if(i1>0)       /* Posterior */
                    {
                        qq = oz0 - dm[0];
                        for(k=0; k<dm[3]; k++) a[k] += qq[k*m]*w[1];
                    }

                    if(i1<dm[1]-1) /* Anterior */
                    {
                        qq = oz0 + dm[0];
                        for(k=0; k<dm[3]; k++) a[k] += qq[k*m]*w[1];
                    }

                    if(i0>0)       /* Left */
                    {
                        qq = oz0 - 1;
                        for(k=0; k<dm[3]; k++) a[k] += qq[k*m]*w[0];
                    }

                    if(i0<dm[0]-1) /* Right */
                    {
                        qq = oz0 + 1;
                        for(k=0; k<dm[3]; k++) a[k] += qq[k*m]*w[0];
                    }

                    /* Responsibility data is uint8 */
                    for(k=0; k<dm[3]; k++)
                        a[k]/=(255.0*6.0);

                    /* Weights are in the form of a matrix,
                       shared among all voxels. */
                    float *g;
                    for(k=0, g=lnG; k<dm[3]; k++)
                    {
                        sa[k] = 0;
                        for(n=0; n<dm[3]; n++, g++)
                            sa[k] += (*g)*a[n];
                    }
                    
                    for(k=0; k<dm[3]; k++)
                        nz0[k*m] = sa[k];
                }
            }
        }
    }
}

/* mex_mrf("lowerbound",uint8(nZ),single(lnG),single(w)) */
void lowerbound(mwSize dm[], unsigned char oZ[], float lnG[], float w[], double lb[])
{
    mwSize i0, i1, i2, k, m, n;
    float a[MAXCLASSES], sa[MAXCLASSES];
    unsigned char *oz0 = NULL, *oz1 = NULL;
    int it;
    double r = 0;
    
    m = dm[0]*dm[1]*dm[2];        
    
    /* Checkerboard loop (only over one color --- it<1) */
    for(it=0; it<1; it++) 
    {
        mwSize i2start = it%2;
        for(i2=0; i2<dm[2]; i2++) /* Inferior -> Superior */
        {
            mwSize i1start = (i2start == (i2%2));
            for(i1=0; i1<dm[1]; i1++) /* Posterior -> Anterior */
            {
                mwSize i0start = (i1start == (i1%2));                
                oz1 = oZ + dm[0]*(i1+dm[1]*i2);
                
                for(i0=i0start; i0<dm[0]; i0+=2) /* Left -> Right */
                {
                    unsigned char *qq = NULL;

                    /* Pointers to current voxel in first volume */                    
                    oz0 = oz1 + i0;

                    /* Initialise neighbour counts to zero */
                    for(k=0; k<dm[3]; k++) a[k] = 0.0;

                    /* Count neighbours of each class */
                    if(i2>0)       /* Inferior */
                    {

                        qq = oz0 - dm[0]*dm[1];
                        for(k=0; k<dm[3]; k++) a[k] += qq[k*m]*w[2];
                    }

                    if(i2<dm[2]-1) /* Superior */
                    {
                        qq = oz0 + dm[0]*dm[1];
                        for(k=0; k<dm[3]; k++) a[k] += qq[k*m]*w[2];
                    }

                    if(i1>0)       /* Posterior */
                    {
                        qq = oz0 - dm[0];
                        for(k=0; k<dm[3]; k++) a[k] += qq[k*m]*w[1];
                    }

                    if(i1<dm[1]-1) /* Anterior */
                    {
                        qq = oz0 + dm[0];
                        for(k=0; k<dm[3]; k++) a[k] += qq[k*m]*w[1];
                    }

                    if(i0>0)       /* Left */
                    {
                        qq = oz0 - 1;
                        for(k=0; k<dm[3]; k++) a[k] += qq[k*m]*w[0];
                    }

                    if(i0<dm[0]-1) /* Right */
                    {
                        qq = oz0 + 1;
                        for(k=0; k<dm[3]; k++) a[k] += qq[k*m]*w[0];
                    }

                    /* Responsibility data is uint8 */
                    for(k=0; k<dm[3]; k++)
                        a[k]/=(255.0*6.0);

                    /* Weights are in the form of a matrix,
                       shared among all voxels. */
                    float *g;
                    for(k=0, g=lnG; k<dm[3]; k++)
                    {
                        sa[k] = 0;
                        for(n=0; n<dm[3]; n++, g++)
                            sa[k] += (*g)*a[n];
                    }
                                        
                    for(k=0; k<dm[3]; k++)
                    {
                        r     = (double)oz0[k*m];
                        r     /= 255;
                        lb[0] += r*(double)sa[k];
                    }
                }
            }
        }
    }
}

/* mex_mrf("update",uint8(nZ),single(G),single(w)) */
void update(mwSize dm[], unsigned char oZ[], float w[], double nG[])
{
    mwSize i0, i1, i2, k, m, n, k1, k2;        
    int it, K;
    double r = 0;
    double s = 0;
    
    K = dm[3];
    m = dm[0]*dm[1]*dm[2];        
        
    for(k1=0; k1<K; k1++) 
    {
        for(k2=0; k2<K; k2++) 
        {
            unsigned char *oz0 = NULL, *oz1 = NULL;
            
            /* Checkerboard loop (only over one color --- it<1) */
            for(it=0; it<1; it++) 
            {    
                mwSize i2start = it%2;
                for(i2=0; i2<dm[2]; i2++) /* Inferior -> Superior */
                {
                    mwSize i1start = (i2start == (i2%2));
                    for(i1=0; i1<dm[1]; i1++) /* Posterior -> Anterior */
                    {
                        mwSize i0start = (i1start == (i1%2));                
                        oz1 = oZ + dm[0]*(i1+dm[1]*i2);

                        for(i0=i0start; i0<dm[0]; i0+=2) /* Left -> Right */
                        {
                            unsigned char *qq = NULL;

                            /* Pointers to current voxel in first volume */                    
                            oz0 = oz1 + i0;
                            s   = 0;
                            
                            /* Count neighbours of each class */
                            if(i2>0)       /* Inferior */
                            {
                                qq = oz0 - dm[0]*dm[1];
                                s += qq[k2*m]*w[2];
                            }

                            if(i2<dm[2]-1) /* Superior */
                            {
                                qq = oz0 + dm[0]*dm[1];
                                s += qq[k2*m]*w[2];
                            }

                            if(i1>0)       /* Posterior */
                            {
                                qq = oz0 - dm[0];
                                s += qq[k2*m]*w[1];
                            }

                            if(i1<dm[1]-1) /* Anterior */
                            {
                                qq = oz0 + dm[0];
                                s += qq[k2*m]*w[1];
                            }

                            if(i0>0)       /* Left */
                            {
                                qq = oz0 - 1;
                                s += qq[k2*m]*w[0];
                            }

                            if(i0<dm[0]-1) /* Right */
                            {
                                qq = oz0 + 1;
                                s += qq[k2*m]*w[0];
                            }
                                                        
                            s /=(255.0*6.0);

                            r = (double)oz0[k1*m];
                            r /= 255;
                            
                            s *= r;
                                                                                    
                            nG[K*k1 + k2] += s;
                        }
                    }
                }
            }
        }
    }
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    if ((nrhs>=1) && mxIsChar(prhs[0]))
    {
        int buflen;
        char *fnc_str;
        buflen = mxGetNumberOfElements(prhs[0]);
        fnc_str = (char *)mxCalloc(buflen+1,sizeof(mxChar));
        mxGetString(prhs[0],fnc_str,buflen+1);

        if (!strcmp(fnc_str,"apply") || !strcmp(fnc_str,"lowerbound") || !strcmp(fnc_str,"update"))
        {	 
            mwSize i;
            
            /* Parse input */
            if (nrhs<2 || nrhs>4)
                mexErrMsgTxt("Incorrect usage");
            
            if (!mxIsNumeric(prhs[1]) || mxIsComplex(prhs[1]) || mxIsSparse(prhs[1]) || !mxIsUint8(prhs[1]))
                mexErrMsgTxt("Second arg must be numeric, real, full and uint8.");

            /* Resonsibilities */
            unsigned char *oZ = NULL;
            
            oZ = (unsigned char *)mxGetData(prhs[1]);
            
            /* Get dimensions */
            mwSize dm[4];
            
            for(i=0; i<mxGetNumberOfDimensions(prhs[1]); i++)
                dm[i] = mxGetDimensions(prhs[1])[i];

            for(i=mxGetNumberOfDimensions(prhs[1]); i<4; i++)
                dm[i] = 1;                                              
                                                      
            if (!strcmp(fnc_str,"apply") || !strcmp(fnc_str,"lowerbound"))
            {
                /* Apply and lower bound */
                
                /* Log of weight matrix */
                if (!mxIsNumeric(prhs[2]) || mxIsComplex(prhs[2]) || mxIsSparse(prhs[2]) || !mxIsSingle(prhs[2]))
                    mexErrMsgTxt("Third arg must be numeric, real, full and single.");  

                float *lnG = NULL;

                lnG = (float *)mxGetData(prhs[2]);

                /* Classes */
                mwSize K[4];

                for(i=0; i<mxGetNumberOfDimensions(prhs[2]); i++)
                    K[i] = mxGetDimensions(prhs[2])[i];

                /* Sanity check */
                if (mxGetNumberOfDimensions(prhs[2]) > 2 || K[0] != K[1] || K[0] != dm[3])
                    mexErrMsgTxt("Error with number of classes");

                /* Get square of each voxel size */
                float w[3];

                if (nrhs==4)
                {
                    if (!mxIsNumeric(prhs[3]) || mxIsComplex(prhs[3]) || mxIsSparse(prhs[3]) || !mxIsSingle(prhs[3]))
                        mexErrMsgTxt("Fourth arg must be numeric, real, full and single.");

                    if (mxGetNumberOfElements(prhs[3]) != 3)
                        mexErrMsgTxt("Fourth arg must contain three elements.");

                    for(i=0; i<3; i++) w[i] = ((float *)mxGetData(prhs[3]))[i];
                }
                else
                {
                    for(i=0; i<3; i++) w[i] = 1.0;
                }
            
                if (!strcmp(fnc_str,"apply"))
                {
                    /* Allocate output */
                    float *nZ;
                    plhs[0] = mxCreateNumericArray(4,dm, mxSINGLE_CLASS, mxREAL);
                    nZ       = (float *)mxGetData(plhs[0]);        

                    /* Call function */
                    apply(dm,oZ,lnG,w,nZ);
                }
                else if (!strcmp(fnc_str,"lowerbound"))
                {
                    /* Allocate output */
                    double *lb;
                    plhs[0] = mxCreateDoubleScalar(0);
                    lb      = (double *)mxGetData(plhs[0]);  

                    /* Call function */
                    lowerbound(dm,oZ,lnG,w,lb);
                }
            }
            else if (!strcmp(fnc_str,"update"))
            {
                /* Update */
                
                /* Get square of each voxel size */
                float w[3];

                if (nrhs==4)
                {
                    if (!mxIsNumeric(prhs[2]) || mxIsComplex(prhs[2]) || mxIsSparse(prhs[2]) || !mxIsSingle(prhs[2]))
                        mexErrMsgTxt("Fourth arg must be numeric, real, full and single.");

                    if (mxGetNumberOfElements(prhs[2]) != 3)
                        mexErrMsgTxt("Fourth arg must contain three elements.");

                    for(i=0; i<3; i++) w[i] = ((float *)mxGetData(prhs[2]))[i];
                }
                else
                {
                    for(i=0; i<3; i++) w[i] = 1.0;
                }
            
                /* Allocate output */
                double *nG;
                plhs[0] = mxCreateNumericMatrix(dm[3],dm[3], mxDOUBLE_CLASS, mxREAL);
                nG      = (double *)mxGetData(plhs[0]);  
                 
                /* Call function */
                update(dm,oZ,w,nG);
            }
        }
        else
        {
            mexErrMsgTxt("Option not recognised.");
        }

        mxFree(fnc_str);
    }
    else
    {
        mexErrMsgTxt("Incorrect usage.");
    }
}
