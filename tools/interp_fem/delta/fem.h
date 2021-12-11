#ifndef FEM
#define FEM

#include <iostream>

class Fem{

private:
    int numgp{11};
    int numparams{5};

public:
    double* gp;


    Fem(){gp = gausspoints();}


    ~Fem(){
        delete[] gp;
        gp  = NULL;
    }


    double* gausspoints(){
        gp = new double[numgp*numparams];

        for (int i=0; i<3; ++i){
            for (int j=0; j<10; j++){
                gp[numgp*i+j] = 0.0;
            }

            gp[numgp*i+i]   = 1.0;
            gp[numgp*i+i+4] = 0.5;
            gp[numgp*i+i+7] = 0.5;
            gp[numgp*i+10]  = 0.25;
        }

        gp[numgp*1+4] = 0.5;
        gp[numgp*2+5] = 0.5;
        gp[numgp*0+9] = 0.5;

        for (int j=0; j<4; ++j){
            gp[4*numgp+j] = 1.0/360.0;
        }

        for (int j=4; j<10; ++j){
            gp[4*numgp+j] = 1.0/90.0;
        }

        gp[4*numgp+10] = 4.0/45.0;
        gp[0*numgp+6]  = 0.5;
        gp[0*numgp+9]  = 0.0;

        //compute fourth point coordinates
        for (int j=0; j<numgp; ++j){
            gp[3*numgp+j] = 1.0 - (gp[0*numgp+j] + gp[1*numgp+j] + gp[2*numgp+j]);
        }

        return gp;
    }


    double tetshp(int p, double* xl){

        double shp[4*4];

        shp[3*4+0] = gp[0*numgp+p];
        shp[3*4+1] = gp[1*numgp+p];
        shp[3*4+2] = gp[2*numgp+p];
        shp[3*4+3] = gp[3*numgp+p];

        double shpxi[3*4];

        shpxi[0*4+0] =  1.0;
        shpxi[0*4+1] =  0.0;
        shpxi[0*4+2] =  0.0;
        shpxi[0*4+3] = -1.0;

        shpxi[1*4+0] =  0.0;
        shpxi[1*4+1] =  1.0;
        shpxi[1*4+2] =  0.0;
        shpxi[1*4+3] = -1.0;

        shpxi[2*4+0] =  0.0;
        shpxi[2*4+1] =  0.0;
        shpxi[2*4+2] =  1.0;
        shpxi[2*4+3] = -1.0;

        double xs[3*3];

        for (int i=0; i<3; ++i){
            for (int j=0; j<3; ++j){
                xs[3*i+j] = 0.0;
                for (int k=0; k<4; ++k){
                    xs[3*i+j] -= shpxi[4*i+k] * xl[4*j+k];
                }
            }
        }

        double xsj = xs[0] * (xs[4]*xs[8] - xs[5]*xs[7]) + xs[1] * (xs[5]*xs[6] - xs[3]*xs[8]) + xs[2] * (xs[3]*xs[7] - xs[4]*xs[6]);

        double xsinv[3*3];

        xsinv[0] = (xs[4]*xs[8] - xs[7]*xs[5]) / xsj;
        xsinv[3] = (xs[7]*xs[2] - xs[1]*xs[8]) / xsj;
        xsinv[6] = (xs[1]*xs[5] - xs[4]*xs[2]) / xsj;

        xsinv[1] = (xs[6]*xs[5] - xs[3]*xs[8]) / xsj;
        xsinv[4] = (xs[0]*xs[8] - xs[6]*xs[2]) / xsj;
        xsinv[7] = (xs[3]*xs[2] - xs[0]*xs[5]) / xsj;

        xsinv[2] = (xs[3]*xs[7] - xs[6]*xs[4]) / xsj;
        xsinv[5] = (xs[6]*xs[1] - xs[0]*xs[7]) / xsj;
        xsinv[8] = (xs[0]*xs[4] - xs[3]*xs[1]) / xsj;

        for (int k=0; k<3; ++k){
            for (int j=0; j<2; ++j){
                for (int i=0; i<2; ++i){
                    shp[4*j+k] += shpxi[4*i+k]*xsinv[3*i+j];
                }
            }
        }

        return xsj;
    }
};

#endif
