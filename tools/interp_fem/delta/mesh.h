#ifndef MESH
#define MESH

#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <vector>

std::vector<std::string> tokenize(std::string input_string){
    std::string buf;
    std::stringstream ss(input_string);
    std::vector<std::string> tokens;

    while(ss >> buf)
        tokens.push_back(buf);

    return tokens;
}

class Mesh{

private:
    int numel;
    int numnp;
    int nen{4};
    int ndm{3};
    int* ix;
    double* xc;

public:
    Mesh(const int& num1, const int& num2)
    : numel(num1), numnp(num2){

        ix = readElemCon();
        xc = readNodeCoords();
    }


    ~Mesh(){
        delete[] ix;
        delete[] xc;

        ix = NULL;
        xc = NULL;
    }


    int* readElemCon(){
        ix = new int[nen*numel];

        std::ifstream elemcon;
        elemcon.open("elemcon.txt");

        std::vector<std::string> tokens;

        for (int i=0; i<numel; ++i){
            std::string current_line;
            getline(elemcon, current_line);
            tokens = tokenize(current_line);

            for (int j=0; j<nen; ++j){
                ix[i*nen+j] = atoi(tokens[j].c_str());
            }
        }

        elemcon.close();
        return ix;
    }


    double* readNodeCoords(){
        xc = new double[ndm*numnp];

        std::ifstream meshpoints;
        meshpoints.open("meshpoints.txt");

        std::vector<std::string> tokens;

        for (int i=0; i<numnp; ++i){
            std::string current_line;
            getline(meshpoints, current_line);
            tokens = tokenize(current_line);

            for (int j=0; j<ndm; ++j){
                xc[i*ndm+j] = atof(tokens[j].c_str());
            }
        }

        meshpoints.close();
        return xc;
    }


    void elementsContainingNode(const int& gid){
        for (int i=0; i<numel; ++i){
            for (int j=0; j<nen; ++j){
                if (ix[i*nen+j] == gid){
                    std::cout << "Node " << gid << " is in element " << i << " with volume " << computeElemVolume(i) << '\n';
                }
            }
        }
    }


    double computeElemVolume(const int& elemId){
        double* xl = new double[ndm*nen];

        Fem fem;

        for (int i=0; i<nen; ++i){
            int ii = ix[elemId*nen+i];
            for (int j=0; j<ndm; ++j){
                xl[nen*j+i] = xc[ndm*ii+j];
            }
        }

        double volel = 0.0;
        double xsj;

        for (int kk=0; kk<11; ++kk){
            xsj = fem.tetshp(kk, xl);
            xsj *= fem.gp[11*4+kk];

            volel += xsj;
        }

        delete[] xl;
        xl = NULL;

        return volel;
    }


    double computeMeshVolume(){
        double* xl  = new double[ndm*nen];

        Fem fem;

        double vol = 0.0;

        for (int n=0; n<numel; ++n){
            for (int i=0; i<nen; ++i){
                int ii = ix[n*nen+i];
                for (int j=0; j<ndm; ++j){
                    xl[nen*j+i] = xc[ndm*ii+j];
                }
            }

            double volel = 0.0;
            double xsj;

            for (int kk=0; kk<11; ++kk){
                xsj = fem.tetshp(kk, xl);
                xsj *= fem.gp[11*4+kk];

                volel += xsj;
            }

            std::cout << volel << '\n';

            vol += volel;
        }

        delete[] xl;
        xl  = NULL;

        return vol;
    }

    friend std::vector<std::string> tokenize(std::string input_string);
};

#endif
