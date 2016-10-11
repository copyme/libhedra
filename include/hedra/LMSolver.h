// This file is part of libhedra, a library for polyhedral mesh processing
//
// Copyright (C) 2016 Amir Vaxman <avaxman@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
#ifndef HEDRA_LEVENBERG_MARQUADT_SOLVER_H
#define HEDRA_LEVENBERG_MARQUADT_SOLVER_H
#include <igl/igl_inline.h>
#include <Eigen/Core>
#include <string>
#include <vector>
#include <cstdio>
#include <iostream>

namespace hedra {
    namespace optimization
    {
        
        template<class LinearSolver, class SolverTraits>
        class LMSolver{
        public:
            Eigen::VectorXd x;      //current solution; always updated
            Eigen::VectorXd prevx;  //the solution of the previous iteration
            Eigen::VectorXd x0;     //the initial solution to the system
            Eigen::VectorXd d;             //the direction taken.
            Eigen::VectorXd currEnergy;    //the current value of the energy
            Eigen::VectorXd prevEnergy;    //the previous value of the energy
            
            Eigen::VectorXi HRows, HCols;  //(row,col) pairs for H=J^T*J matrix
            Eigen::VectorXd HVals;      //values for H matrix
            Eigen::MatrixXi S2D;        //single J to J^J indices
            
            LinearSolver* LS;
            SolverTraits* ST;
            int maxIterations;
            double xTolerance;
            double fooTolerance;
            
            //Input: pattern of matrix M by (iI,iJ) representation
            //Output: pattern of matrix M^T*M by (oI, oJ) representation
            //        map between values in the input to values in the output (Single2Double). The map is aggregating values from future iS to oS
            //prerequisite: iI are sorted by rows (not necessary columns)
            void MatrixPattern(const Eigen::VectorXi& iI,
                               const Eigen::VectorXi& iJ,
                               Eigen::VectorXi& oI,
                               Eigen::VectorXi& oJ,
                               Eigen::MatrixXi& S2D)
            {
                int CurrTri=0;
                using namespace Eigen;
                std::vector<int> oIlist;
                std::vector<int> oJlist;
                std::vector<std::pair<int, int> > S2Dlist;
                do{
                    int CurrRow=iI(CurrTri);
                    int NumCurrTris=0;
                    while ((CurrTri+NumCurrTris<iI.size())&&(iI(CurrTri+NumCurrTris)==CurrRow))
                        NumCurrTris++;
                    
                    for (int i=CurrTri;i<CurrTri+NumCurrTris;i++){
                        for (int j=CurrTri;j<CurrTri+NumCurrTris;j++){
                            if (iJ(j)>=iJ(i)){
                                oIlist.push_back(iJ(i));
                                oJlist.push_back(iJ(j));
                                S2Dlist.push_back(std::pair<int,int>(i,j));
                            }
                        }
                    }
                    CurrTri+=NumCurrTris;
                    //std::cout<<"CurrTri, NumCurrTris: "<<CurrTri<<","<<NumCurrTris<<std::endl;
                }while (CurrTri!=iI.size());
                
                oI.resize(oIlist.size());
                oJ.resize(oJlist.size());
                S2D.resize(S2Dlist.size(),2);
                
                for (int i=0;i<oIlist.size();i++){
                    oI(i)=oIlist[i];
                    oJ(i)=oJlist[i];
                    S2D.row(i)<<S2Dlist[i].first, S2Dlist[i].second;
                }
            }
            
            //returns the values of M^T*M+miu*I by multiplication and aggregating from Single2double list.
            //prerequisite - oS is allocated
            void MatrixValues(const Eigen::VectorXi& oI,
                              const Eigen::VectorXi& oJ,
                              const Eigen::VectorXd& iS,
                              const Eigen::MatrixXi& S2D,
                              Eigen::VectorXd& oS,
                              double miu)
            {
                for (int i=0;i<S2D.rows();i++)
                    oS(i)=iS(S2D(i,0))*iS(S2D(i,1));
                
                //adding miu*I
                for (int i=0;i<S2D.rows();i++)
                    if (oI(i)==oJ(i))
                        oS(i)+=miu;
            }
            
            //returns M^t*ivec by (I,J,S) representation
            void MultiplyAdjointVector(const Eigen::VectorXi& iI,
                                       const Eigen::VectorXi& iJ,
                                       const Eigen::VectorXd& iS,
                                       const Eigen::VectorXd& iVec,
                                       Eigen::VectorXd& oVec)
            {
                oVec.setZero();
                for (int i=0;i<iI.size();i++)
                    oVec(iJ(i))+=iS(i)*iVec(iI(i));
            }
            
            
        public:
            
            LMSolver(){};
            
            void init(LinearSolver* _LS,
                      SolverTraits* _ST,
                      int _maxIterations=100,
                      double _xTolerance=10e-9,
                      double _fooTolerance=10e-9){
                
                LS=_LS;
                ST=_ST;
                maxIterations=_maxIterations;
                xTolerance=_xTolerance;
                fooTolerance=_fooTolerance;
                //analysing pattern
                MatrixPattern(ST->JRows, ST->JCols,HRows,HCols,S2D);
                HVals.resize(HRows.size());
                
                LS->analyze(HRows,HCols);
                
                d.resize(ST->xSize);
                x.resize(ST->xSize);
                x0.resize(ST->xSize);
                prevx.resize(ST->xSize);
                currEnergy.resize(ST->EVec.size());
                prevEnergy.resize(ST->EVec.size());
            }
            
            
            bool solve(const bool verbose) {
                
                using namespace Eigen;
                using namespace std;
                ST->initial_solution(x0);
                prevx<<x0;
                int currIter=0;
                bool stop=false;
                double currError, prevError;
                VectorXd rhs(ST->xSize);
                VectorXd direction;
                if (verbose)
                    cout<<"******Beginning Optimization******"<<endl;
                
                double tau=10e-3;
                
                //estimating initial miu
                double miu=0.0;
                ST->update_jacobian(prevx);
                MatrixValues(HRows, HCols, ST->JVals, S2D, HVals, miu);
                for (int i=0;i<HRows.size();i++)
                    if (HRows(i)==HCols(i))  //on the diagonal
                        miu=(miu < HVals(i) ? HVals(i) : miu);
                miu*=tau;
                cout<<"initial miu: "<<miu<<endl;
                double beta=2.0;
                double nu=beta;
                double gamma=3.0;
                do{
                    currIter=0;
                    stop=false;
                    do{
                        ST->pre_iteration(prevx);
                        ST->update_energy(prevx);
                        ST->update_jacobian(prevx);
                        if (verbose)
                            cout<<"Initial Energy for Iteration "<<currIter<<": "<<ST->EVec.template squaredNorm()<<endl;
                        MatrixValues(HRows, HCols, ST->JVals, S2D, HVals, miu);
                        MultiplyAdjointVector(ST->JRows, ST->JCols, ST->JVals, -ST->EVec, rhs);
                        
                        double firstOrderOptimality=rhs.template lpNorm<Infinity>();
                        cout<<"firstOrderOptimality: "<<firstOrderOptimality<<endl;
                        
                        if (firstOrderOptimality<fooTolerance){
                            if (verbose){
                                cout<<"First-order optimality has been reached"<<endl;
                                break;
                            }
                        }
                        
                        //solving to get the GN direction
                        if(!LS->factorize(HVals)) {
                            // decomposition failed
                            cout<<"Solver Failed to factorize! "<<endl;
                            return false;
                        }
                        
                        LS->solve(rhs,direction);
                        cout<<"direction magnitude: "<<direction.norm()<<endl;
                        if (direction.norm() < xTolerance * prevx.norm()){
                            if (verbose)
                                cout<<"Stopping since direction magnitude small."<<endl;
                            break;
                        }
                        VectorXd tryx=prevx+direction;
                        ST->update_energy(prevx);
                        double prevE=ST->EVec.squaredNorm();
                        ST->update_energy(tryx);
                        double currE=ST->EVec.squaredNorm();
                        
                        double rho=(prevE-currE)/(direction.dot(miu*direction+rhs));
                        if (rho>0){
                            x=tryx;
                            cout<<"Energy: "<<currE<<endl;
                            cout<<"1.0-pow(2.0*rho-1.0,3.0): "<<1.0-pow(2.0*rho-1.0,3.0)<<endl;
                            miu*=(1.0/gamma > 1.0-(beta-1.0)*pow(2.0*rho-1.0,3) ? 1.0/gamma : 1.0-(beta-1.0)*pow(2.0*rho-1.0,3));
                            nu=beta;
                            cout<<"rho, miu, nu: "<<rho<<","<<miu<<","<<nu<<endl;
                        } else {
                            miu = miu*nu;
                            nu=2*nu;
                            cout<<"rho, miu, nu: "<<rho<<","<<miu<<","<<nu<<endl;
                        }
                        
                        //The SolverTraits can order the optimization to stop by giving "true" of to continue by giving "false"
                        if (ST->post_iteration(x)){
                            if (verbose)
                                cout<<"ST->Post_iteration() gave a stop"<<endl;
                            break;
                        }
                        currIter++;
                        prevx=x;
                    }while (currIter<=maxIterations);
                }while (!ST->post_optimization(x));
                return true;
            }
        };
        
    }
}


#endif
