// This file is part of libhedra, a library for polyhedral mesh processing //
// Copyright (C) 2019 Amir Vaxman <avaxman@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
#ifndef HEDRA_COPYLEFT_CGAL_EXTRACT_MESH_H
#define HEDRA_COPYLEFT_CGAL_EXTRACT_MESH_H

#include <stdexcept>
#include <vector>
#include <set>


#include <boost/config.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/connected_components.hpp>

#include <igl/igl_inline.h>
#include <Eigen/Dense>
#include <CGAL/basic.h>
#include <CGAL/Cartesian.h>
#include <CGAL/Polygon_2.h>

#include <CGAL/Arr_circle_segment_traits_2.h>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Arr_linear_traits_2.h>

#include <CGAL/Arr_geometry_traits/Circle_segment_2.h>

#include <CGAL/Arrangement_2.h>
#include <CGAL/Arr_extended_dcel.h>
#include <CGAL/Arr_overlay_2.h>
#include <CGAL/Arr_default_overlay_traits.h>
#include <hedra/copyleft/cgal/basic_cgal_definitions.h>
#include <hedra/dcel.h>
#include <CGAL/to_rational.h>
#include <CGAL/Exact_integer.h>


namespace hedra
{
  namespace copyleft
  {
    namespace cgal
    {

      const int PARAM_LINE_VERTEX = -2;
      const int ORIGINAL_VERTEX = -1;

      struct ArrEdgeData
          {
        bool isParam;
        int origEdge;
        int newHalfedge;
      };

      //f1: original mesh
      //f2: parameter lines
      template<class ArrangementA, class ArrangementB, class ArrangementR>
      class Arr_mesh_generation_overlay_traits :
          public CGAL::_Arr_default_overlay_traits_base<ArrangementA, ArrangementB, ArrangementR> {
      public:

        typedef typename ArrangementA::Face_const_handle Face_handle_A;
        typedef typename ArrangementB::Face_const_handle Face_handle_B;
        typedef typename ArrangementR::Face_handle Face_handle_C;

        typedef typename ArrangementA::Vertex_const_handle Vertex_handle_A;
        typedef typename ArrangementB::Vertex_const_handle Vertex_handle_B;
        typedef typename ArrangementR::Vertex_handle Vertex_handle_C;

        typedef typename ArrangementA::Halfedge_const_handle Halfedge_handle_A;
        typedef typename ArrangementB::Halfedge_const_handle Halfedge_handle_B;
        typedef typename ArrangementR::Halfedge_handle Halfedge_handle_C;

      public:

        virtual void create_face(Face_handle_A f1, Face_handle_B f2, Face_handle_C f) const
        {
          // Overlay the data objects associated with f1 and f2 and store the result with f.
          f->set_data(f1->data());
        }

        virtual void create_vertex(Vertex_handle_A v1, Vertex_handle_B v2, Vertex_handle_C v)
        {
          v->set_data(PARAM_LINE_VERTEX);
        }

        virtual void create_vertex(Vertex_handle_A v1, Halfedge_handle_B e2, Vertex_handle_C v)
        {
          v->set_data(ORIGINAL_VERTEX);
        }

        virtual void create_vertex(Vertex_handle_A v1, Face_handle_B f2, Vertex_handle_C v)
        {
          v->set_data(ORIGINAL_VERTEX);
        }

        virtual void create_vertex(Halfedge_handle_A e1, Vertex_handle_B v2, Vertex_handle_C v)
        {
          v->set_data(PARAM_LINE_VERTEX);
        }

        virtual void create_vertex(Face_handle_A f1, Vertex_handle_B v2, Vertex_handle_C v)
        {
          v->set_data(PARAM_LINE_VERTEX);
        }

        virtual void create_vertex(Halfedge_handle_A e1, Halfedge_handle_B e2, Vertex_handle_C v)
        {
          v->set_data(ORIGINAL_VERTEX);
        }

        virtual void create_edge(Halfedge_handle_A e1, Halfedge_handle_B e2, Halfedge_handle_C e)
        {
          ArrEdgeData aed;
          aed.isParam = true;
          aed.origEdge = e1->data().origEdge;
          aed.newHalfedge = -1;
          e->set_data(aed);
          e->twin()->set_data(aed);
        }

        virtual void create_edge(Halfedge_handle_A e1, Face_handle_B f2, Halfedge_handle_C e)
        {
          ArrEdgeData aed;
          aed.isParam = false;
          aed.origEdge = e1->data().origEdge;
          aed.newHalfedge = -1;
          e->set_data(aed);
          e->twin()->set_data(aed);
        }

        virtual void create_edge(Face_handle_A f1, Halfedge_handle_B e2, Halfedge_handle_C e)
        {
          ArrEdgeData aed;
          aed.isParam = true;
          aed.origEdge = -1;
          aed.newHalfedge = -1;
          e->set_data(aed);
          e->twin()->set_data(aed);
        }
      };

      typedef CGAL::Arr_linear_traits_2<Kernel>           Traits2;
      typedef Traits2::Point_2                             Point2;
      typedef Traits2::Segment_2                           Segment2;
      typedef Traits2::Line_2                              Line2;
      typedef Traits2::X_monotone_curve_2                  X_monotone_curve_2;

      typedef CGAL::Arr_extended_dcel<Traits2, int,ArrEdgeData,int>  Dcel;
      typedef CGAL::Arrangement_2<Traits2, Dcel>                  Arr_2;
      typedef Arr_2::Face_iterator                          Face_iterator;
      typedef Arr_2::Face_handle                            Face_handle;
      typedef Arr_2::Edge_iterator                          Edge_iterator;
      typedef Arr_2::Halfedge_iterator                      Halfedge_iterator;
      typedef Arr_2::Vertex_iterator                        Vertex_iterator;
      typedef Arr_2::Vertex_handle                          Vertex_handle;
      typedef Arr_2::Halfedge_handle                        Halfedge_handle;
      typedef Arr_2::Ccb_halfedge_circulator                Ccb_halfedge_circulator;
      typedef Arr_mesh_generation_overlay_traits <Arr_2, Arr_2,Arr_2>  Overlay_traits;


      //for now doing quad (u,v,-u -v) only!
      Point2 paramCoord2texCoord(const Eigen::RowVectorXd & paramCoord)
      {
        return Point2(Number(paramCoord(0)), Number(paramCoord(1)));
      }


      // Input:

      // UV           #V by 2, vertex coordinates in the parameter space
      // FUV          #F by 3, map between faces' vertices and coordinates in the parameter space
      // resolution   # rational number accuracy
      // ti           # index of the current tirangle

      // Output:

      // paramArr     # square grid pattern in the CGAL representation
//      void square_grid_pattern(const Eigen::MatrixXd & UV, const Eigen::MatrixXi & FUV, const int resolution, const int ti, Arr_2 & paramArr)
//      {
//        //creating an arrangement of parameter lines
//        Eigen::MatrixXd facePC(3, UV.cols()); // PC.cols == 2
//        for (int i = 0; i < 3; i++)
//          facePC.row(i) = UV.row(FUV(ti, i));
//
//        for (int i = 0; i < facePC.cols(); i++)
//        {
//          //inserting unbounded lines
//          int coordMin = (int) std::floor(facePC.col(i).minCoeff() - 1.0);
//          int coordMax = (int) std::ceil(facePC.col(i).maxCoeff() + 1.0);
//          std::vector<Line2> lineCurves;
//
//          for (int coordIndex = coordMin; coordIndex <= coordMax; coordIndex++)
//          {
//            //The line coord = coordIndex
//            Eigen::RowVectorXd LineCoord1 = Eigen::RowVectorXd::Zero(facePC.cols());
//            Eigen::RowVectorXd LineCoord2 = Eigen::RowVectorXd::Ones(facePC.cols());
//            LineCoord1(i) = coordIndex;
//            LineCoord2(i) = coordIndex;
//            lineCurves.emplace_back(paramCoord2texCoord(LineCoord1, resolution), paramCoord2texCoord(LineCoord2, resolution));
//          }
//          insert(paramArr, lineCurves.begin(), lineCurves.end());
//        }
//      }

      // the signature is the same as for the suare grid pattern function
      void tri_grid_pattern(const Eigen::MatrixXd & UV, const Eigen::MatrixXi & FUV, const int resolution, const int ti, Arr_2 & paramArr)
      {
        //creating an arrangement of parameter lines
        // try to round on the grid

        Eigen::Matrix2d cH;
        cH << std::sqrt(3.), -std::sqrt(3.) / 2., 0., -3. / 2.;
        std::vector<Number> coordsX(3);
        std::vector<Number> coordsY(3);

        auto sqrt_3 = CGAL::sqrt(Number(3));

        Eigen::MatrixXd facePC(3, UV.cols()); // PC.cols == 2
        for (int i = 0; i < 3; i++) {
          facePC.row(i) = cH.inverse() * UV.row(FUV(ti, i)).transpose(); // back to the axial coordinates
          // round in the cube coordinates

            Eigen::Vector3d cube(facePC(i, 0), -facePC(i, 1), -facePC(i, 0) + facePC(i, 1));
            Eigen::Vector3d cubeR(std::round(facePC(i, 0)), std::round(facePC(i, 1)), std::round(-facePC(i, 0) + facePC(i, 1)));
            Eigen::Vector3d diff(std::fabs(cubeR(0) - cube(0)), std::fabs(cubeR(1) - cube(1)), std::fabs(cubeR(2) - cube(2)));

            if(diff(0) > diff(1) && diff(0) > diff(2))
            {
              facePC(i, 0) = (cubeR(1) - cubeR(2));
              facePC(i, 1) = cubeR(1);
            }
            else if (diff(1) > diff(2))
            {
              facePC(i, 0) = cubeR(0);
              facePC(i, 1) = (cubeR(0) + cubeR(2));
            }
            else
            {
              facePC(i, 0) = cubeR(0);
              facePC(i, 1) = cubeR(1);
            }
            Eigen::Vector2d p = facePC.row(i);
            coordsX[i] = Number((int)p(0)) * sqrt_3 - Number((int)p(1)) * sqrt_3 / Number(2);
            coordsY[i] = Number((int)p(1)) * Number(-3) / Number(2);
        }

        // find min and max x
        Number coordMinY = *(std::min_element(coordsY.cbegin(), coordsY.cend()));
        Number coordMaxY = *(std::max_element(coordsY.cbegin(), coordsY.cend()));
        Number coordMinX = *(std::min_element(coordsX.cbegin(), coordsX.cend()));
        Number coordMaxX = *(std::max_element(coordsX.cbegin(), coordsX.cend()));

        //inserting unbounded lines -- vertical
          std::vector<Line2> lineCurves;
        Number incX = sqrt_3 / Number(4);
        Number incY = Number(3) / Number(2);

        Number yShift = Number(1) / Number(2);
        Number xShift = sqrt_3 / Number(2);

        int c = 0;

          for (Number coordIndexX = coordMinX - sqrt_3; coordIndexX <= coordMaxX + sqrt_3; coordIndexX += incX) {
            lineCurves.emplace_back(Point2(coordIndexX, Number(0)), Point2(coordIndexX, Number(1)));
            Number coordIndexY  = coordMinY - Number(6);
            if (c % 2 != 0)
              coordIndexY = coordMinY - Number(3) / Number(4);
            for (; coordIndexY <= coordMaxY + Number(6); coordIndexY += incY) {
              lineCurves.emplace_back(Point2(coordIndexX, coordIndexY), Point2(coordIndexX + xShift, coordIndexY + yShift));
              lineCurves.emplace_back(Point2(coordIndexX, coordIndexY), Point2(coordIndexX + xShift, coordIndexY - yShift));
            }

            c++;
          }
        insert(paramArr, lineCurves.begin(), lineCurves.end());
      }


    struct PointPair{
        int Index1, Index2;
        double Distance;

        PointPair(int i1, int i2, double d):Index1(i1), Index2(i2), Distance(d){}
        ~PointPair(){}

        const bool operator<(const PointPair& pp) const {
          if (Distance>pp.Distance) return false;
          if (Distance<pp.Distance) return true;

          if (Index1>pp.Index1) return false;
          if (Index1<pp.Index1) return true;

          if (Index2>pp.Index2) return false;
          if (Index2<pp.Index2) return true;

          return false;

        }
      };

     void vertex_sets_match(const vector<Point3D>& set1, const vector<Point3D>& set2, std::vector<std::pair<int,int> > & result)
      {
        std::set<PointPair> pairSet;
        for (int i = 0; i < set1.size(); i++)
          for (int j = 0; j < set2.size(); j++)
            pairSet.insert(PointPair(i, j, Norm(set1[i] - set2[j])));

        //adding greedily legal connections until graph is full
        vector<bool> set1Connect(set1.size(), false);
        vector<bool> set2Connect(set2.size(), false);

        int numConnected = 0;
        for (auto ppi = pairSet.begin(); ppi != pairSet.end(); ppi++)
        {
          PointPair currPair = *ppi;
          //checking legality - if any of one's former are connected to ones latters or vice versa
          bool foundConflict = false;
          for (int i = 0; i < result.size(); i++) {
            if (((result[i].first > currPair.Index1) && (result[i].second < currPair.Index2)) ||
                ((result[i].first < currPair.Index1) && (result[i].second > currPair.Index2))) {
              foundConflict = true;
              break;
            }
          }
          if (foundConflict) {
            continue;
          }

          //std::cout << "Match distance " << currPair.Distance << std::endl;

          //otherwise this edge is legal, so add it
          result.push_back(pair<int, int>(currPair.Index1, currPair.Index2));
          if (!set1Connect[currPair.Index1]) numConnected++;
          if (!set2Connect[currPair.Index2]) numConnected++;
          set1Connect[currPair.Index1] = set2Connect[currPair.Index2] = true;
          if (numConnected == set1.size() + set2.size())
            break;  //all nodes are connected
        }
      }

    void find_maching_vertices(const std::vector<std::vector<int> > & boundEdgesLeft,
                               const std::vector<std::vector<int> > & boundEdgesRight,
                               const Eigen::VectorXi & nextH,
                               const std::vector<Point3D> & HE3D,
                               const Eigen::VectorXi & HV,
                               std::vector<std::pair<int, int> > & vertexMatches)
    {
      std::vector<std::vector<int> > vertexSetL(boundEdgesLeft.size()), vertexSetR(boundEdgesRight.size());
        for (size_t i = 0; i < boundEdgesLeft.size(); i++) {
          for (auto eid : boundEdgesLeft[i])
            vertexSetL[i].push_back(HV(eid));
          if (!boundEdgesLeft[i].empty())
            vertexSetL[i].push_back(HV(nextH(boundEdgesLeft[i].back())));
        }

        for (size_t i = 0; i < boundEdgesRight.size(); i++) {
          for (auto eid : boundEdgesRight[i])
            vertexSetR[i].push_back(HV(eid));
          if (!boundEdgesRight[i].empty())
            vertexSetR[i].push_back(HV(nextH(boundEdgesRight[i].back())));
          std::reverse(vertexSetR[i].begin(), vertexSetR[i].end());
        }


        for (size_t i = 0; i < boundEdgesRight.size(); i++) {
          std::vector<Point3D> pointSetL(vertexSetL[i].size());
          std::vector<Point3D> pointSetR(vertexSetR[i].size());

          for (size_t j = 0; j < pointSetL.size(); j++)
            pointSetL[j] = HE3D[vertexSetL[i][j]];
          for (size_t j = 0; j < pointSetR.size(); j++)
            pointSetR[j] = HE3D[vertexSetR[i][j]];

          std::vector<std::pair<int, int> > currMatches;
          if ((!pointSetL.empty()) && (!pointSetR.empty()))
            vertex_sets_match(pointSetL, pointSetR, currMatches);

          for (size_t j = 0; j < currMatches.size(); j++) {
            currMatches[j].first = vertexSetL[i][currMatches[j].first];
            currMatches[j].second = vertexSetR[i][currMatches[j].second];
          }
          vertexMatches.insert(vertexMatches.end(), currMatches.begin(), currMatches.end());
        }
    }

      struct TwinFinder{
        int index;
        int v1, v2;

        TwinFinder(int i, int vv1, int vv2): index(i), v1(vv1), v2(vv2){}
        ~TwinFinder(){}

        const bool operator<(const TwinFinder& tf) const
        {
          if (v1<tf.v1) return false;
          if (v1>tf.v1) return true;

          if (v2<tf.v2) return false;
          if (v2>tf.v2) return true;

          return false;
        }
      };


      void removeFace(const int findex,
                      const int heindex,
                      Eigen::VectorXi & HV,
                      Eigen::VectorXi & VH,
                      Eigen::VectorXi & HF,
                      Eigen::VectorXi & FH,
                      Eigen::VectorXi & twinH,
                      Eigen::VectorXi & nextH,
                      Eigen::VectorXi & prevH,
                      std::vector<bool> & validHE,
                      std::vector<bool> & validV,
                      std::vector<bool> & validF,
                      const std::vector<int> HE2origEdges,
                      std::ostream & log = std::cout)
      {
        int leftVertex = HV(heindex);

        int hebegin = heindex;
        int heiterate = hebegin;

        log << "Removing Face " << findex << " with initial edge "<< heindex << "\n";
        log << "Leaving Vertex " << leftVertex<<"\n";

        validF[findex] = false;
        std::vector<int> replaceOrigins;
        do{
          log << "Removing edge part " << heiterate << " with twin " << twinH(heiterate) <<" and origin " << HV(heiterate) <<"\n";
          validHE[heiterate] = false;
          if (HF(heiterate) != findex)
            throw std::runtime_error("Wrong face!");
          if (HV(heiterate) != leftVertex) {
            validV[HV(heiterate)] = false;
            replaceOrigins.push_back(HV(heiterate));
          }
          if (twinH[heiterate] == -1 ) {
            heiterate = nextH(heiterate);
            continue;
          }
          int reduceEdge = twinH(heiterate);
          validHE[reduceEdge] = false;
          prevH(nextH(reduceEdge)) = prevH(reduceEdge);
          nextH(prevH(reduceEdge)) = nextH(reduceEdge);
          log << "Connecting " << prevH(reduceEdge) << " -> " << nextH(reduceEdge) << "\n";

          FH(HF(reduceEdge)) = nextH(reduceEdge);
          VH(leftVertex) = nextH(reduceEdge);
          log << "Vertex " << leftVertex << " now points to " << nextH(reduceEdge) <<"\n";

          heiterate = nextH(heiterate);
        } while (heiterate != hebegin);

        for (int i = 0; i < HV.rows(); i++)
          for (size_t j = 0; j < replaceOrigins.size(); j++)
            if (HV(i) == replaceOrigins[j]) {
              log << "Now halfedge " << i << " has vertex " << leftVertex << " as origin instead of " << replaceOrigins[j] <<"\n";
              HV(i) = leftVertex;
            }

        if (!validHE[VH(leftVertex)]) {
          log << "Discarting Vertex " << leftVertex << "\n";
          validV[leftVertex] = false;
        }
      }

      void removeEdge(const int heindex,
                      Eigen::VectorXi & HV,
                      Eigen::VectorXi & VH,
                      Eigen::VectorXi & HF,
                      Eigen::VectorXi & FH,
                      Eigen::VectorXi & twinH,
                      Eigen::VectorXi & nextH,
                      Eigen::VectorXi & prevH,
                      std::vector<bool> & validHE,
                      std::vector<bool> & validV,
                      std::vector<bool> & validF,
                      const std::vector<int> HE2origEdges,
                      std::ostream & log = std::cout
      )
      {
        if (twinH(heindex) != -1) {
          int ebegin = twinH(heindex);
          int ecurr = ebegin;
          int counter = 0;
          do {
            counter++;
            ecurr = nextH(ecurr);
          } while (ebegin != ecurr);
          if(counter <= 3) {
            removeFace(HF(twinH(heindex)), twinH(heindex), HV, VH, HF, FH, twinH, nextH, prevH, validHE, validV, validF, HE2origEdges);
            return;
          }
        }
        log << "Removing Edge " << heindex << " with Twin " << twinH(heindex) << std::endl;
        validHE[heindex] = false;
        prevH(nextH(heindex)) = prevH(heindex);
        nextH(prevH(heindex)) = nextH(heindex);
        log << "Connecting " << prevH(heindex) << " -> " << nextH(heindex) << std::endl;

        VH(HV(heindex)) = nextH(heindex);
        log << "Vertex " << HV(heindex) << " points to " << nextH(heindex) <<"\n";
        int leftVertex = HV(heindex);
        int removeVertex = HV(nextH(heindex));
        log << "Vertex " << removeVertex << " has been removed" << "\n";
        validV[removeVertex] = false;
        HV(nextH(heindex)) = HV(heindex);
        log << "half-edge " << nextH(heindex) << " has the vertex " << HV(heindex) << " as its origin \n";
        FH(HF(heindex)) = nextH(heindex);

        if (twinH(heindex) != -1) {
          validHE[twinH(heindex)] = false;
          prevH(nextH(twinH(heindex))) = prevH(twinH(heindex));
          nextH(prevH(twinH(heindex))) = nextH(twinH(heindex));
          log << "Connecting " << prevH(twinH(heindex)) << "->" << nextH(twinH(heindex)) <<"\n";
          HV(nextH(twinH(heindex))) = leftVertex;
          FH(HF(twinH(heindex))) = nextH(twinH(heindex));
        }

        for (int i = 0; i < HV.rows(); i++) {
          if (HV(i) == removeVertex)
            HV(i) = leftVertex;
          log << "now halfedge " << i << " has the vertex " << leftVertex << " as its origin \n";
        }
      }

      bool checkMesh(const Eigen::VectorXi & HV,
                     const Eigen::VectorXi & VH,
                     const Eigen::VectorXi & HF,
                     const Eigen::VectorXi & FH,
                     const Eigen::VectorXi & twinH,
                     const Eigen::VectorXi & nextH,
                     const Eigen::VectorXi & prevH,
                     const std::vector<bool> & isParamHE,
                     const std::vector<bool> & validHE,
                     const std::vector<bool> & validV,
                     const std::vector<bool> & validF)
    {
      for (size_t i = 0; i < validV.size(); i++) {
      if (!validV[i])
      continue;

      if (!validHE[VH(i)])
        throw std::runtime_error("Valid vertex adjacent to an invalid edge!");

      if (HV(VH(i)) != i)
        throw std::runtime_error("Edge adjacent to the vertex has no the vertex as its source!");
      }

      for (size_t i = 0; i < validHE.size(); i++) {
        if (!validHE[i])
          continue;

        if (prevH(nextH(i)) != i)
          throw std::runtime_error("PrevHE of nextHE is not correct!");

        if (nextH(prevH(i)) != i)
          throw std::runtime_error("NextHE of prevHE is not correct!");

        if (!validV[HV(i)])
          throw std::runtime_error("The source vertex is not correct!");

        if (!validF[HF(i)])
          throw std::runtime_error("The adjacent face is invalid!");

        if (twinH(i) != -1)
          if (twinH(twinH(i)) != i)
            throw std::runtime_error("TwinHE of twinHE is not correct!");

        if (isParamHE[i]) {  //checking that it is not left alone
            if (prevH(i) == twinH(i))
              throw std::runtime_error("PrevHE is the twinHE!");

          if (nextH(i) == twinH(i))
            throw std::runtime_error("NextHE is the twinHE!");
        }
        return true;
      }

        for (size_t i = 0; i < validF.size(); i++) {
          if (!validF[i])
            continue;

          int hebegin = FH(i);
          int heiterate = hebegin;
          int NumEdges = 0;

          do{
            if (!validHE[heiterate])
              throw std::runtime_error("Invalid HE adjacent to the face!");

            if (HF(heiterate) != i)
              throw std::runtime_error("Adjacent face is not correct!");

            heiterate= nextH(heiterate);
            NumEdges++;
            if (NumEdges > 1000)
              throw std::runtime_error("The face has no correct curcuite!");
          }while (heiterate!=hebegin);
        }
      }

     void graph_verification(const std::vector<std::vector<int> > & boundEdgesLeft,
                             const std::vector<std::vector<int> > & boundEdgesRight,
                             Eigen::MatrixXd & currV,
                             Eigen::VectorXi & nextH,
                             Eigen::VectorXi & prevH,
                             Eigen::VectorXi & twinH,
                             Eigen::VectorXi & HV,
                             Eigen::VectorXi & VH,
                             Eigen::VectorXi & HF,
                             Eigen::VectorXi & FH,
                             std::vector<std::pair<int, int> > & vertexMatches,
                             std::vector<Point3D> & HE3D,
                             std::vector<int> & transVertices,
                             std::vector<bool> & isParamVertex,
                             std::vector<bool> & validHE, std::vector<bool> & validV, std::vector<bool> & validF,
                             const std::vector<bool> & isParamHE,
                             const std::vector<int> & HE2origEdges,
                             std::ostream & log = std::cout)
     {
       typedef boost::adjacency_list <boost::vecS, boost::vecS, boost::undirectedS> Graph;
       Graph MatchGraph;

       for (size_t i = 0; i < HE3D.size(); i++)
         boost::add_vertex(MatchGraph);
       for (size_t i = 0; i < vertexMatches.size(); i++)
         boost::add_edge(vertexMatches[i].first, vertexMatches[i].second, MatchGraph);

       int numNewVertices = boost::connected_components(MatchGraph, &transVertices[0]);

       //adding other vertex to the degeneration if needed
       bool thereisChange;
       do{
         thereisChange = false;
         for (int i = 0; i < twinH.size(); i++) {
           if ((!validHE[i]) || (twinH(i) != -1))
             continue;
           if (transVertices[HV(i)] != transVertices[HV(nextH(i))])
             continue;  //this edge is OK

           int ebegin = i;
           int ecurr = ebegin;
           std::vector<int> FV;
           do {
             FV.push_back(HV(ecurr));
             ecurr = nextH(ecurr);
           } while (ebegin != ecurr);

           if (FV.size() <= 3) {
             for (size_t j = 0; j < FV.size(); j++) {
               if (transVertices[FV[j]] != transVertices[HV(i)]){
                 boost::add_edge(FV[j], HV(i), MatchGraph);
                 thereisChange = true;
               }
             }
           }
         }

         if (thereisChange){
           transVertices.clear();
           transVertices.resize(HE3D.size());
           numNewVertices = boost::connected_components(MatchGraph, &transVertices[0]);
         }
       }while (thereisChange);

       //removing edges (and consequent faces) which will degenerate
       for (int i = 0; i < twinH.rows(); i++) {
         if (!validHE[i] || twinH(i) != -1)
           continue;
         if (transVertices[HV(i)] != transVertices[HV(nextH(i))])
           continue;  //this edge is OK

         int ebegin = i;
         int ecurr = ebegin;
         int counter = 0;
         do {
           counter++;
           ecurr = nextH(ecurr);
         } while (ebegin != ecurr);

         if (counter <= 3) {
           std::cout << "check before remove face!" << std::endl;
           checkMesh(HV, VH, HF, FH, twinH, nextH, prevH, isParamHE, validHE, validV, validF);
           removeFace(HF(i), i, HV, VH, HF, FH, twinH, nextH, prevH, validHE, validV, validF, HE2origEdges, log);
           std::cout << "check after remove face!" << std::endl;
           checkMesh(HV, VH, HF, FH, twinH, nextH, prevH, isParamHE, validHE, validV, validF);
         }
         else {
           std::cout << "check before remove edge!" << std::endl;
           checkMesh(HV, VH, HF, FH, twinH, nextH, prevH, isParamHE, validHE, validV, validF);
           removeEdge(i, HV, VH, HF, FH, twinH, nextH, prevH, validHE, validV, validF, HE2origEdges, log);
           std::cout << "check after remove edge!" << std::endl;
           checkMesh(HV, VH, HF, FH, twinH, nextH, prevH, isParamHE, validHE, validV, validF);
         }
       }

       std::vector<Point3D> newVertices(numNewVertices);
       std::vector<bool> newValidV(numNewVertices);
       std::vector<bool> newIsParamV(numNewVertices);
       currV.resize(numNewVertices, 3);
       Eigen::VectorXi newVH(numNewVertices);
       for (size_t i = 0; i < HE3D.size(); i++) {
         if(!validV[i])
           continue;
         Point3D newVertex = HE3D[i];
         currV.row(transVertices[i]) = Eigen::RowVector3d(CGAL::to_double(newVertex.x()),
                                                          CGAL::to_double(newVertex.y()),
                                                          CGAL::to_double(newVertex.z()));
         newVertices[transVertices[i]] = newVertex;
         newValidV[transVertices[i]] = validV[i];
         newVH(transVertices[i]) = VH(i);
         newIsParamV[transVertices[i]] = isParamVertex[i];
       }
       HE3D = newVertices;
       VH = newVH;
       isParamVertex = newIsParamV;
       validV = newValidV;

       for (int i = 0; i < HV.rows(); i++) {
         HV(i) = transVertices[HV(i)];
         if (validHE[i])
           VH(HV(i)) = i;
       }
       std::cout << "check after vertex adjust!" << std::endl;
       checkMesh(HV, VH, HF, FH, twinH, nextH, prevH, isParamHE, validHE, validV, validF);
     }

     void joinFace(const int heindex,
      Eigen::VectorXi & twinH,
      Eigen::VectorXi & prevH,
      Eigen::VectorXi & nextH,
      Eigen::VectorXi & HF,
      Eigen::VectorXi & FH,
      Eigen::VectorXi & HV,
      Eigen::VectorXi & VH,
      std::vector<bool> & validHE, std::vector<bool> & validV, std::vector<bool> & validF)
    {
      if (twinH(heindex) == -1)
          return;

        int fid0 = HF(heindex);
        int fid1 = HF(twinH(heindex));

        //check if a lonely edge
        if ((prevH(heindex) == twinH(heindex)) && (nextH(heindex) == twinH(heindex))) {
          validHE[heindex] = validHE[twinH(heindex)] = false;
          validV[HV(heindex)] = validV[HV(twinH(heindex))] = false;
          if ((FH(HF(heindex)) == heindex) || (twinH((FH(HF(heindex)))) == heindex)) {
            for (size_t i = 0; i < validHE.size(); i++) {
              if (!validHE[i])
                continue;
              if (HF(i) == HF(heindex)) {
                FH(HF(heindex)) = i;
                break;
              }
            }
          }
          return;
        }

        //check if spike edge
        if ((prevH(heindex) == twinH(heindex)) || (nextH(heindex) == twinH(heindex))) {
          int closeEdge = heindex;
          if (prevH(heindex) == twinH(heindex))
            closeEdge = twinH(heindex);

          validHE[closeEdge] = validHE[twinH(closeEdge)] = false;
          VH(HV(closeEdge)) = nextH(twinH(closeEdge));
          FH(fid0) = prevH(closeEdge);

          nextH(prevH(closeEdge)) = nextH(twinH(closeEdge));
          prevH(nextH(twinH(closeEdge))) = prevH(closeEdge);
          validV[HV(twinH(closeEdge))] = false;
          return;
        }

        FH(fid0) = nextH(heindex);
        if (fid0 != fid1)
          validF[fid1] = false;

        validHE[heindex] = validHE[twinH(heindex)] = false;

        prevH(nextH(heindex)) = prevH(twinH(heindex));
        nextH(prevH(twinH(heindex))) = nextH(heindex);
        prevH(nextH(twinH(heindex))) = prevH(heindex);
        nextH(prevH(heindex)) = nextH(twinH(heindex));

        VH(HV(heindex)) = nextH(twinH(heindex));
        VH(HV(nextH(heindex))) = nextH(heindex);

        //all other floating halfedges should renounce this one
        for (int i = 0; i < HF.rows(); i++)
          if (HF(i) == fid1)
            HF(i) = fid0;
      }

//edge has to be sourcing a 2-valence vertex!!
      void unifyEdges(const int heindex,
                      Eigen::VectorXi & twinH,
                      Eigen::VectorXi & prevH,
                      Eigen::VectorXi & nextH,
                      Eigen::VectorXi & HF,
                      Eigen::VectorXi & FH,
                      Eigen::VectorXi & HV,
                      Eigen::VectorXi & VH,
                      std::vector<bool> & validHE, std::vector<bool> & validV)
      {
        //adjusting source
        validV[HV(heindex)] = false;
        HV(heindex) = HV(prevH(heindex));
        VH(HV(heindex)) = heindex;
        FH(HF(heindex)) = nextH(heindex);

        //adjusting halfedges
        validHE[prevH(heindex)] = false;
        prevH(heindex) = prevH(prevH(heindex));
        nextH(prevH(heindex)) = heindex;

        //adjusting twin, if exists
        if (twinH(heindex) != -1) {
          validHE[nextH(twinH(heindex))] = false;
          nextH(twinH(heindex)) = nextH(nextH(twinH(heindex)));
          prevH(nextH(twinH(heindex))) = twinH(heindex);
          FH(HF(twinH(heindex))) = nextH(twinH(heindex));
        }
      }

      void testUnmatchedTwins(const Eigen::VectorXi & twinH, const Eigen::VectorXi & nextH, const Eigen::VectorXi & HV, const std::vector<bool> & validHE, const std::vector<Point3D> & HE3D, const std::vector<int> & HE2origEdges, std::ostream & log = std::cout)
    {
      vector<int> untwinned;
        for (int i=0; i < validHE.size(); i++)
          if ((twinH(i) == -1) && (validHE[i]))
            untwinned.push_back(i);

        for (size_t i = 0; i < untwinned.size(); i++) {
          for (size_t j = i + 1; j < untwinned.size(); j++) {
            Vector3D diff1 = HE3D[HV(untwinned[i])] - HE3D[HV(nextH(untwinned[j]))];
            Vector3D diff2 = HE3D[HV(untwinned[j])] - HE3D[HV(nextH(untwinned[i]))];
            if ((CGAL::to_double(CGAL::sqrt(diff1.squared_length())) < 10e-4) &&(CGAL::to_double(CGAL::sqrt(diff2.squared_length())) < 10e-4)) {
              log << "Halfedge " << untwinned[i] << " (org: " << HE2origEdges[untwinned[i]] << "), :(" << HV(untwinned[i]) << ", " << HV(nextH(untwinned[i])) << ") is untwinned to ";
              log << "Halfedge " << untwinned[j]  << " (org: " << HE2origEdges[untwinned[j]] << "), :(" << HV(untwinned[j]) <<", " << HV(nextH(untwinned[j])) << ")\n";
              log << HE3D[HV(untwinned[i])] << " and " << HE3D[HV(nextH(untwinned[i]))] << std::endl;
              log << HE3D[HV(untwinned[j])] << " and "<< HE3D[HV(nextH(untwinned[j]))] << std::endl;
            }
          }
        }
      }


      void cleanMesh(const std::vector<bool> & validHE,
                     const std::vector<bool> & validV,
                     const std::vector<bool> & validF,
                     Eigen::VectorXi & twinH,
                     Eigen::VectorXi & prevH,
                     Eigen::VectorXi & nextH,
                     Eigen::MatrixXd & currV,
                     Eigen::VectorXi & HF,
                     Eigen::VectorXi & FH,
                     Eigen::VectorXi & HV,
                     Eigen::VectorXi & VH,
                     std::vector<bool> & isParamVertex,
                     std::vector<int> & HE2origEdges,
                     std::vector<bool> & isParamHE)
      {
        //edges
        for (int hid = validHE.size() - 1; hid >= 0; hid--) {
          // skipe a rediscovered useless edge, well normally it should not happen
          if (validHE[hid])
            continue;
          int numRows = HE2origEdges.size() - 1;

          if (hid < numRows) {
            HF.segment(hid, numRows - hid) = HF.segment(hid + 1, numRows - hid).eval();
            HV.segment(hid, numRows - hid) = HV.segment(hid + 1, numRows - hid).eval();
            nextH.segment(hid, numRows - hid) = nextH.segment(hid + 1, numRows - hid).eval();
            twinH.segment(hid, numRows - hid) = twinH.segment(hid + 1, numRows - hid).eval();
            prevH.segment(hid, numRows - hid) = prevH.segment(hid + 1, numRows - hid).eval();
          }

          HF.conservativeResize(numRows);
          HV.conservativeResize(numRows);
          nextH.conservativeResize(numRows);
          twinH.conservativeResize(numRows);
          prevH.conservativeResize(numRows);

          HE2origEdges.erase(HE2origEdges.cbegin() + hid);
          isParamHE.erase(isParamHE.cbegin() + hid);

          //update IDs
          for (int k = 0; k < FH.rows(); k++)
            if (FH(k) > hid)
              FH(k)--;

          for (int k = 0; k < VH.rows(); k++)
            if (VH(k) > hid)
              VH(k)--;

          for (int k = 0; k < nextH.rows(); k++)
            if (nextH(k) > hid)
              nextH(k)--;

          for (int k = 0; k < prevH.rows(); k++)
            if (prevH(k) > hid)
              prevH(k)--;

          for (int k = 0; k < twinH.rows(); k++)
            if (twinH(k) > hid)
              twinH(k)--;
      }

        //removed unreferenced faces
        for (int fid = validF.size() - 1; fid >= 0; fid--) {
          if (validF[fid])
            continue;
          //remove the row
          int numRows = FH.rows() - 1;
          if (fid < numRows)
            FH.segment(fid, numRows - fid) = FH.segment(fid + 1, numRows - fid).eval();
          FH.conservativeResize(numRows);
          //update IDs
          for (int k = 0; k < HF.rows(); k++)
            if (HF(k) > fid)
              HF(k)--;
        }

        //removed unreferenced faces
        std::vector<int> hitFaces(FH.rows(), 0);
        for (int k = 0; k < HF.rows(); k++) {
          hitFaces[HF(k)]++;
        }

        for (int fid = hitFaces.size() - 1; fid >= 0; fid--) {
          if (hitFaces[fid])
            continue;
          //remove the row
          int numRows = FH.rows() - 1;
          if (fid < numRows)
            FH.segment(fid, numRows - fid) = FH.segment(fid + 1, numRows - fid).eval();
          FH.conservativeResize(numRows);
          //update IDs
          for (int k = 0; k < HF.rows(); k++)
            if (HF(k) > fid)
              HF(k)--;
        }

        //removed non-valid vertices
        for (int vid = validV.size() - 1; vid >= 0; vid--) {
          if (validV[vid])
            continue;
          //remove the row
          int numRows = VH.rows() - 1;
          if (vid < numRows)
          {
            currV.block(vid, 0, numRows - vid, 3) = currV.block(vid + 1, 0, numRows - vid, 3).eval();
            VH.segment(vid, numRows - vid) = VH.segment(vid + 1, numRows - vid).eval();
          }
          currV.conservativeResize(numRows, 3);
          VH.conservativeResize(numRows);
          isParamVertex.erase(isParamVertex.begin() + vid);
          //update IDs
          for (int k = 0; k < HV.rows(); k++)
            if (HV(k) > vid)
              HV(k)--;
        }

        //removed unreferenced vertices
        std::vector<int> hitVers(VH.rows(), 0);
        for (int k = 0; k < HV.rows(); k++) {
          hitVers[HV(k)]++;
        }

        for (int vid = hitVers.size() - 1; vid >= 0; vid--) {
          if (hitVers[vid] > 0)
            continue;
          //remove the row
          int numRows = VH.rows() - 1;
          if (vid < numRows)
          {
            currV.block(vid, 0, numRows - vid, 3) = currV.block(vid + 1, 0, numRows - vid, 3).eval();
            VH.segment(vid, numRows - vid) = VH.segment(vid + 1, numRows - vid).eval();
          }
          currV.conservativeResize(numRows, 3);
          VH.conservativeResize(numRows);
          isParamVertex.erase(isParamVertex.begin() + vid);
          //update IDs
          for (int k = 0; k < HV.rows(); k++)
            if (HV(k) > vid)
              HV(k)--;
        }
      }

      void WalkBoundary(int &CurrEdge,
                        Eigen::VectorXi & nextH,
                        Eigen::VectorXi & twinH)
    {
        do{
          CurrEdge = nextH(CurrEdge);
          if (twinH(CurrEdge) == -1)
            break;  //next boundary over a 2-valence vertex
          CurrEdge = twinH(CurrEdge);
        }while (twinH(CurrEdge) != -1);
      }


      IGL_INLINE void stitch_boundaries2(
          std::vector<Point3D> & HE3D,
          int resolution,
          const Eigen::MatrixXd & V,
          const Eigen::MatrixXi & triEF,
          const Eigen::VectorXi & triInnerEdges,
          Eigen::MatrixXd & currV,
          const Eigen::MatrixXi & EV,
          Eigen::VectorXi & VH,
          Eigen::VectorXi & HV,
          Eigen::VectorXi & HF,
          Eigen::VectorXi & FH,
          Eigen::VectorXi & nextH,
          Eigen::VectorXi & prevH,
          Eigen::VectorXi & twinH,
          std::vector<bool> & isParamVertex,
          std::vector<int> & HE2origEdges,
          std::vector<bool> & isParamHE,
          std::vector<int> & overlayFace2Tri,
          const double closeTolerance = 10e-6) {
        // create a map from the original edges to the half-edges
        std::vector<std::vector<int> > origEdges2HE(triEF.rows());
        for (int i = 0; i < HE2origEdges.size(); i++) {
          if (HE2origEdges[i] < 0)
            continue;
          origEdges2HE[HE2origEdges[i]].push_back(i);
        }
        std::fstream log;
        log.open("/home/kacper/log.txt", std::ios::out);

        std::vector< std::vector<int> > boundEdgesLeft(triEF.rows()), boundEdgesRight(triEF.rows());

        std::vector<bool> validHE(HV.rows(), true), validV(HE3D.size(), true), validF(FH.rows(), true);

        std::cout << "check before everything" << std::endl;
        checkMesh(HV, VH, HF, FH, twinH, nextH, prevH, isParamHE, validHE, validV, validF);

          //finding out vertex correspondence along twin edges of the original mesh by walking on boundaries
          std::vector<bool> Marked(HE2origEdges.size(), false);
          for (int i = 0; i < HV.rows(); i++) {
            if ((HE2origEdges[i] < 0) ||(Marked[i]))
              continue;
            //find the next beginning of a boundary
            int PrevOrig;
            int CurrEdge = i;
            do{
              PrevOrig = HE2origEdges[i];
              WalkBoundary(CurrEdge, nextH, twinH);
            }while(PrevOrig == HE2origEdges[CurrEdge]);

            //filling out strips of boundary with the respective attached original halfedges
            int BeginEdge = CurrEdge;
            vector<pair<int,int> > CurrEdgeCollect;
            do{
              CurrEdgeCollect.push_back(pair<int, int> (HE2origEdges[CurrEdge], CurrEdge));
              Marked[CurrEdge]=true;
              WalkBoundary(CurrEdge, nextH, twinH);
            }while (CurrEdge != BeginEdge);

            PrevOrig=-1000;
            bool In1 = true;
            for (int j = 0; j < CurrEdgeCollect.size(); j++) {
              if (CurrEdgeCollect[j].first != PrevOrig)
                In1 = boundEdgesLeft[CurrEdgeCollect[j].first].empty();

              if (In1)
                boundEdgesLeft[CurrEdgeCollect[j].first].push_back(CurrEdgeCollect[j].second);
              else
                boundEdgesRight[CurrEdgeCollect[j].first].push_back(CurrEdgeCollect[j].second);
              PrevOrig = CurrEdgeCollect[j].first;
            }
          }

        // match vertives
        std::vector<std::pair<int, int> > vertexMatches;
        find_maching_vertices(boundEdgesLeft, boundEdgesRight, nextH, HE3D,HV, vertexMatches);

        std::vector<int> transVertices(HE3D.size());
        graph_verification(boundEdgesLeft, boundEdgesRight, currV, nextH, prevH, twinH, HV, VH, HF, FH, vertexMatches, HE3D, transVertices, isParamVertex, validHE, validV, validF, isParamHE, HE2origEdges, log);


        //remove non-valid faces
        for(size_t i = 0; i < validF.size(); i++)
        {
          if(!validF[i] ||nextH(nextH(FH(i))) != FH(i))
            continue;

          validF[i] = false;
          int he = FH(i);

          if (twinH(nextH(he)) == -1 && twinH(he) != -1)
            he = nextH(he);

          if(twinH(nextH(he)) != -1 && twinH(he) != -1) {
            validHE[nextH(he)] = false;
            validHE[he] = false;
            twinH(twinH(nextH(he))) = twinH(he);
            twinH(twinH(he)) = twinH(nextH(he));
            VH(HV(nextH(he))) = twinH(he);
            VH(HV(he)) = twinH(nextH(he));

            if(isParamHE[he] || isParamHE[nextH(he)])
            {
              isParamHE[twinH(he)] = true;
              isParamHE[twinH(nextH(he))] = true;
            }
          }
          else if(twinH(he) == -1 && twinH(nextH(he)) != -1)
          {
            if(isParamHE[twinH(nextH(he))] || isParamHE[nextH(he)])
            {
              std::cout << "oups! " << std::endl;
            }

            int tmp = twinH(nextH(he));
            validHE[tmp] = false;
            validHE[nextH(he)] = false;
            prevH(nextH(tmp)) = he;
            nextH(prevH(tmp)) = he;

            VH(HV(nextH(he))) = nextH(tmp);
            VH(HV(tmp)) = he;

            nextH(he) = nextH(tmp);
            prevH(he) = prevH(tmp);
            HF(he) = HF(tmp);
          }
        }


        std::cout << "check after graph" << std::endl;
        checkMesh(HV, VH, HF, FH, twinH, nextH, prevH, isParamHE, validHE, validV, validF);


        //twinning up edges
        std::set<TwinFinder> twinning;
        for (int i = 0; i < twinH.rows(); i++) {
          if (twinH(i) != -1 || !validHE[i])
            continue;

          auto twinit = twinning.find(TwinFinder(0, HV(nextH(i)), HV(i)));
          if (twinit != twinning.end())  {
            twinH(twinit->index) = i;
            twinH(i) = twinit->index;
            log << "Twinning " << i << " and " << twinit->index << "\n";
            twinning.erase(*twinit);
          } else {
            twinning.insert(TwinFinder(i,HV(i), HV(nextH(i))));
          }
        }

        testUnmatchedTwins(twinH, nextH, HV, validHE, HE3D, HE2origEdges, log);

        std::cout << "check after twins" << std::endl;
        checkMesh(HV, VH, HF, FH, twinH, nextH, prevH, isParamHE, validHE, validV, validF);

        for (size_t i = 0; i < isParamHE.size(); i++){
          if(!validHE[i])
            continue;
          if (isParamHE[i])
            log <<"Hex edge "<< i << "\n";
          else
            log << "Triangle edge " << i << "\n";

          log << "Origin: " << HV(i) << "\n";
          log << "Prev: " << prevH(i) << "\n";
          log << "Next: " << nextH(i) <<"\n";
          log << "Twin: "<< twinH(i) << "\n";
          log << "Face: "<< HF(i) << "\n";
        }

//        for (size_t i = 0;i < validHE.size(); i++){
//          if (!isParamHE[i] && validHE[i]) {
//            joinFace(i, twinH, prevH, nextH, HF, FH, HV, VH, validHE, validV, validF);
//          }
//        }

        //unifying chains of edges
        //counting valences
//        std::vector<int> valences(HE3D.size(), 0);
//
//        for (size_t i = 0; i < validHE.size(); i++) {
//          if (validHE[i]) {
//            valences[HV(i)]++;
//            if (twinH(i) == -1)  //should account for the target as well
//              valences[HV(nextH(i))]++;
//          }
//        }
//
//        for (size_t i = 0; i < valences.size(); i++)
//          if (validV[i] && valences[i] < 2)
//            validV[i] = false;
//
//        for (size_t i = 0; i < valences.size(); i++)
//          if (validV[i] && valences[i] <= 2)
//            unifyEdges(VH(i), twinH, prevH, nextH, HF, FH, HV, VH, validHE, validV);


        std::cout << "check after degenerate faces removed" << std::endl;
        checkMesh(HV, VH, HF, FH, twinH, nextH, prevH, isParamHE, validHE, validV, validF);

        cleanMesh(validHE, validV, validF, twinH, prevH, nextH, currV, HF, FH, HV, VH, isParamVertex, HE2origEdges, isParamHE);
        log.close();
      }


      // Generates a new mesh from the parametrization over a given 2D grid (supported: square and hexagonal)

      // Input:

      // V           #V by 3, mesh vertices
      // F           #F by 3, vertex indices in face (it works only with triangles!)
      // EV          #E by 2, stores the edge description as pair of indices to vertices
      // FE          #F by 3, stores the Face-Edge relation
      // EF          #E by 2, stores the Edge-Face relation
      // InnerEdges, indices into EV of which edges are internal (not boundary)
      // UV           #V by 2, vertex coordinates in the parameter space
      // FUV          #F by 3, map between faces' vertices and coordinates in the parameter space

      // Output:

      // newV        #newV by 3, new mesh vertices
      // newD        #newF, number of vertices of each face
      // newF        #newF by 3, vertex indices in the new faces
      IGL_INLINE void generate_mesh(int N,
                                    const Eigen::MatrixXd & V,
                                    const Eigen::MatrixXi & F,
                                    const Eigen::MatrixXi & EV,
                                    const Eigen::MatrixXi & FE,
                                    const Eigen::MatrixXi & EF,
                                    const Eigen::VectorXi & innerEdges,
                                    const Eigen::MatrixXd & UV,
                                    const Eigen::MatrixXi & FUV,
                                    Eigen::MatrixXd & newV,
                                    Eigen::VectorXi & newD,
                                    Eigen::MatrixXi & newF)
                                    {
        Eigen::VectorXi VH; // map from vertices to the half-edges which start at these vertices
        Eigen::VectorXi HV; // map from half-edges to their source vertices
        Eigen::VectorXi HF; // map from half-edges to the corresponding faces
        Eigen::VectorXi FH; // map of each face with one of the half-edges
        Eigen::VectorXi nextH; // next half-edge
        Eigen::VectorXi prevH; // previous half-edge
        Eigen::VectorXi twinH; // twin half-edge (if -1 then an edge is a boundery edge)
        Eigen::MatrixXd currV; // current vertices (old and new vertices together ???)

        std::vector<bool> isParamVertex; // information of a given vertex is from the parametrization
        std::vector<int> HE2origEdges; // map between half-edges and original edges
        std::vector<bool> isParamHE; // information if a given half-edge is from the parametrization
        std::vector<int> overlayFace2Triangle; // triangle face ID or -1 when a face is unbounded

        std::vector<Point3D> HE3D;

        double minrange = (UV.colwise().maxCoeff() - UV.colwise().minCoeff()).minCoeff();
        // find the denominator for the  rational number representation
        int resolution = std::pow(10., std::ceil(std::log10(100000. / minrange)));

        if(F.cols() != 3)
          throw std::runtime_error("libhedra::generate_mesh: For now, it works only with triangular faces!");

        //creating an single-triangle arrangement
        //Intermediate growing DCEL
        for (int ti = 0; ti < F.rows(); ti++)
        {
          Arr_2 paramArr, triangleArr, overlayArr;

          std::cout << "Triangle: " << ti << std::endl;

          /* for all vertices of each face take UVs
           * and add edges of the face to the arrangment
           */
          for (int j = 0; j < 3; j++)
          {
            Eigen::RowVectorXd UV1 = UV.row(FUV(ti, j));
            Eigen::RowVectorXd UV2 = UV.row(FUV(ti, (j + 1) % 3));

            //avoid degenerate cases in non-bijective parametrizations
            if(paramCoord2texCoord(UV1) == paramCoord2texCoord(UV2))
              throw std::runtime_error("libhedra::generate_mesh: Only bijective parametrizations are supported, sorry!");
            Halfedge_handle he=CGAL::insert_non_intersecting_curve(triangleArr, Segment2D(paramCoord2texCoord(UV1), paramCoord2texCoord(UV2)));

            ArrEdgeData aed;
            aed.isParam = false;
            aed.origEdge = FE(ti, j);
            he->set_data(aed);
            he->twin()->set_data(aed);
          }

          // step up meta data for the faces in the arrangment of the mesh
          for (Face_iterator fi = triangleArr.faces_begin(); fi != triangleArr.faces_end(); fi++)
          {
            if (fi->is_unbounded())
              fi->data() = -1;
            else
              fi->data() = ti;
          }

          // generate a respective grid pattern
//          if (N == 4)
//            square_grid_pattern(UV, FUV, resolution, ti, paramArr);
//          else if (N == 6)
            tri_grid_pattern(UV, FUV, resolution, ti, paramArr);
          std::cout << "after " << std::endl;
//          else
//            throw std::runtime_error("libhedra::generate_mesh: Only the square and hexagonal grids are supported!");

          //Constructing the overlay arrangement
          Overlay_traits ot;
          overlay(triangleArr, paramArr, overlayArr, ot);


          //creating new halfedge structure from given mesh
          int formerNumVertices = currV.rows();
          int formerNumHalfedges = nextH.rows();
          int formerNumFaces = FH.rows();

          /*
           * collect information about new faces, vertices and edges and find the number of them which is used later to
           * resize the containers
           */
          int currFace = 0, currVertex = 0, currHalfedge = 0;
          for (Face_iterator fi = overlayArr.faces_begin(); fi != overlayArr.faces_end(); fi++)
          {
            if (fi->data() == -1)
              continue;  //one of the outer faces, i.e., unbounded

            overlayFace2Triangle.push_back(fi->data());
            fi->data() = formerNumFaces + currFace; // the id of a new face?
            currFace++;

            Ccb_halfedge_circulator hebegin = fi->outer_ccb();
            Ccb_halfedge_circulator heiterate = hebegin;
            //iterate over edges of the face
            do
            {
              if (heiterate->source()->data() < 0)  //new vertex
              {
                // information if an edge starting from this vertex is the parameter line
                isParamVertex.push_back(heiterate->source()->data() == PARAM_LINE_VERTEX); //source is the source vertex
                heiterate->source()->data() = formerNumVertices + currVertex; // the id of the new vertex ????
                currVertex++;
              }

              if (heiterate->data().newHalfedge < 0)  //new halfedge
              {
                // only edges which co-exist with original edges have values >= 0?????
                HE2origEdges.push_back(heiterate->data().origEdge);
                isParamHE.push_back(heiterate->data().isParam); // check if the edge is from the parameter lines
                heiterate->data().newHalfedge = formerNumHalfedges + currHalfedge; // the id of the new edge???
                currHalfedge++;
              }
              heiterate++;
            } while (heiterate != hebegin);
          }

          /*
           * build the initial mesh => original with seperated triangles plus the parametrization lines
           */
          currV.conservativeResize(currV.rows() + currVertex, 3);
          VH.conservativeResize(VH.size() + currVertex);
          HV.conservativeResize(HV.size() + currHalfedge);
          HF.conservativeResize(HF.size() + currHalfedge);
          FH.conservativeResize(FH.size() + currFace);
          nextH.conservativeResize(nextH.size() + currHalfedge);
          prevH.conservativeResize(prevH.size() + currHalfedge);
          twinH.conservativeResize(twinH.size() + currHalfedge);
          HE3D.resize(currV.rows());

          // we use the CGAL data structure in order to update the libheadra Egien containers which represent the mesh
          for (Face_iterator fi = overlayArr.faces_begin(); fi != overlayArr.faces_end(); fi++)
          {
            if (fi->data() == -1)
              continue;  //one of the outer faces, i.e., unbounded

            Ccb_halfedge_circulator hebegin = fi->outer_ccb();
            Ccb_halfedge_circulator heiterate = hebegin;
            //now assigning nexts and prevs
            do
            {
              // newHalfedge is the id of the halfEdge == formerNumHalfedges + currHalfedge
              // here we use the CGAL data structure in order to update the libredra data struture, which represent the halfedges
              nextH(heiterate->data().newHalfedge) = heiterate->next()->data().newHalfedge; // map the edge with its next
              prevH(heiterate->data().newHalfedge) = heiterate->prev()->data().newHalfedge; // map the edge with its prev
              twinH(heiterate->data().newHalfedge) = heiterate->twin()->data().newHalfedge; // map the edge with its twin

              if (heiterate->twin()->data().newHalfedge != -1) // no twin => the edge is a boundary edge
                twinH(heiterate->twin()->data().newHalfedge) = heiterate->data().newHalfedge;

              // as above but here we update source vertices of the half-edges
              HV(heiterate->data().newHalfedge) = heiterate->source()->data(); // map an edge with its source
              VH(heiterate->source()->data()) = heiterate->data().newHalfedge; // map a vertex with one of its edges
              HF(heiterate->data().newHalfedge) = fi->data(); // map an edge with the face
              FH(fi->data()) = heiterate->data().newHalfedge; // map the face with an edge
              heiterate++;
            } while (heiterate != hebegin);
          }

          /* constructing the actual vertices
           * The vertices from the arrangment are 2D and we need to project them on the 3D mesh,
           * and this is what is going on in this loop.
           */
          for (Vertex_iterator vi = overlayArr.vertices_begin(); vi != overlayArr.vertices_end(); vi++)
          {
            /* this are the vertices that are not sources of the parameter lines but their status seems to be -2,
             * otherwise we would have changed their data() to a current index. So what are they?
             * Do we have more vertices then we should?
             */
            if (vi->data() < 0)
              continue;

            Number BaryValues[3];
            Number Sum = 0;

            for (int i = 0; i < 3; i++)
            {
              //finding out barycentric coordinates
              // ti -- current face of the original mesh (the main loop)
              /* we start from + 1 and + 2 and not 0 and 1 because the weight of vertex v0 is the area/sum of the piece orthogonal to it
               * see the later loop where we multiply with BaryValues
               */
              Eigen::RowVectorXd UV2 = UV.row(FUV(ti, (i + 1) % 3));
              Eigen::RowVectorXd UV3 = UV.row(FUV(ti, (i + 2) % 3));

              Triangle2D t(vi->point(), paramCoord2texCoord(UV2),  paramCoord2texCoord(UV3));
              BaryValues[i] = t.area();
              Sum += BaryValues[i];
            }

            for (int i = 0; i < 3; i++)
              BaryValues[i] /= Sum;

            Point3D ENewPosition(0, 0, 0);
            //find the weighted position of the vertex inside the face, i.e., the 3D position of the vertex lifted to 3D
            for (int i = 0; i < 3; i++)
            {
              Point3D vertexCoord(V(F(ti, i), 0), V(F(ti, i), 1), V(F(ti, i), 2));

              ENewPosition = ENewPosition + (vertexCoord - CGAL::ORIGIN) * BaryValues[i];
            }
            HE3D[vi->data()] = ENewPosition;
            currV.row(vi->data()) = Eigen::RowVector3d(CGAL::to_double(ENewPosition.x()),
                                                       CGAL::to_double(ENewPosition.y()),
                                                       CGAL::to_double(ENewPosition.z()));
          }
        }

        //mesh unification
        stitch_boundaries2(HE3D, resolution, V, EF, innerEdges, currV, EV, VH, HV, HF, FH, nextH, prevH, twinH, isParamVertex, HE2origEdges, isParamHE, overlayFace2Triangle);

        //consolidation
        newV = currV;

        newD.conservativeResize(FH.rows());
        for (int i = 0; i < newD.size(); i++)
        {
          int ebegin = FH(i);
          int ecurr = ebegin;
          newD(i) = 0;
          do
          {
            newD(i)++;
            ecurr = nextH(ecurr);
          } while (ebegin != ecurr);
        }

        newF.conservativeResize(newD.rows(), newD.maxCoeff());
        for (int i = 0; i < newF.rows(); i++)
        {
          int ebegin = FH(i);
          int ecurr = ebegin;
          int currIndex = 0;
          do
          {
            newF(i, currIndex++) = HV(ecurr);
            ecurr = nextH(ecurr);
          } while (ebegin != ecurr);
        }
      }
    }
  }
}

#endif
  
  
