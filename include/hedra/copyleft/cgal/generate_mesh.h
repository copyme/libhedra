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

#include <CGAL/Root_of_traits.h>

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

        auto sqrt_3 = CGAL::make_sqrt(Number(3));

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
        Number coordMinX = *(std::min_element(coordsX.cbegin(), coordsX.cend())) /*- shift*/;
        Number coordMaxX = *(std::max_element(coordsX.cbegin(), coordsX.cend())) /*+ shift*/;

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
        //Number shift((int)std::round(CGAL::to_double(std::max(coordMaxX - coordMinX, *(std::max_element(coordsY.cbegin(), coordsY.cend())) - *(std::min_element(coordsY.cbegin(), coordsY.cend()))))));

//        inc = Number(3) / Number(2);
//        Number xShift = inc;
//        Number yShift = sqrt_3 / Number(2);
//
//        for (Number coordIndex = coordMinX; coordIndex <= coordMaxX; coordIndex += inc) {
//          //The line coord = coordIndex
//          lineCurves.emplace_back(Point2(coordIndex, coordMinY), Point2(coordIndex + xShift, coordMinY + yShift));
//          lineCurves.emplace_back(Point2(coordIndex, coordMinY), Point2(coordIndex + xShift, coordMinY - yShift));
//
//        }
        insert(paramArr, lineCurves.begin(), lineCurves.end());
      }


      void merge_edges_left(
          const int idx,
          std::vector<int> &leftHE,
          Eigen::VectorXi &HV,
          Eigen::VectorXi &HF,
          Eigen::VectorXi & FH,
          Eigen::VectorXi &nextH,
          Eigen::VectorXi &prevH,
          Eigen::VectorXi &twinH,
          std::set<int> & removedHE,
          std::set<int> & removedV) {

        prevH(nextH(leftHE[idx])) = prevH(leftHE[idx]);
        nextH(prevH(leftHE[idx])) = nextH(leftHE[idx]);

        int ebegin = twinH(prevH(leftHE[idx]));
        int ecurr = ebegin;
        int counter = 0;
        do {
          if (counter > 30)
            throw std::runtime_error("libhedra:stitch_boundaries: This should not happened!\n "
                                     "Verify if your UV coordinates match at the cut and if so then report a bug at: https://github.com/avaxman/libhedra/issues");
          if(ecurr == -1)
            break;

          HV(ecurr) = HV(nextH(leftHE[idx]));
          ecurr = twinH(prevH(ecurr));
          counter++;
        } while (ebegin != ecurr);


        removedHE.insert(leftHE[idx]);
        removedV.insert(HV(leftHE[idx]));
        FH(HF(leftHE[idx])) = nextH(leftHE[idx]);

        leftHE.erase(leftHE.begin() + idx);
      }


      void kill_triangle_left(
          const int idx,
          std::vector<int> &leftHE,
          Eigen::MatrixXd & currV,
          Eigen::VectorXi &HV,
          Eigen::VectorXi & HF,
          Eigen::VectorXi & FH,
          Eigen::VectorXi &nextH,
          Eigen::VectorXi &prevH,
          Eigen::VectorXi &twinH,
          std::set<int> & removedHE,
          std::set<int> & removedV) {

        //removed tri edges
        removedHE.insert(leftHE[idx]);
        removedHE.insert(nextH(leftHE[idx]));
        removedHE.insert(prevH(leftHE[idx]));


        if(twinH(prevH(leftHE[idx])) != -1 && twinH(nextH(leftHE[idx])) != -1)
        {
          return;
          std::cout << "middle tri" << std::endl;
          Eigen::Vector3d mid = currV.row(HV(leftHE[idx])) + (currV.row(HV(nextH(leftHE[idx]))) - currV.row(HV(leftHE[idx]))) / 2.;

          removedV.insert(HV(leftHE[idx]));
          removedV.insert(HV(nextH(leftHE[idx])));
          //bottom
          prevH(nextH(twinH(prevH(leftHE[idx])))) = prevH(twinH(prevH(leftHE[idx])));
          nextH(prevH(twinH(prevH(leftHE[idx])))) = nextH(twinH(prevH(leftHE[idx])));

          currV.row(HV(twinH(nextH(leftHE[idx])))) = mid;

          removedHE.insert(twinH(prevH(leftHE[idx])));
          removedHE.insert(twinH(nextH(leftHE[idx])));

          FH(HF(twinH(nextH(leftHE[idx])))) =  nextH(twinH(nextH(leftHE[idx])));
          FH(HF(twinH(prevH(leftHE[idx])))) =  nextH(twinH(prevH(leftHE[idx])));

          int ebegin = twinH(prevH(leftHE[idx]));
          int ecurr = ebegin;
          int counter = 0;
          do {
            if (counter > 30)
              throw std::runtime_error("libhedra:stitch_boundaries: This should not happened!\n "
                                       "Verify if your UV coordinates match at the cut and if so then report a bug at: https://github.com/avaxman/libhedra/issues");
            HV(ecurr) = HV(nextH(twinH(prevH(leftHE[idx]))));
            ecurr = twinH(prevH(ecurr));
            counter++;
          } while (ecurr != -1);

          // top
          prevH(nextH(twinH(nextH(leftHE[idx])))) = prevH(twinH(nextH(leftHE[idx])));
          nextH(prevH(twinH(nextH(leftHE[idx])))) = nextH(twinH(nextH(leftHE[idx])));

          ebegin = nextH(twinH(nextH(leftHE[idx])));
          ecurr = ebegin;
          counter = 0;
          do {
            if (counter > 30)
              throw std::runtime_error("libhedra:stitch_boundaries: This should not happened!\n "
                                       "Verify if your UV coordinates match at the cut and if so then report a bug at: https://github.com/avaxman/libhedra/issues");
            HV(ecurr) = HV(nextH(twinH(prevH(leftHE[idx]))));
            if(twinH(ecurr) == -1)
              break;
            ecurr = nextH(twinH(ecurr));
            counter++;
          } while (ecurr != ebegin);

        }
        // top triangle
        else if(twinH(prevH(leftHE[idx])) != -1)
        {
          std::cout << "top tri" << std::endl;
          prevH(nextH(twinH(prevH(leftHE[idx])))) = prevH(twinH(prevH(leftHE[idx])));
          nextH(prevH(twinH(prevH(leftHE[idx])))) = nextH(twinH(prevH(leftHE[idx])));

          removedHE.insert(twinH(prevH(leftHE[idx])));

          FH(HF(twinH(prevH(leftHE[idx])))) =  nextH(twinH(prevH(leftHE[idx])));

          int ebegin = twinH(prevH(leftHE[idx]));
          int ecurr = ebegin;
          int counter = 0;
          do {
            if (counter > 30)
              throw std::runtime_error("libhedra:stitch_boundaries: This should not happened!\n "
                                       "Verify if your UV coordinates match at the cut and if so then report a bug at: https://github.com/avaxman/libhedra/issues");
            HV(ecurr) = HV(nextH(leftHE[idx]));
            ecurr = twinH(prevH(ecurr));
            counter++;
          } while (ecurr != -1);

          ebegin = nextH(twinH(prevH(leftHE[idx])));
          ecurr = ebegin;
          counter = 0;
          do {
            if (counter > 30)
              throw std::runtime_error("libhedra:stitch_boundaries: This should not happened!\n "
                                       "Verify if your UV coordinates match at the cut and if so then report a bug at: https://github.com/avaxman/libhedra/issues");
            HV(ecurr) = HV(nextH(leftHE[idx]));
            if(twinH(ecurr) == -1)
              break;
            ecurr = nextH(twinH(ecurr));
            counter++;
          } while (ecurr != ebegin);
        }
        else {
          std::cout << "bottom tri" << std::endl;

          prevH(nextH(twinH(nextH(leftHE[idx])))) = prevH(twinH(nextH(leftHE[idx])));
          nextH(prevH(twinH(nextH(leftHE[idx])))) = nextH(twinH(nextH(leftHE[idx])));

          removedHE.insert(twinH(nextH(leftHE[idx])));

          FH(HF(twinH(nextH(leftHE[idx])))) =  nextH(twinH(nextH(leftHE[idx])));

          int ebegin = nextH(twinH(nextH(leftHE[idx])));
          int ecurr = ebegin;
          int counter = 0;
          do {
            if (counter > 30)
              throw std::runtime_error("libhedra:stitch_boundaries: This should not happened!\n "
                                       "Verify if your UV coordinates match at the cut and if so then report a bug at: https://github.com/avaxman/libhedra/issues");
            HV(ecurr) = HV(leftHE[idx]);
            if(twinH(ecurr) == -1)
              break;
            ecurr = nextH(twinH(ecurr));
            counter++;
          } while (ecurr != -1);

          ebegin = twinH(nextH(leftHE[idx]));
          ecurr = ebegin;
          counter = 0;
          do {
            if (counter > 30)
              throw std::runtime_error("libhedra:stitch_boundaries: This should not happened!\n "
                                       "Verify if your UV coordinates match at the cut and if so then report a bug at: https://github.com/avaxman/libhedra/issues");
            HV(ecurr) = HV(leftHE[idx]);
            ecurr = twinH(prevH(ecurr));
            counter++;
          } while (ecurr != -1);
        }

        leftHE.erase(leftHE.begin() + idx);
      }

      bool edge_reduction(std::vector<int> & leftHE,
                          std::vector<int> & rightHE,
                        Eigen::MatrixXd & currV,
                        Eigen::VectorXi & HV,
                        Eigen::VectorXi & HF,
                          Eigen::VectorXi & FH,
                          Eigen::VectorXi & nextH,
                        Eigen::VectorXi & prevH,
                        Eigen::VectorXi & twinH,
                        std::set<int> & removedHE,
                        std::set<int> & removedV,
                        const double closeTolerance = 10e-30) {
        bool status = false;


//          for (int j = 0; j < rightHE.size(); j++) {
//            std::multimap<int, int> right2left;
//            Eigen::Vector3d M = currV.row(HV(rightHE[j])) - currV.row(HV(nextH(rightHE[j])));
//            Eigen::Vector3d B = currV.row(HV(nextH(rightHE[j])));
//
//            for (int k = 0; k < leftHE.size(); k++) {
//              Eigen::Vector3d N = (currV.row(HV(nextH(leftHE[k]))) - currV.row(HV(leftHE[k]))) / 2.;
//              // take middle point to avoid mismachings
//              Eigen::Vector3d P = currV.row(HV(leftHE[k])).transpose() + N;
//              double t = M.dot((P - B)) / (M.dot(M));
//              std::cout << t << std::endl;
//              if (t > 0. && t < 1.)
//                right2left.insert(std::pair<int, int>(j, k));
//            }
//
//            typedef std::multimap<int, int>::iterator MMAPIterator;
//            std::pair<MMAPIterator, MMAPIterator> result = right2left.equal_range(j);
//            int idx = -1;
//            double N = 300000;
//
//            if (std::distance(result.first, result.second) < 2)
//              continue;
//
//            for (auto it = result.first; it != result.second; it++)
//            {
//              double tmp = (currV.row(HV(nextH(leftHE[it->second]))) - currV.row(HV(leftHE[it->second]))).norm();
//              if(tmp < N)
//              {
//                idx = it->second;
//                N = tmp;
//              }
//            }
//
//            if(nextH(nextH(leftHE[idx])) == prevH(leftHE[idx])) {
//            kill_triangle_left(idx, leftHE, HV, HF, FH, nextH, prevH, twinH, removedHE, removedV);
//             return true;
//            }
//           else {
//              merge_edges_left(idx, leftHE, HV, HF, FH, nextH, prevH, twinH, removedHE, removedV);
//              return true;
//            }
//
//            status = true;
//          }
//          return status;


            int idx = -1;
            double N = 300000;

        for (int i = 0; i < leftHE.size(); i++)
        {
          double tmp = (currV.row(HV(nextH(leftHE[i]))) - currV.row(HV(leftHE[i]))).norm();
          if(tmp < N)
          {
            idx = i;
            N = tmp;
          }
        }

        if(nextH(nextH(leftHE[idx])) == prevH(leftHE[idx])) {
          kill_triangle_left(idx, leftHE, currV, HV, HF, FH, nextH, prevH, twinH, removedHE, removedV);
          return true;
        }
        else {
         // merge_edges_left(idx, leftHE, HV, HF, FH, nextH, prevH, twinH, removedHE, removedV);
          return true;
        }
      }


      // Connects disconnected pieces of the mesh

      // Input:
      // EF               #E by 2, stores the Edge-Face relation of the original triangle mesh
      // triInnerEdges    #E, stores information about inner edges of the original triagle mesh i.e., no boundary edges
      // currV            current vertices (old and new vertices together ???)
      // VH               map from vertices to the half-edges which start at these vertices
      // HV               map from half-edges to their source vertices
      // HF               map from half-edges to the corresponding faces
      // FH               map of each face with one of the half-edges
      // nextH            next half-edge
      // prevH            previous half-edge
      // twinH            twin half-edge (if -1 then an edge is a boundery edge)
      // isParamVertex    information of a given vertex is from the parametrization
      // HE2origEdges     map between half-edges and original edges
      // ParamHE          information if a given half-edge is from the parametrization
      // overlayFace2Tri  triangle face ID or -1 when a face is unbounded
      // closeTolerance   value with which the new vertices are stiched

      // Output:
      IGL_INLINE void stitch_boundaries(
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
                                        const double closeTolerance = 10e-30) {
        // create a map from the original edges to the half-edges
        std::vector<std::vector<int> > origEdges2HE(triEF.rows());
        for (int i = 0; i < HE2origEdges.size(); i++) {
          if (HE2origEdges[i] < 0)
            continue;
          origEdges2HE[HE2origEdges[i]].push_back(i);
        }
        Eigen::VectorXi oldHF = HF;

        // bind faces with parameter lines
        for (size_t k = 0; k < FH.size(); k++) {
          int ebegin = FH(k);
          int ecurr = ebegin;
          do {
            if (isParamHE[ecurr]) {
              FH(k) = ecurr;
              break;
            }
            ecurr = nextH(ecurr);
          } while (ebegin != ecurr);
        }

        //for every original inner edge, stitching up boundary (original boundary edges don't have any action item)
        for (int i = 0; i < triInnerEdges.size(); i++) {
          //first sorting to left and right edges according to faces
          int currEdge = triInnerEdges(i);
          int leftFace = triEF(currEdge, 1);
          int rightFace = triEF(currEdge, 0);

          std::vector<int> leftHE, rightHE;
          for (size_t k = 0; k < origEdges2HE[currEdge].size(); k++) {
            if (overlayFace2Tri[oldHF(origEdges2HE[currEdge][k])] == leftFace) {
              leftHE.push_back(origEdges2HE[currEdge][k]);
            } else if (overlayFace2Tri[oldHF(origEdges2HE[currEdge][k])] == rightFace) {
              rightHE.push_back(origEdges2HE[currEdge][k]);
            } else
              throw std::runtime_error(
                  "libhedra:stitch_boundaries: This should not happened! Report a bug at: https://github.com/avaxman/libhedra/issues");
          }


          //sort left and right edges
          Eigen::RowVector3d refV = V.row(EV(currEdge, 0)) - (V.row(EV(currEdge, 1)) - V.row(EV(currEdge, 0))) * 2.;
          Point3D ref(Number(refV(0)), Number(refV(1)), Number(refV(2)));

          std::stable_sort(leftHE.begin(), leftHE.end(),
                           [&ref, &HE3D, &HV](const int &a, const int &b) -> bool {
                             Point3D A = HE3D[HV(a)];
                             Point3D B = HE3D[HV(b)];
                             return CGAL::has_smaller_distance_to_point(ref, A, B);
                           }
          );

          std::stable_sort(rightHE.begin(), rightHE.end(),
                           [&ref, &HE3D, &HV](const int &a, const int &b) -> bool {
                             Point3D A = HE3D[HV(a)];
                             Point3D B = HE3D[HV(b)];
                             return CGAL::has_smaller_distance_to_point(ref, A, B);

                           }
          );

          Point3D A = HE3D[HV(leftHE[0])];
          Point3D B = HE3D[HV(rightHE[0])];
          // swap if the right is really left and vice versa
          if(CGAL::has_smaller_distance_to_point(ref, B, A))
          {
            for(size_t k = 0; k < leftHE.size(); k++)
              std::swap(leftHE[k], rightHE[k]);
          }

          // garbage collector
          std::set<int> removedHE, removedV;


          if(rightHE.size() > leftHE.size())
            continue;

          while(leftHE.size() > rightHE.size()) {
            std::cout << leftHE.size() << " before " << rightHE.size() << " l " << leftFace << " r " << rightFace  << std::endl;
            edge_reduction(leftHE, rightHE, currV, HV, HF, FH, nextH, prevH, twinH, removedHE, removedV, closeTolerance);
            std::cout << leftHE.size() << " after " << rightHE.size() << std::endl;

          }

          for (size_t j = 0; j < leftHE.size(); j++) {
            if (!(isParamHE[leftHE[j]] && isParamHE[rightHE[j]])) {
              int ebegin = rightHE[j];
              int ecurr = ebegin;
              int counter = 0;
              do {
                if (counter > 30)
                  throw std::runtime_error("libhedra:stitch_boundaries: This should not happened!\n "
                                           "Verify if your UV coordinates match at the cut and if so then report a bug at: https://github.com/avaxman/libhedra/issues");
                HF(ecurr) = HF(leftHE[j]);
                ecurr = nextH(ecurr);
                counter++;
              } while (ebegin != ecurr);
            }
          }

          //find maching source vertices from left to right
          for (size_t j = 1; j < leftHE.size(); j++) {
            /* consider a case when both edges from the pair are not parameter lines, i.e.,
             * the edges have to be removed.
             */
            // remove the pair from the orphanage
            if (!isParamHE[leftHE[j]] && !isParamHE[rightHE[j]] && !isParamVertex[HV(leftHE[j])]) {
              nextH(prevH(leftHE[j])) = nextH(nextH(rightHE[j]));
              prevH(nextH(nextH(rightHE[j]))) = prevH(leftHE[j]);

              nextH(twinH(nextH(rightHE[j]))) = nextH(twinH(prevH(leftHE[j])));
              prevH(nextH(twinH(prevH(leftHE[j])))) = twinH(nextH(rightHE[j]));

              // garbage collector
              removedHE.insert(twinH(prevH(leftHE[j])));
              removedHE.insert(nextH(rightHE[j]));

              // stich twins
              twinH(prevH(leftHE[j])) = twinH(nextH(rightHE[j]));
              twinH(twinH(nextH(rightHE[j]))) = prevH(leftHE[j]);

              // garbage collector
              removedHE.insert(leftHE[j]);
              removedHE.insert(rightHE[j]);
              removedV.insert(HV(leftHE[j]));
              removedV.insert(HV(nextH(rightHE[j])));

              //ensure that a face is not refered to a removed edge
              FH(HF(twinH(prevH(leftHE[j])))) = twinH(nextH(rightHE[j]));
              FH(HF(nextH(rightHE[j]))) = prevH(leftHE[j]);
            }
              //rotated cross case
            else if (!isParamHE[leftHE[j]] && !isParamHE[rightHE[j]] && isParamVertex[HV(leftHE[j])]) {

              // stitch f0
              nextH(prevH(leftHE[j])) = nextH(rightHE[j]);
              prevH(nextH(rightHE[j])) = prevH(leftHE[j]);

              //garnage collector
              removedHE.insert(leftHE[j]);
              removedHE.insert(rightHE[j]);
              removedV.insert(HV(nextH(rightHE[j])));

              // update vertex
              int ebegin = rightHE[j];
              int ecurr = ebegin;
              int counter = 0;
              do {
                if (counter > 30)
                  throw std::runtime_error("libhedra:stitch_boundaries: This should not happened!\n "
                                           "Verify if your UV coordinates match at the cut and if so then report a bug at: https://github.com/avaxman/libhedra/issues");
                HV(nextH(ecurr)) = HV(leftHE[j]);
                ecurr = twinH(nextH(ecurr));
                counter++;
              } while (twinH(nextH(ecurr)) != -1);

              // stitch f1
              prevH(twinH(prevH(twinH(prevH(leftHE[j]))))) = twinH(nextH(twinH(nextH(rightHE[j]))));
              nextH(twinH(nextH(twinH(nextH(rightHE[j]))))) = twinH(prevH(twinH(prevH(leftHE[j]))));

            } else if (isParamHE[leftHE[j]] && isParamHE[rightHE[j]]) {
              twinH(leftHE[j]) = rightHE[j];
              twinH(rightHE[j]) = leftHE[j];
              removedV.insert(HV(nextH(rightHE[j])));
              HV(nextH(rightHE[j])) = HV(leftHE[j]);
              HV(nextH(twinH(nextH(rightHE[j])))) = HV(leftHE[j]);
            }
            else
              throw std::runtime_error(
                  "libhedra:stitch_boundaries: This should not happened! Report a bug at: https://github.com/avaxman/libhedra/issues");
          }

          // merge edge ends
          if (!isParamHE[leftHE[0]]) {
            nextH(prevH(leftHE[0])) = nextH(rightHE[0]); //
            prevH(nextH(rightHE[0])) = prevH(leftHE[0]);
            //garbage collector
            removedHE.insert(leftHE[0]);
            removedHE.insert(rightHE[0]);
            HV(nextH(rightHE[0])) = HV(leftHE[0]);

          } else if (isParamHE[leftHE[0]]) {
            twinH(leftHE[0]) = rightHE[0];
            twinH(rightHE[0]) = leftHE[0];
            removedV.insert(HV(nextH(rightHE[0])));
            HV(nextH(rightHE[0])) = HV(leftHE[0]);

          } else
            throw std::runtime_error(
                "libhedra:stitch_boundaries: This should not happened! Report a bug at: https://github.com/avaxman/libhedra/issues");

          int last = leftHE.size() - 1;
          if (!isParamHE[rightHE[last]]) {
            nextH(prevH(rightHE[last])) = nextH(leftHE[last]);
            prevH(nextH(leftHE[last])) = prevH(rightHE[last]);
            removedHE.insert(rightHE[last]);
            removedHE.insert(leftHE[last]);

            HV(rightHE[last]) = HV(nextH(leftHE[last]));
          } else if (isParamHE[rightHE[last]]) {
            twinH(leftHE[last]) = rightHE[last];
            twinH(rightHE[last]) = leftHE[last];

            removedV.insert(HV(rightHE[last]));
            HV(rightHE[last]) = HV(nextH(leftHE[last]));
          } else
            throw std::runtime_error(
                "libhedra:stitch_boundaries: This should not happened! Report a bug at: https://github.com/avaxman/libhedra/issues");

          /* removed virtual objects
           *
           */

          //edges
          for (auto he = removedHE.rbegin(); he != removedHE.rend(); he++) {
            // skipe a rediscovered useless edge, well normally it should not happen
            if (*he > HF.rows()) {
              continue;
            }
            int numRows = HE2origEdges.size() - 1;

            if (*he < numRows) {
              HF.segment(*he, numRows - *he) = HF.segment(*he + 1, numRows - *he).eval();
              oldHF.segment(*he, numRows - *he) = oldHF.segment(*he + 1, numRows - *he).eval();
              HV.segment(*he, numRows - *he) = HV.segment(*he + 1, numRows - *he).eval();
              nextH.segment(*he, numRows - *he) = nextH.segment(*he + 1, numRows - *he).eval();
              twinH.segment(*he, numRows - *he) = twinH.segment(*he + 1, numRows - *he).eval();
              prevH.segment(*he, numRows - *he) = prevH.segment(*he + 1, numRows - *he).eval();
            }

            HF.conservativeResize(numRows);
            oldHF.conservativeResize(numRows);
            HV.conservativeResize(numRows);
            nextH.conservativeResize(numRows);
            twinH.conservativeResize(numRows);
            prevH.conservativeResize(numRows);

            HE2origEdges.erase(HE2origEdges.cbegin() + (*he));
            isParamHE.erase(isParamHE.cbegin() + (*he));

            //update IDs
            for (int k = 0; k < FH.rows(); k++)
              if (FH(k) > *he)
                FH(k)--;

            for (int k = 0; k < VH.rows(); k++)
              if (VH(k) > *he)
                VH(k)--;
              else if (VH(k) == *he && removedHE.find(*he) == removedHE.end())
                throw std::runtime_error("libhedra:stitch_boundaries: bad ref. VH! Report a bug at: https://github.com/avaxman/libhedra/issues");

            for (int k = 0; k < nextH.rows(); k++)
              if (nextH(k) > *he)
                nextH(k)--;
              else if (nextH(k) == *he && removedHE.find(*he) == removedHE.end())
                throw std::runtime_error("libhedra:stitch_boundaries: bad ref. nextH! Report a bug at: https://github.com/avaxman/libhedra/issues");

            for (int k = 0; k < prevH.rows(); k++)
              if (prevH(k) > *he)
                prevH(k)--;
              else if (prevH(k) == *he && removedHE.find(*he) == removedHE.end())
                throw std::runtime_error("libhedra:stitch_boundaries: bad ref. prevH! Report a bug at: https://github.com/avaxman/libhedra/issues");

            for (int k = 0; k < twinH.rows(); k++)
              if (twinH(k) > *he)
                twinH(k)--;
              else if (twinH(k) == *he && removedHE.find(*he) == removedHE.end())
                throw std::runtime_error("libhedra:stitch_boundaries: bad ref. twinH! Report a bug at: https://github.com/avaxman/libhedra/issues");

            for (size_t k = 0; k < origEdges2HE.size(); k++) {
              auto it = std::find(origEdges2HE[k].begin(), origEdges2HE[k].end(), *he);
              if (it != origEdges2HE[k].end())
                origEdges2HE[k].erase(it);
              for (size_t h = 0; h < origEdges2HE[k].size(); h++)
                if (origEdges2HE[k][h] > *he)
                  origEdges2HE[k][h]--;
            }
          }

          //vertices
          for (auto vi = removedV.rbegin(); vi != removedV.rend(); vi++) {
            //remove the row
            int numRows = currV.rows() - 1;
            if (*vi < numRows) {
              currV.block(*vi, 0, numRows - *vi, 3) = currV.block(*vi + 1, 0, numRows - *vi, 3).eval();
              VH.segment(*vi, numRows - *vi) = VH.segment(*vi + 1, numRows - *vi).eval();
            }
            currV.conservativeResize(numRows, 3);
            VH.conservativeResize(numRows);
            isParamVertex.erase(isParamVertex.begin() + *vi);

            HE3D.erase(HE3D.begin() + *vi);
            //update IDs
            for (int k = 0; k < HV.rows(); k++) {
              if (HV(k) > *vi)
                HV(k)--;
            }
          }

        }

        //removed unreferenced faces
        std::vector<int> hitFaces(FH.rows(), 0);
        for (int k = 0; k < HF.rows(); k++)
          hitFaces[HF(k)]++;

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

        //merge vertives
        for(int k = 0; k < HV.rows(); k++)
        {
          for(int j = 0; j < HV.rows(); j++)
          {
            if(k == j)
              continue;
            if((currV.row(HV(j)) - currV.row(HV(k))).norm() < closeTolerance)
              HV(j) = HV(k);
          }
        }

        //removed unreferenced vertices
        std::vector<int> hitVers(VH.rows(), 0);
        for (int k = 0; k < HV.rows(); k++)
          hitVers[HV(k)]++;

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
          const double closeTolerance = 10e-30) {
        // create a map from the original edges to the half-edges
        std::vector<std::vector<int> > origEdges2HE(triEF.rows());
        for (int i = 0; i < HE2origEdges.size(); i++) {
          if (HE2origEdges[i] < 0)
            continue;
          origEdges2HE[HE2origEdges[i]].push_back(i);
        }
        Eigen::VectorXi oldHF = HF;


        //for every original inner edge, stitching up boundary (original boundary edges don't have any action item)
        for (int i = 0; i < triInnerEdges.size(); i++) {
          //first sorting to left and right edges according to faces
          int currEdge = triInnerEdges(i);
          int leftFace = triEF(currEdge, 1);
          int rightFace = triEF(currEdge, 0);

          std::cout << leftFace << " " << rightFace << std::endl;

          std::vector<int> leftHE, rightHE;
          for (size_t k = 0; k < origEdges2HE[currEdge].size(); k++) {
            if (overlayFace2Tri[oldHF(origEdges2HE[currEdge][k])] == leftFace) {
              leftHE.push_back(origEdges2HE[currEdge][k]);
            } else if (overlayFace2Tri[oldHF(origEdges2HE[currEdge][k])] == rightFace) {
              rightHE.push_back(origEdges2HE[currEdge][k]);
            } else
              throw std::runtime_error(
                  "libhedra:stitch_boundaries: This should not happened! Report a bug at: https://github.com/avaxman/libhedra/issues");
          }

          //sort left and right edges
          Eigen::RowVector3d refV = V.row(EV(currEdge, 0)) - (V.row(EV(currEdge, 1)) - V.row(EV(currEdge, 0))) * 2.;
          Point3D ref(Number(refV(0)), Number(refV(1)), Number(refV(2)));

          std::stable_sort(leftHE.begin(), leftHE.end(),
                           [&ref, &HE3D, &HV](const int &a, const int &b) -> bool {
                             Point3D A = HE3D[HV(a)];
                             Point3D B = HE3D[HV(b)];
                             return CGAL::has_smaller_distance_to_point(ref, A, B);
                           }
          );

          std::stable_sort(rightHE.begin(), rightHE.end(),
                           [&ref, &HE3D, &HV](const int &a, const int &b) -> bool {
                             Point3D A = HE3D[HV(a)];
                             Point3D B = HE3D[HV(b)];
                             return CGAL::has_smaller_distance_to_point(ref, A, B);

                           }
          );

          Point3D A = HE3D[HV(leftHE[0])];
          Point3D B = HE3D[HV(rightHE[0])];
          // swap if the right is really left and vice versa
          if (CGAL::has_smaller_distance_to_point(ref, B, A)) {
            for (size_t k = 0; k < leftHE.size(); k++)
              std::swap(leftHE[k], rightHE[k]);
          }

          // garbage collector
          std::set<int> removedHE, removedV;

          //if(rightHE.size() != leftHE.size())
         //   continue;

          if(rightHE.size() == leftHE.size() - 1) {
            edge_reduction(leftHE, rightHE, currV, HV, HF, FH, nextH, prevH, twinH, removedHE, removedV, closeTolerance);
          }
          if(rightHE.size() < leftHE.size())
            continue;

//
//          if(rightHE.size() > leftHE.size()) {
//            edge_reduction(rightHE, leftHE, currV, HV, HF, FH, nextH, prevH, twinH, removedHE, removedV, closeTolerance);
//          }
//
//          if(rightHE.size() > leftHE.size())
//            continue;

          //find maching source vertices from left to right
          for (size_t j = 1; j < leftHE.size(); j++) {
            twinH(leftHE[j]) = rightHE[j];
            twinH(rightHE[j]) = leftHE[j];
            removedV.insert(HV(nextH(rightHE[j])));
            HV(nextH(rightHE[j])) = HV(leftHE[j]);
            if(twinH(nextH(rightHE[j])) != -1)
              HV(nextH(twinH(nextH(rightHE[j])))) = HV(leftHE[j]);
          }
          HV(rightHE.back()) = HV(nextH(leftHE.back()));

          /* removed virtual objects
           *
           */

          //edges
          for (auto he = removedHE.rbegin(); he != removedHE.rend(); he++) {
            // skipe a rediscovered useless edge, well normally it should not happen
            if (*he > HF.rows()) {
              continue;
            }
            int numRows = HE2origEdges.size() - 1;

            if (*he < numRows) {
              HF.segment(*he, numRows - *he) = HF.segment(*he + 1, numRows - *he).eval();
              oldHF.segment(*he, numRows - *he) = oldHF.segment(*he + 1, numRows - *he).eval();
              HV.segment(*he, numRows - *he) = HV.segment(*he + 1, numRows - *he).eval();
              nextH.segment(*he, numRows - *he) = nextH.segment(*he + 1, numRows - *he).eval();
              twinH.segment(*he, numRows - *he) = twinH.segment(*he + 1, numRows - *he).eval();
              prevH.segment(*he, numRows - *he) = prevH.segment(*he + 1, numRows - *he).eval();
            }

            HF.conservativeResize(numRows);
            oldHF.conservativeResize(numRows);
            HV.conservativeResize(numRows);
            nextH.conservativeResize(numRows);
            twinH.conservativeResize(numRows);
            prevH.conservativeResize(numRows);

            HE2origEdges.erase(HE2origEdges.cbegin() + (*he));
            isParamHE.erase(isParamHE.cbegin() + (*he));

            //update IDs
            for (int k = 0; k < FH.rows(); k++)
              if (FH(k) > *he)
                FH(k)--;

            for (int k = 0; k < VH.rows(); k++)
              if (VH(k) > *he)
                VH(k)--;
              else if (VH(k) == *he && removedHE.find(*he) == removedHE.end())
                throw std::runtime_error("libhedra:stitch_boundaries: bad ref. VH! Report a bug at: https://github.com/avaxman/libhedra/issues");

            for (int k = 0; k < nextH.rows(); k++)
              if (nextH(k) > *he)
                nextH(k)--;
              else if (nextH(k) == *he && removedHE.find(*he) == removedHE.end())
                throw std::runtime_error("libhedra:stitch_boundaries: bad ref. nextH! Report a bug at: https://github.com/avaxman/libhedra/issues");

            for (int k = 0; k < prevH.rows(); k++)
              if (prevH(k) > *he)
                prevH(k)--;
              else if (prevH(k) == *he && removedHE.find(*he) == removedHE.end())
                throw std::runtime_error("libhedra:stitch_boundaries: bad ref. prevH! Report a bug at: https://github.com/avaxman/libhedra/issues");

            for (int k = 0; k < twinH.rows(); k++)
              if (twinH(k) > *he)
                twinH(k)--;
              else if (twinH(k) == *he && removedHE.find(*he) == removedHE.end())
                throw std::runtime_error("libhedra:stitch_boundaries: bad ref. twinH! Report a bug at: https://github.com/avaxman/libhedra/issues");

            for (size_t k = 0; k < origEdges2HE.size(); k++) {
              auto it = std::find(origEdges2HE[k].begin(), origEdges2HE[k].end(), *he);
              if (it != origEdges2HE[k].end())
                origEdges2HE[k].erase(it);
              for (size_t h = 0; h < origEdges2HE[k].size(); h++)
                if (origEdges2HE[k][h] > *he)
                  origEdges2HE[k][h]--;
            }
          }

          //vertices
          for (auto vi = removedV.rbegin(); vi != removedV.rend(); vi++) {
            //remove the row
            int numRows = currV.rows() - 1;
            if (*vi < numRows) {
              currV.block(*vi, 0, numRows - *vi, 3) = currV.block(*vi + 1, 0, numRows - *vi, 3).eval();
              VH.segment(*vi, numRows - *vi) = VH.segment(*vi + 1, numRows - *vi).eval();
            }
            currV.conservativeResize(numRows, 3);
            VH.conservativeResize(numRows);
            isParamVertex.erase(isParamVertex.begin() + *vi);

            HE3D.erase(HE3D.begin() + *vi);
            //update IDs
            for (int k = 0; k < HV.rows(); k++) {
              if (HV(k) > *vi)
                HV(k)--;
            }
          }

        }

        //removed unreferenced faces
        std::vector<int> hitFaces(FH.rows(), 0);
        for (int k = 0; k < HF.rows(); k++)
          hitFaces[HF(k)]++;

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
  
  
