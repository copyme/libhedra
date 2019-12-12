#include <igl/opengl/glfw/Viewer.h>
#include <igl/readOBJ.h>
#include <hedra/polygonal_edge_topology.h>
#include <hedra/copyleft/cgal/generate_mesh.h>
#include <hedra/polygonal_write_OFF.h>
#include <hedra/dual_mesh.h>

Eigen::MatrixXd V, newV;
Eigen::MatrixXi F, newF, FTC;
Eigen::VectorXi newD;

bool key_down(igl::opengl::glfw::Viewer& viewer, unsigned char key, int modifiers)
{
  return true;
}

int main(int argc, char *argv[])
{
  using namespace Eigen;
  using namespace std;

  Eigen::MatrixXd TC;
  Eigen::MatrixXd N;
  Eigen::MatrixXd FN, FEs;
  Eigen::MatrixXi EF, FE, EV, EFi;
  Eigen::VectorXi innerEdges;

  // first testing case
//  igl::readOBJ("/home/kacper/Projects/Directional/tutorial/shared/ellipsoid_3-param-full-seamless.obj", V, TC, N, F, FTC, FN);
//  hedra::polygonal_edge_topology(VectorXi::Constant(F.rows(),3), F, EV, FE, EF, EFi, FEs, innerEdges);
//  Eigen::RowVector3d spans = V.colwise().maxCoeff() - V.colwise().minCoeff();
//  hedra::copyleft::cgal::generate_mesh(4, V, F, EV, FE, EF, innerEdges, TC, FTC, newV, newD, newF);
//  hedra::polygonal_write_OFF(std::string("ellipsoid_3-param-full-seamless.off"), newV, newD, newF);
  Eigen::MatrixXd fineV;
  Eigen::VectorXi fineD;
  Eigen::MatrixXi fineF;

  igl::readOFF("/home/kacper/tt.off", V, F);
  Eigen::VectorXi D = Eigen::VectorXi::Constant(F.rows(), 3);
  hedra::dual_mesh(V, D, F, hedra::LINEAR_SUBDIVISION, fineV, fineD, fineF);
  hedra::polygonal_write_OFF(std::string("/home/kacper/hh.off"), fineV, fineD, fineF);

//  igl::readOBJ("/home/kacper/Projects/PHex2/data/single_flap_grid_line_on_boundary.obj", V, TC, N, F, FTC, FN);
//  hedra::polygonal_edge_topology(VectorXi::Constant(F.rows(),3), F, EV, FE, EF, EFi, FEs, innerEdges);
//  spans = V.colwise().maxCoeff() - V.colwise().minCoeff();
//  hedra::copyleft::cgal::generate_mesh(4, V, F, EV, FE, EF, innerEdges, TC, FTC, newV, newD, newF);
//  hedra::polygonal_write_OFF(std::string("single_flap_grid_line_on_boundary.off"), newV, newD, newF);

  //second testing case
//  igl::readOBJ("/home/kacper/Projects/PHex2/data/single_flap_grid_line_off_boundary.obj", V, TC, N, F, FTC, FN);
//  hedra::polygonal_edge_topology(VectorXi::Constant(F.rows(),3), F, EV, FE, EF, EFi, FEs, innerEdges);
//  spans = V.colwise().maxCoeff() - V.colwise().minCoeff();
//  hedra::copyleft::cgal::generate_mesh(4, V, F, EV, FE, EF, innerEdges, TC, FTC, newV, newD, newF);
//  hedra::polygonal_write_OFF(std::string("single_flap_grid_line_off_boundary.off"), newV, newD, newF);
//
//  //third testing case
//  igl::readOBJ("/home/kacper/Projects/PHex2/data/single_flap_grid_virtual_vertex_on_boundary.obj", V, TC, N, F, FTC, FN);
//  hedra::polygonal_edge_topology(VectorXi::Constant(F.rows(),3), F, EV, FE, EF, EFi, FEs, innerEdges);
//  spans = V.colwise().maxCoeff() - V.colwise().minCoeff();
//  hedra::copyleft::cgal::generate_mesh(4, V, F, EV, FE, EF, innerEdges, TC, FTC, newV, newD, newF);
//  hedra::polygonal_write_OFF(std::string("single_flap_grid_virtual_vertex_on_boundary.off"), newV, newD, newF);
//
//  //fourth testing case
//  igl::readOBJ("/home/kacper/Projects/PHex2/data/single_flap_grid_vertex_on_boundary.obj", V, TC, N, F, FTC, FN);
//  hedra::polygonal_edge_topology(VectorXi::Constant(F.rows(),3), F, EV, FE, EF, EFi, FEs, innerEdges);
//  spans = V.colwise().maxCoeff() - V.colwise().minCoeff();
//  hedra::copyleft::cgal::generate_mesh(4, V, F, EV, FE, EF, innerEdges, TC, FTC, newV, newD, newF);
//  hedra::polygonal_write_OFF(std::string("single_flap_grid_vertex_on_boundary.off"), newV, newD, newF);
}
