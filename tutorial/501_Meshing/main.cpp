#include <algorithm>
#include <math.h>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/readOBJ.h>
#include <hedra/polygonal_edge_topology.h>
#include <hedra/polygonal_write_OFF.h>

#include <hedra/triangulate_mesh.h>
#include <hedra/point_spheres.h>
#include <hedra/polygonal_edge_lines.h>

#include <hedra/copyleft/cgal/generate_mesh.h>

Eigen::MatrixXd V, TV, newV;
Eigen::MatrixXi F, newF, FTC;
Eigen::VectorXi newD;

Eigen::MatrixXi TMesh;
Eigen::VectorXi TF;

bool key_down(igl::opengl::glfw::Viewer& viewer, unsigned char key, int modifiers)
{
  return true;
}

int main(int argc, char *argv[])
{
  using namespace Eigen;
  using namespace std;
  
  cout<<"0: show original mesh with parameterization"<<endl<<
  "1: show UV coordinates"<<endl<<
  "2: show arrangement with tiling"<<endl<<
  "3: Show final mesh"<<endl;

  Eigen::MatrixXd TC;
  Eigen::MatrixXd N;
  Eigen::MatrixXd FN, FEs;
  Eigen::MatrixXi EF, FE, EV, EFi;
  Eigen::VectorXi innerEdges;
  // the obj file needs to have texture coordinates
  igl::readOBJ(TUTORIAL_SHARED_PATH "/horsers-param-full-seamless.obj", V, TC, N, F, FTC, FN);
  
  //OBJ values are skewed to fit PNG as follows
  //OutputCornerValues=[CornerValues(:,1)/3 CornerValues(:,2)/sqrt(3)];
  //TC.col(0).array()*=3;
  //TC.col(1).array()*=SQRT3;
  
  hedra::polygonal_edge_topology(VectorXi::Constant(F.rows(),3), F, EV, FE, EF, EFi, FEs, innerEdges);
  std::cout << V.rows() << " " << V.cols() << std::endl;
  Eigen::RowVector3d spans = V.colwise().maxCoeff() - V.colwise().minCoeff(); // ? diagonal
  hedra::copyleft::cgal::generate_mesh(4, V, F, EV, FE, EF, EFi, innerEdges, TC, FTC, newV, newD, newF); // ? re-meshing?
  hedra::polygonal_write_OFF(std::string(TUTORIAL_SHARED_PATH "/horsers-quads.off"),newV, newD, newF);

  Eigen::MatrixXi newEV;
  Eigen::MatrixXi newFE;
  Eigen::MatrixXi newEF;
  Eigen::MatrixXi newEFi;
  Eigen::MatrixXd newFEs;
  Eigen::VectorXi newInnerEdges;
  hedra::polygonal_edge_topology(newD, newF, newEV, newFE, newEF, newEFi, newFEs, newInnerEdges);

  igl::opengl::glfw::Viewer viewer;
  viewer.callback_key_down = &key_down;

 //set mesh color with the marked faces: 3954 and 5584.
  Eigen::MatrixXd MColors = hedra::default_mesh_color().replicate(F.rows(), 1);
  MColors.row(3954) = Eigen::RowVector3d(0, 0, 1);
  MColors.row(5584) = Eigen::RowVector3d(0, 0, 1);
  viewer.data_list[0].clear();
  viewer.data_list[0].set_mesh(V, F);
  viewer.data_list[0].set_colors(MColors);
  viewer.data_list[0].set_face_based(true);
  viewer.data_list[0].show_lines=false;

  //edges mesh
  Eigen::MatrixXd VEdges;
  Eigen::MatrixXi TEdges;
  Eigen::MatrixXd CEdges;
  viewer.append_mesh();
  hedra::polygonal_edge_lines(newV, newF, newEV, 0.0001, 10, VEdges, TEdges, CEdges);
  viewer.data_list[1].set_mesh(VEdges, TEdges);
  viewer.data_list[1].set_colors(CEdges);
  viewer.data_list[1].show_lines=false;

  //edges mesh
  Eigen::MatrixXd VEdgesO;
  Eigen::MatrixXi TEdgesO;
  Eigen::MatrixXd CEdgesO;
  viewer.append_mesh();
  hedra::polygonal_edge_lines(V, F, EV, 0.0001, 10, VEdgesO, TEdgesO, CEdgesO, Eigen::RowVector3d(1,0,0));
  viewer.data_list[2].set_mesh(VEdgesO, TEdgesO);
  viewer.data_list[2].set_colors(CEdgesO);
  viewer.data_list[2].show_lines=false;


  //spheres mesh
  Eigen::MatrixXd VSpheres;
  Eigen::MatrixXi TSpheres;
  Eigen::MatrixXd CSpheres;
  Eigen::MatrixXd sphereColors;

  sphereColors.resize(V.rows(),3);
  for (int i=0;i< V.rows();i++)
    sphereColors.row(i) << 0, 1, 0;

  viewer.append_mesh();
  hedra::point_spheres(V, 0.0002, sphereColors, 10, VSpheres, TSpheres, CSpheres);
  viewer.data_list[3].set_mesh(VSpheres, TSpheres);
  viewer.data_list[3].set_colors(CSpheres);
  viewer.data_list[3].show_lines=true;
  viewer.data_list[3].show_faces=true;

  viewer.selected_data_index=0;

  viewer.launch();
}
