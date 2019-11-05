#include <igl/opengl/glfw/Viewer.h>
#include <igl/readOBJ.h>
#include <hedra/polygonal_edge_topology.h>
#include <hedra/copyleft/cgal/generate_mesh.h>
#include <hedra/polygonal_write_OFF.h>

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
  
  cout<<"0: show original mesh with parameterization"<<endl<<
  "1: show UV coordinates"<<endl<<
  "2: show arrangement with tiling"<<endl<<
  "3: Show final mesh"<<endl;

  Eigen::MatrixXd TC;
  Eigen::MatrixXd N;
  Eigen::MatrixXd FN, FEs;
  Eigen::MatrixXi EF, FE, EV, EFi;
  Eigen::VectorXi innerEdges;
  igl::readOBJ("/home/kacper/Projects/PHex2/data/single_flap.obj", V, TC, N, F, FTC, FN);
  //igl::readOBJ("/home/kacper/Projects/Directional/tutorial/shared/horsers-param-full-seamless.obj", V, TC, N, F, FTC, FN);

  //for(int i = 0; i < TC.rows(); i++)
  //  std::cout << TC.row(i) << std::endl;

  std::cout << "FTC # " << FTC.rows() << std::endl;

  //for(int i = 0; i < FTC.rows(); i++)
  //  std::cout << FTC.row(i) << std::endl;
  //exit(1);

  //OBJ values are skewed to fit PNG as follows
  //OutputCornerValues=[CornerValues(:,1)/3 CornerValues(:,2)/sqrt(3)];
  //TC.col(0).array()*=3;
  //TC.col(1).array()*=SQRT3;
  
  hedra::polygonal_edge_topology(VectorXi::Constant(F.rows(),3), F, EV, FE, EF, EFi, FEs, innerEdges);
  
  Eigen::RowVector3d spans=V.colwise().maxCoeff()-V.colwise().minCoeff();
  hedra::copyleft::cgal::generate_mesh(4, V, F, EV, FE, EF, innerEdges, TC, FTC, newV, newD, newF);
  hedra::polygonal_write_OFF(std::string("flap-quads.off"),newV,newD,newF);
  
  igl::opengl::glfw::Viewer viewer;
  viewer.callback_key_down = &key_down;
  
  //viewer.core.background_color<<0.75,0.75,0.75,1.0;
  
  //edges mesh
  viewer.append_mesh();
  
  //control polygon mesh
  viewer.append_mesh();
  
  viewer.selected_data_index=0;
  //update_mesh(viewer);
  //viewer.launch();
}
