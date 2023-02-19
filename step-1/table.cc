/* ---------------------------------------------------------------------
 *
 * Copyright (C) 1999 - 2022 by the deal.II authors
 *
 * This file is part of the deal.II library.
 *
 * The deal.II library is free software; you can use it, redistribute
 * it, and/or modify it under the terms of the GNU Lesser General
 * Public License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 * The full text of the license can be found in the file LICENSE.md at
 * the top level directory of deal.II.
 *
 * ---------------------------------------------------------------------

 */
// #include <vtkActor.h>
// #include <vtkImageData.h>
// #include <vtkImageDataGeometryFilter.h>
// #include <vtkNamedColors.h>
// #include <vtkNew.h>
// #include <vtkPolyDataMapper.h>
// #include <vtkProperty.h>
// #include <vtkRenderWindow.h>
// #include <vtkRenderWindowInteractor.h>
// #include <vtkRenderer.h>
// #include <vtkXMLImageDataReader.h>
// #include <vtkXMLImageDataWriter.h>
// @sect3{Include files}
// The most fundamental class in the library is the Triangulation class, which
// is declared here:
#include <deal.II/grid/tria.h>
// Here are some functions to generate standard grids:
#include <deal.II/grid/grid_generator.h>
// Output of grids in various graphics formats:
#include <deal.II/grid/grid_out.h>

// This is needed for C++ output:
#include <iostream>
#include <fstream>
#include <cmath>
#include <random>

using namespace dealii;

enum class DataInterpretation{T,U,V
};
std::string to_string(DataInterpretation& d){
  switch(d){
    case DataInterpretation::T:
      return "T";
    case DataInterpretation::U:
      return "U";
    case DataInterpretation::V:
      return "V";

  }
}
std::ostream& operator<<(std::ostream& stm, DataInterpretation& d){
  switch(d){
    case DataInterpretation::T:
      return stm<<"T";
    case DataInterpretation::U:
      return stm<<"U";
    case DataInterpretation::V:
      return stm<<"V";

  }
}

void file_generator(int size1, int size2, int size3,std::string name){
  double lower_bound=1;
  double upper_bound=100;
  

  std::ofstream file;
  file.open(name);
  for(int z=0;z<size3;++z){

    for(int y=0;y<size2;++y){
      for(int x=0;x<size1;++x){
          file<<std::to_string((double)rand()/RAND_MAX*(upper_bound-lower_bound)+lower_bound)+" ";
      }
      file<<"\n";
    }
    file<<"\n";
  }
}
Table<3,double> table_generator(const int n1,const int n2, const int n3, std::ifstream& myFile){
  Table<3,double> t(n1,n2,n3);
  for(int z=0;z<n3;++z){
    for(int y=0;y<n2;++y){
      for(int x=0;x<n1;++x){
        double d;
        myFile>>d;
        t[x][y][z]=d;
      }
    }
  }

  // testing indexing
 
  // for(int y=0;y<n2;y++){
  //   for(int x=0;x<n1;x++){
  //     t[x][0][0]=t[x][0][0]+t[0][y][0];
  //   }
  // }
    // for(int i=0;i<n1;++i){
    //   t[i][0][0]=5*t[i][0][0]+5*t[0][i][0];
    // }
   

    
  return t;
}


void to_vtk(Table<3,double> &t, Point<3,int> min,
Point<3,int> max,const std::string filename,
std::vector<DataInterpretation> component_name){
  TableIndices<3> dim=t.size();
  int n1=dim[0];
  int n2=dim[1];
  int n3=dim[2];
  std::ofstream file;
  file.open(filename);
  file<<"<VTKFile type=\"ImageData\" "
  "byte_order=\"LittleEndian\"> \n \t"
  "<ImageData WholeExtent=\""+std::to_string(min[0])+" "+std::to_string(max[0])+""
  " "+std::to_string(min[1])+" "+std::to_string(max[1])+" "+std::to_string(min[2])+" "+std::to_string(max[2])+""
  "\" Origin=\"0 0 0\" Spacing=\"1 1 1\">"
  "\n\t <Piece Extent=\""+std::to_string(min[0])+" "+std::to_string(max[0])+" "
  ""+std::to_string(min[1])+" "+std::to_string(max[1])+" "+std::to_string(min[2])+" "+std::to_string(max[2])+"\">"
  "\n\t<PointData scalars=\"T\">";
  for(int i=0;i<component_name.size();++i){
    file<<"\n\t <DataArray Name=\""+to_string(component_name[i])+"\" type=\"Float32\" format=\"ascii\">\n";

    for(int z=0;z<n3;z++){

      for(int y=0;y<n2;y++){
        for (int x=0;x<n1;x++){
          file<<std::to_string(t[x][y][z])+" ";
        }
        file<<"\n";
      }
      file<< "\n";
    }
  
  file<<"</DataArray> \n ";
  }

  file<<"</PointData> \n "
  "<CellData> </CellData> \n"
  "</Piece> \n </ImageData> \n </VTKFile>";


}
int main(){
  srand(time(NULL));
  file_generator(10,9,8,"data.txt");
  std::ifstream myFile("data.txt");
  Table<3,double> t=table_generator(10,9,8,myFile);
  std::vector<DataInterpretation> v{DataInterpretation::V,DataInterpretation::U};
  to_vtk(t,Point<3,int>(0,0,0),Point<3,int>(9,8,7),"image.vti",v);
  

}