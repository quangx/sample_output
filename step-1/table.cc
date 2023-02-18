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
// And this for the declarations of the `std::sqrt` and `std::fabs` functions:
#include <cmath>
#include <random>
// The final step in importing deal.II is this: All deal.II functions and
// classes are in a namespace <code>dealii</code>, to make sure they don't
// clash with symbols from other libraries you may want to use in conjunction
// with deal.II. One could use these functions and classes by prefixing every
// use of these names by <code>dealii::</code>, but that would quickly become
// cumbersome and annoying. Rather, we simply import the entire deal.II
// namespace for general use:
using namespace dealii;

// @sect3{Creating the first mesh}

// In the following, first function, we simply use the unit square as domain
// and produce a globally refined grid from it.
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
  return t;
}


void to_vtk(Table<3,double> t,std::string name){
  TableIndices<3> dim=t.size();
  int n1=dim[0];
  int n2=dim[1];
  int n3=dim[2];
  std::ofstream file;
  file.open(name);
  file<<"<VTKFile type=\"ImageData\" "
  "byte_order=\"LittleEndian\"> \n \t"
  "<ImageData WholeExtent=\"0 " +std::to_string(dim[0]-1)+
  " 0 "+std::to_string(dim[1]-1)+" 0 "+std::to_string(dim[2]-1)
  +"\" Origin=\"0 0 0\" Spacing=\".1 .1 .1\">"
  "\n\t <Piece Extent=\"0 "+std::to_string(dim[0]-1)+
  " 0 "+std::to_string(dim[1]-1)+" 0 "+std::to_string(dim[2]-1)+"\">"

  "\n\t<PointData scalars=\"T\">"
  "\n\t <DataArray Name=\"T\" type=\"Float32\" format=\"ascii\">";

  for(int z=0;z<n3;z++){

    for(int y=0;y<n2;y++){
      for (int x=0;x<n1;x++){
        file<<std::to_string(t[x][y][z])+" ";
      }
    }
  }
  file<<"</DataArray> \n </PointData> \n "
  "<CellData> </CellData> \n"
  "</Piece> \n </ImageData> \n </VTKFile>";


}
int main(){
  srand(time(NULL));
  file_generator(100,100,100,"data.txt");
  std::ifstream myFile("data.txt");
  Table<3,double> t=table_generator(90,100,100,myFile);
  to_vtk(t,"image.vti");
  std::cout<<t[1][2][3];

}