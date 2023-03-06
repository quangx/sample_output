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

enum class DataInterpretation{
  compoenent_is_scalar,component_is_vector
};
std::string to_string(DataInterpretation& d){
  switch(d){
    case DataInterpretation::compoenent_is_scalar:
      return "Scalar";
    case DataInterpretation::component_is_vector:
      return "Vector";
    

  }
}
std::ostream& operator<<(std::ostream& stm, DataInterpretation& d){
  switch(d){
    case DataInterpretation::compoenent_is_scalar:
      return stm<<"Scalar";
    case DataInterpretation::component_is_vector:
      return stm<<"Vector";
    

  }
}
// generates data file to read in to table
// void file_generator(int size1, int size2, int size3,std::string name){
//   double lower_bound=1;
//   double upper_bound=100;
  

//   std::ofstream file;
//   file.open(name);
//   for(int z=0;z<size3;++z){

//     for(int y=0;y<size2;++y){
//       for(int x=0;x<size1;++x){
//           file<<std::to_string((double)rand()/RAND_MAX*(upper_bound-lower_bound)+lower_bound)+" ";
//       }
//       file<<"\n";
//     }
//     file<<"\n";
//   }
// }


//testing indexing
void file_generator(int size1, int size2, int size3,std::string name){
  double lower_bound=1;
  double upper_bound=100;
  

  std::ofstream file;
  file.open(name);
  
    for(int z=0;z<size3;++z){

      for(int y=0;y<size2;++y){
        for(int x=0;x<size1;++x){
            int t=5*y+5*z;
            file<<std::to_string(t)+" ";
        }
        file<<"\n";
      }
      file<<"\n";
    }
  
}
void vector_file_generator(int size1,int size2,int size3,std::string name){
  double lower_bound=1;
  double upper_bound=100;
  std::ofstream file;
  file.open(name);
  for(int i=0;i<4;++i){
    for(int z=0;z<size3;++z){
      for(int y=0;y<size1;++y){
        for(int x=0;x<size1;++x){

            file<<std::to_string((double)rand()/RAND_MAX*(upper_bound-lower_bound)+lower_bound)+" ";
          
          file<<"\n";
        }
        file<<"\n";
      }
      file<<"\n";
    }
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
Table<4,double> table_generator(std::vector<int> dims,std::ifstream& myFile,std::vector<DataInterpretation> types){   
  Table<4,double> tab(dims[0],dims[1],dims[2],4);
  for(int i=0;i<4;++i){
    for(int z=0;z<dims[2];++z){
      for(int y=0;y<dims[1];++y){
        for(int x=0;x<dims[0];++x){
          double d;
          myFile>>d;
          tab[x][y][z][i]=d;
        }
      }
    }
  }
  return tab;
}
void to_vtk(Table<4,double> &t,const Point<3,int> min,
const Point<3,int> max,const std::string filename,
std::vector<DataInterpretation> component_type,const std::vector<std::string> names){
  int n1;
  int n2;
  int n3;
  std::vector<double> spacing;
  TableIndices<4> table_dim=t.size();
  
  n1=table_dim[0];
  n2=table_dim[1];
  n3=table_dim[2];
  for(int i=0;i<3;++i){
    spacing.push_back((double)(max[i]-min[i])/(double)(table_dim[i]-1));
  }
  
  // else if(dim==3){
  //   TableIndices<3> table_dim=t.size();
  //   n1=table_dim[0];
  //   n2=table_dim[1];
  //   n3=table_dim[2];
    
  //   for(int i=0;i<3;++i){
  //     spacing.push_back((double)(max[i]-min[i])/(double)(table_dim[i]-1));
  //   }
  // }
  
  std::ofstream file;
  file.open(filename);
  file<<"<VTKFile type=\"ImageData\" "
  "byte_order=\"LittleEndian\"> \n \t"
  "<ImageData WholeExtent=\""+std::to_string(min[0])+" "+std::to_string(max[0])+""
  " "+std::to_string(min[1])+" "+std::to_string(max[1])+" "+std::to_string(min[2])+" "+std::to_string(max[2])+""
  "\" Origin=\"0 0 0\" Spacing=\""+std::to_string(spacing[0])+" "+std::to_string(spacing[1])+" "+std::to_string(spacing[2])+"\">"
  "\n\t <Piece Extent=\""+std::to_string(min[0])+" "+std::to_string(max[0])+" "
  ""+std::to_string(min[1])+" "+std::to_string(max[1])+" "+std::to_string(min[2])+" "+std::to_string(max[2])+"\">"
  "\n\t<PointData Scalars=\"T\">";
  for(int i=0;i<component_type.size();++i){
    file<<"\n\t <DataArray type=\"Float32\" Name=\""+names[i]+"\" ";
    if(component_type[i]==DataInterpretation::component_is_vector){
      file<<" NumberOfComponents=\"3\" ";

    }
    file<<"format=\"ascii\">\n";
    if(component_type[i]==DataInterpretation::component_is_vector){
      for(int z=0;z<n3;++z){
        for(int y=0;y<n2;++y){
          for(int x=0;x<n1;++x){
            for(int j=0;j<3;++j){
              file<<std::to_string(t[x][y][z][i+j])+" ";
            }
          }
        }
      }
      i=i+2;
    }
    else{
      for(int z=0;z<n3;++z){
        for(int y=0;y<n2;++y){
          for(int x=0;x<n1;++x){
            file<<std::to_string(t[x][y][z][i])+" ";
          }
        }
      }
    }
   
   
  
    file<<"</DataArray> \n ";
  }

  file<<"</PointData> \n "
  "<CellData> </CellData> \n"
  "</Piece> \n </ImageData> \n </VTKFile>";


}

void to_vtk_cell(Table<4,double> &t,const Point<3,int> min,
const Point<3,int> max,const std::string filename,
std::vector<DataInterpretation> component_type,const std::vector<std::string> names){
  int n1;
  int n2;
  int n3;
  std::vector<double> spacing;
  TableIndices<4> table_dim=t.size();
  
  n1=table_dim[0];
  n2=table_dim[1];
  n3=table_dim[2];
  for(int i=0;i<3;++i){
    spacing.push_back((double)(max[i]-min[i])/(double)(table_dim[i]));
  }
  
 
  std::ofstream file;
  file.open(filename);
  file<<"<VTKFile type=\"ImageData\" "
  "byte_order=\"LittleEndian\"> \n \t"
  "<ImageData WholeExtent=\""+std::to_string(min[0])+" "+std::to_string(max[0])+""
  " "+std::to_string(min[1])+" "+std::to_string(max[1])+" "+std::to_string(min[2])+" "+std::to_string(max[2])+""
  "\" Origin=\"0 0 0\" Spacing=\""+std::to_string(spacing[0])+" "+std::to_string(spacing[1])+" "+std::to_string(spacing[2])+"\">"
  "\n\t <Piece Extent=\""+std::to_string(min[0])+" "+std::to_string(max[0])+" "
  ""+std::to_string(min[1])+" "+std::to_string(max[1])+" "+std::to_string(min[2])+" "+std::to_string(max[2])+"\">"
  "\n\t<PointData Scalars=\"T\">";
  

  file<<"</PointData> \n "
  "<CellData>";
  for(int i=0;i<component_type.size();++i){
    file<<"\n\t <DataArray type=\"Float32\" Name=\""+names[i]+"\" ";
    if(component_type[i]==DataInterpretation::component_is_vector){
      file<<" NumberOfComponents=\"3\" ";

    }
    file<<"format=\"ascii\">\n";
    if(component_type[i]==DataInterpretation::component_is_vector){
      for(int z=0;z<n3;++z){
        for(int y=0;y<n2;++y){
          for(int x=0;x<n1;++x){
            for(int j=0;j<3;++j){
              file<<std::to_string(t[x][y][z][i+j])+" ";
            }
          }
        }
      }
      i=i+2;
    }
    else{
      for(int z=0;z<n3;++z){
        for(int y=0;y<n2;++y){
          for(int x=0;x<n1;++x){
            file<<std::to_string(t[x][y][z][i])+" ";
          }
        }
      }
    }
   
   
  
    file<<"</DataArray> \n ";
  }

  file<<"</CellData> \n"
  "</Piece> \n </ImageData> \n </VTKFile>";


}
int main(){
  srand(time(NULL));
  // file_generator(10,9,8,"data.txt"); //make sure dimensino matches t
  // std::ifstream myFile("data.txt");
  // Table<3,double> t=table_generator(10,9,8,myFile);
  vector_file_generator(10,9,8,"vector_data.txt");
  std::ifstream myFile("vector_data.txt");
  std::vector<int> dims{10,9,8};
  std::vector<DataInterpretation> types{DataInterpretation::component_is_vector,DataInterpretation::component_is_vector
  ,DataInterpretation::component_is_vector,DataInterpretation::compoenent_is_scalar
  };
  Table<4,double> t=table_generator(dims,myFile,types);

  std::vector<std::string> names{"velocity","y_velocity","z_velocity","viscosity"};
  to_vtk(t,Point<3,int>(0,0,0),Point<3,int>(9,8,7),"vector_test.vti",types,names);

}