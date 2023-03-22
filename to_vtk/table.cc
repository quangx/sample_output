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


struct StructuredData{
  private: 
     Table<3,double> priorities; 
  public:
    Table<4,double> data;
    Point<3,double> min, max;
    std::array<double,3> spacing;  

    StructuredData(Table<4,double> &d, const Point<3,double> min,const Point<3,double> max,std::array<double,3> spacing){
        data=d;
        this->min=min;
        this->max=max;
        this->spacing=spacing;
    }
    std::array<unsigned int,3> location_to_index(const Point<3,double> &p){
      std::array<unsigned int,3> index;
      for(int i=0;i<3;++i){
        unsigned int temp_index=std::floor((p(i)-min[i])/spacing[i]);
        if(((1.0*temp_index*spacing[i]+min[i])+(1.0*(temp_index+1)*spacing[i]+min[i]))/2 <= p(i)){
          ++temp_index;
        }
        index[i]=temp_index;
      }
      
      return index;
    }

    Point<3,double> index_to_location(const
    std::array<unsigned int,3> &idx){
      Point<3,double> location;
      for(int i=0;i<3;++i){
        location(i)=min[i]+idx[i]*spacing[i];
      }
      return location;
    }
    void set_values(const std::array<unsigned int,3> &idx,             //unsure if this is best way to check closest
    double distance, std::vector<double> &values){
      if(1./(1e-20+distance)> priorities[idx[0]][idx[1]][idx[2]]){
        for(int i=0;i<values.size();++i){
          data[idx[0]][idx[1]][idx[2]][i]=values[i];
        }
        priorities[idx[0]][idx[1]][idx[2]]=1./(1e-20+distance);
      }
    }
    std::array<unsigned int,3> approximate_extent(const Point<3,double> &position,
    double radius){
      //implement later
      std::array<unsigned int,3> result;
      for(unsigned int d=0;d<3;++d){
        result[d]=3;
      }
      return result;
    }
    void splat(const Point<3,double> &p,std::vector<double> &values,const double radius){   //todo: implement
      const std::array<unsigned int,3> idx=location_to_index(p);
      std::array<unsigned int,3> extent=approximate_extent(p,radius);

    }
    


};
enum class DataInterpretation{
  component_is_scalar,component_is_vector
};
std::string to_string(DataInterpretation& d){
  switch(d){
    case DataInterpretation::component_is_scalar:
      return "Scalar";
    case DataInterpretation::component_is_vector:
      return "Vector";
    

  }
}
std::ostream& operator<<(std::ostream& stm, DataInterpretation& d){
  switch(d){
    case DataInterpretation::component_is_scalar:
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
void vector_file_generator(int size1,int size2,int size3,std::string name,int num_types){
  double lower_bound=1;
  double upper_bound=100;
  std::ofstream file;
  file.open(name);
  for(int i=0;i<num_types;++i){   // depends on how many components you want to pass in
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
  int num_types=types.size();
  Table<4,double> tab(dims[0],dims[1],dims[2],num_types);
  for(int i=0;i<num_types;++i){
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
//the mesh is implied by the uniform spacing (b-a)/(n-1). pts_direction tells us what n1,n2,n3 are for the smaller grid 
Table<4,double> sampler(Table<4,double> &t,const Point<3,double> min, Point<3,double> max,std::vector<int> pts_direction){     
  TableIndices<4> table_dim=t.size();
  int n1=table_dim[0];
  int n2=table_dim[1];
  int n3=table_dim[2];
  Table<4,double> sampled_table(n1,n2,n3,table_dim[3]);
  std::vector<double> spacing;
  for(int i=0;i<3;++i){
    spacing.push_back((double)(max[i]-min[i])/(double)(table_dim[i]));
  }

}
//returns [(x1 y1 z1)  (x2 y2 z2)] which is the rectangular prism the given point is contained in

int main(){
  srand(time(NULL));
  std::vector<int> dims{19,17,15};
  std::vector<DataInterpretation> types{DataInterpretation::component_is_vector,DataInterpretation::component_is_vector
  ,DataInterpretation::component_is_vector,DataInterpretation::component_is_scalar,DataInterpretation::component_is_vector,
  DataInterpretation::component_is_vector,DataInterpretation::component_is_vector
  };
  vector_file_generator(19,17,15,"find_location.txt",types.size());
  std::ifstream myFile("find_location.txt");

  Table<4,double> t=table_generator(dims,myFile,types);

  // std::vector<std::string> names{"x_velocity","y_velocity","z_velocity","viscosity","mu1","mu2","mu3"};
  // to_vtk(t,Point<3,int>(0,0,0),Point<3,int>(9,8,7),"vector_data_large.vti",types,names);
  Triangulation<3> tria;
  Point<3,double> p1(0,0,0);
  Point<3,double> p2(9,8,7);
  std::vector<int> num_dir{19,17,15};
  std::array<double,3> spacing{.5,.5,.5};
  StructuredData s(t,p1,p2,spacing);
  // GridGenerator::hyper_rectangle(tria,p1,p2);
  // tria.refine_global(2);
  // int i=0;
  // for(auto &cell:tria.active_cell_iterators()){
  //   if(i%3==0){
  //     cell->set_refine_flag();
  //   }
  //   ++i;
  // }
  // tria.execute_coarsening_and_refinement();
  
  // std::ofstream out("grid.vtu");
  // GridOut gridout;
  // gridout.write_vtu(tria,out);

  //test location_to_index
  Point<3,double> p_test(.24,.5,.25);
  std::array<unsigned int,3> indices=s.location_to_index(p_test);
  for(int i=0;i<3;++i){
    std::cout<<std::to_string(indices[i])+" ";
  }
  //test index_to_location
  std::array<unsigned int,3> index_test{2,3,4};
  Point<3,double> itl=s.index_to_location(index_test);
  for(int i=0;i<3;++i){
    std::cout<<std::to_string(itl(i))+" ";
  }
  

}