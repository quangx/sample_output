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
    std::array<unsigned int,3> num_values;


    StructuredData(const Point<3,int> &min,const Point<3,int> &max,const std::array<double,3> &spacing,const int num_components){

        this->min=min;
        this->max=max;
        this->spacing=spacing;
        for(int i=0;i<3;++i){
          num_values[i]=((max[i]-min[i])/spacing[i])+1;
          // num_values[i]=max[i]+1;
        }
        TableIndices<4> t_ind(num_values[0],num_values[1],num_values[2],num_components);
        data.reinit(t_ind);
        priorities.reinit(num_values[0],num_values[1],num_values[2]);

        // priorities(num_values[0],num_values[1],num_values[2]);
        // TableIndices<3> p_size=priorities.size();
        // for(int i=0;i<3;++i){
        //   std::cout<<std::to_string(p_size[i])+" ";
        // }
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
      for(unsigned int iz=0;iz<2*extent[2];++iz){
        for(unsigned int iy=0;iy<2*extent[1];++iy){
          for(unsigned int ix=0;ix<2*extent[0];++ix){
            std::array<unsigned int,3> current_index;
            current_index[0] = std::max(0l,static_cast<long>(idx[0] + ix) - static_cast<long>(extent[0]));
            current_index[1] = std::max(0l,static_cast<long>(idx[1] + iy) - static_cast<long>(extent[1]));
            current_index[2] = std::max(0l,static_cast<long>(idx[2] + iz) - static_cast<long>(extent[2]));
            if(current_index[0]>=num_values[0]||
            current_index[1]>=num_values[1] ||
            current_index[2]>=num_values[2]){
              continue;
            }
            
            const double distance=index_to_location(current_index).distance(p);
            set_values(current_index,distance,values);
          
        }
      }

    }
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
      for(int y=0;y<size2;++y){
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
  
  
  
  std::ofstream file;
  file.open(filename);
  file<<"<VTKFile type=\"ImageData\" "
  "byte_order=\"LittleEndian\" header_type=\"UInt64\"> \n \t"
  "<ImageData WholeExtent=\""+std::to_string(0)+" "+std::to_string(n1-1)+""
  " "+std::to_string(0)+" "+std::to_string(n2-1)+" "+std::to_string(0)+" "+std::to_string(n3-1)+""
  "\" Origin=\"0 0 0\" Spacing=\""+std::to_string(spacing[0])+" "+std::to_string(spacing[1])+" "+std::to_string(spacing[2])+"\" Direction=\"1 0 0 0 1 0 0 0 1\">"
  "\n\t <Piece Extent=\""+std::to_string(0)+" "+std::to_string(n1-1)+" "
  ""+std::to_string(0)+" "+std::to_string(n2-1)+" "+std::to_string(0)+" "+std::to_string(n3-1)+"\">"
  "\n\t<PointData Scalars=\"T\">";
  for(unsigned int i=0;i<component_type.size();++i){
    file<<"\n\t <DataArray type=\"Float32\" Name=\""+names[i]+"\" ";
    if(component_type[i]==DataInterpretation::component_is_vector){
      file<<" NumberOfComponents=\"3\" ";

    }
    file<<"format=\"ascii\">\n";
    if(component_type[i]==DataInterpretation::component_is_vector){
      for(int z=0;z<n3;++z){
        file<<"\t";
        for(int y=0;y<n2;++y){
          for(int x=0;x<n1;++x){
            for(int j=0;j<3;++j){
              file<<std::to_string(t[x][y][z][i+j])+" ";
            }
            file<<"\n\t";
          }
          file<<"\n\t";
        }
        file<<"\n";
      }
      i=i+2;
    }
    else{
      for(int z=0;z<n3;++z){
        file<<"\t";
        for(int y=0;y<n2;++y){
          for(int x=0;x<n1;++x){
            file<<std::to_string(t[x][y][z][i])+" ";
          }
          file<<"\n\t";
        }
        file<<"\n";
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
  for(unsigned int i=0;i<component_type.size();++i){
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
Table<4,double> file_to_table(std::ifstream & myFile){
  
}

int main(){
  srand(time(NULL));
  Point<3,int> p1(0,0,0);
  Point<3,int> p2(9,8,7);
  std::vector<int> dims{19,17,15};
  
  std::vector<DataInterpretation> components{DataInterpretation::component_is_scalar,DataInterpretation::component_is_vector,
  DataInterpretation::component_is_vector,DataInterpretation::component_is_vector};
  std::vector<std::string> names{"mu"};
  // vector_file_generator(19,17,15,"retest.txt",components.size());
  std::ifstream myFile("retest.txt");
  // to_vtk(T,p1,p2,"retest.vti",components,names);
  std::vector<Point<3,double>> pts;
  std::vector<double> data;
  std::ifstream corner;
  corner.open("corners.txt");
  while(!corner.eof()){
    double a;
    corner>> a;
    double b;
    corner>>b;
    double c;
    corner>> c;
    double d;
    corner>>d;
    Point<3,double> temp(a,b,c);
    pts.push_back(temp);
    data.push_back(d);
  }
  StructuredData s(p1,p2,{.5,.5,.5},1);
  for(int i=0;i<pts.size();++i){
    std::vector<double> temp;
    temp.push_back(data[i]);
    s.splat(pts[i],temp,2);
  }
  Table<4,double> T=s.data;
  
 
  // Table<4,double> T=s.data;
  to_vtk(T,p1,p2,"sample2.vti",{DataInterpretation::component_is_scalar},{"mu"});

  

  // for(unsigned int z=0;z<n3;++z){
  //   for(unsigned int y=0;y<n2;++y){
  //     for(unsigned int x=0;x<n1;++x){
  //       std::vector<double> data;
  //       std::array<unsigned int,3> idx{x,y,z};
  //       for(unsigned int i=0;i<n4;++i){
  //         data.push_back(t_small[x][y][z][i]);
  //       }
  //       Point<3,double> p_temp=small.index_to_location(idx);

  //       large.splat(p_temp,data,2);
  //     }
  //   }
  // }



  
  
  
  // std::ofstream out("grid.vtu");
  // GridOut gridout;
  // gridout.write_vtu(tria,out);

  //test location_to_index
  // Point<3,double> p_test(.24,.5,.25);
  // std::array<unsigned int,3> indices=s.location_to_index(p_test);
  // for(int i=0;i<3;++i){
  //   std::cout<<std::to_string(indices[i])+" ";
  // }
  // //test index_to_location
  // std::array<unsigned int,3> index_test{2,3,4};
  // Point<3,double> itl=s.index_to_location(index_test);
  // for(int i=0;i<3;++i){
  //   std::cout<<std::to_string(itl(i))+" ";
  // }

  

}