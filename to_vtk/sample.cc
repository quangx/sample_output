// ---------------------------------------------------------------------
//
// Copyright (C) 2022 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE.md at
// the top level directory of deal.II.
//
// ---------------------------------------------------------------------

// Test DataOut::write_deal_II_intermediate_in_parallel() and
// DataOutReader::read_whole_parallel_file() with compression

#include <deal.II/base/mpi.h>

#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>

#include <deal.II/grid/filtered_iterator.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_iterator.h>

#include <deal.II/lac/vector.h>

#include <deal.II/numerics/data_out.h>
#include "structured.h"

using namespace dealii;


template <int dim>
class MyReader: public DataOutReader<dim,dim>
{
public:
  StructuredData inspect(const Point<3,int> &min,const Point<3,int> &max,const std::array<unsigned int,3> &num_pts,const int &num_components)
  {

    StructuredData structured_data(min,max,num_pts,num_components);
    for (const auto & patch: this->get_patches()){
      for(unsigned int k=0;k<8;++k){ //8 vertices in 3d
        Point<3,double> vertex=patch.vertices[k];
        std::vector<double> data;
        for(unsigned int i=0;i<patch.data.n_rows();++i){
          data.push_back(patch.data(i,k));
        }
        structured_data.splat(vertex,data,3);

      }
    }
    return structured_data;

      
  }

  // StructuredData structured_data;
};

void
sample(const std::string &myFile,const std::string &outputName,
Point<3,int> &p1,Point<3,int> &p2,std::array<unsigned int,3> pts_dir,
unsigned int num_components,std::vector<DataInterpretation> datatypes,
std::vector<std::string> names)
{
  // Read the data back in and dump it into the deallog:
  std::ifstream in(myFile);
  Assert(in, dealii::ExcIO());
  MyReader<3> reader;
  reader.read_whole_parallel_file(in);
  
  StructuredData s=reader.inspect(p1,p2,pts_dir,num_components);
  Table<4,double> T=s.data;
  s.to_vtk(T,p1,p2,outputName,datatypes,names);


  std::cout << "OK" << std::endl;
}

int
main(int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
  Point<3,int> p1(0,0,0);
  Point<3,int> p2(9,8,7);
  std::array<unsigned int,3> pts_dir{10,9,8};
  const std::string infile="test.pd2";
  const std::string outfile="sample6.vti";
  unsigned  int num_components=1;
  const std::vector<DataInterpretation> d{DataInterpretation::component_is_scalar};
  const std::vector<std::string> names{"mu"};
  sample(infile,outfile,p1,p2,pts_dir,num_components,d,names);
  
  return 0;
}
