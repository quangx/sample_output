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

using namespace dealii;


template <int dim>
class MyReader: public DataOutReader<dim,dim>
{
public:
  void inspect()
  {
    for (const auto & patch: this->get_patches())
      {
	std::cout << "vertex: " << patch.vertices[0]
		  << " values: ";

	for (unsigned int i=0;i<patch.data.n_cols();++i)
	  std::cout << ' ' << patch.data(0,i);
	std::cout << std::endl;
      }
  }

};

template <int dim>
void
check()
{
  // Read the data back in and dump it into the deallog:
  std::ifstream in("test.pd2");
  Assert(in, dealii::ExcIO());
  MyReader<dim> reader;
  reader.read_whole_parallel_file(in);
  reader.inspect();

  std::cout << "OK" << std::endl;
}

int
main(int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);

  check<3>();

  return 0;
}
