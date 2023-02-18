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
#include <deal.II/base/function_lib.h>
#include <deal.II/numerics/vector_tools.h>

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
void
check()
{
  Triangulation<dim> tria;
  GridGenerator::hyper_cube(tria, 0., 1.);
  tria.refine_global(2);

  Vector<double> cell_data(tria.n_active_cells());
  for (unsigned int i = 0; i < tria.n_active_cells(); ++i)
    cell_data(i) = i * 1.0;

  FE_Q<dim>       fe(1);
  DoFHandler<dim> dof_handler(tria);
  dof_handler.distribute_dofs(fe);

  Vector<double> x(dof_handler.n_dofs());
  VectorTools::interpolate(dof_handler,
			   Functions::PillowFunction<dim>(),
			   x);

  DataOut<dim> data_out;
  data_out.attach_dof_handler(dof_handler);
  data_out.add_data_vector(cell_data,
                           "cell_data",
                           DataOut<dim>::type_cell_data);
  data_out.add_data_vector(x, "solution");

  data_out.build_patches();

  data_out.write_deal_II_intermediate_in_parallel(
    "test.pd2",
    MPI_COMM_WORLD,
		    DataOutBase::CompressionLevel::best_compression);

  std::cout << "OK" << std::endl;
}

int
main(int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);

  check<3>();

  return 0;
}
