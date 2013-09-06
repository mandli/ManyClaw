#include "data_structures.h"
#include <string>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <fstream>
#include <cmath>

Grid::Grid(int nx, int ny)
{
  // TODO make this a function of the encapsulated states
  // Cell centered values
  num_cells[0] = nx;
  num_cells[1] = ny;
}

State::State(Grid& grid, int num_eqn, int num_aux, int num_ghost, const void* aux_global) :
  num_eqn(num_eqn), num_aux(num_aux), num_ghost(num_ghost), aux_global(aux_global), grid(grid)
{
  const int nx = grid.num_cells[0];
  const int ny = grid.num_cells[1];
  q.resize((nx + num_ghost * 2) * (ny + num_ghost * 2) * num_eqn);
  aux.resize((nx + num_ghost * 2) * (ny + num_ghost * 2) * num_aux);
}

void Solution::write(int frame, std::string output_path)
{
  const char file_prefix[] = "fort.";

  // TODO: Create output directory if not present

  std::ostringstream q_file_path;
  q_file_path << output_path << "/" << file_prefix << "q"
              << std::setfill('0') << std::setw(4) << frame;
  std::cout << "Writing out frame: " << frame << " to " << q_file_path.str() << '\n';

  std::fstream q_file;
  q_file.open(q_file_path.str().c_str(), std::ios::out);

  // Write header for fort.q file
  q_file << std::setiosflags(std::ios::fixed) << std::setprecision(5) 
              << 1 << "                 grid_number\n";
  q_file << std::setiosflags(std::ios::fixed) << std::setprecision(5)
              << 1 << "                 AMR_level\n";
  q_file << std::setiosflags(std::ios::fixed) << std::setprecision(5)
              << grid.num_cells[0] << "                 mx\n";
  q_file << std::setiosflags(std::ios::fixed) << std::setprecision(5)
              << grid.num_cells[1] << "                 my\n";
  q_file << std::setiosflags(std::ios::fixed) << std::setprecision(8)
              << std::setw(18) << grid.lower[0] << "    xlow\n";
  q_file << std::setiosflags(std::ios::fixed) << std::setprecision(8)
              << std::setw(18) << grid.lower[1] << "    ylow\n";
  q_file << std::setiosflags(std::ios::fixed) << std::setprecision(8)
              << std::setw(18) << grid.dx[0] << "    dx\n\n";
  q_file << std::setiosflags(std::ios::fixed) << std::setprecision(8)
              << std::setw(18) << grid.dx[1] << "    dy\n\n";

  // Write out q data
  FieldIndexer fi(grid.num_cells[0], grid.num_cells[1], state.num_ghost,
                  state.num_eqn);
  for (int row = state.num_ghost; row < grid.num_cells[1] + state.num_ghost; ++row)
  {
    for (int col = state.num_ghost; col < grid.num_cells[0] + state.num_ghost; ++col)
    {
      for (int eqn = 0; eqn < state.num_eqn; ++eqn)
      {
        q_file << std::setiosflags(std::ios::fixed) 
               << std::setprecision(8)
               << std::setw(16) 
               << state.q[eqn + fi.idx(row, col)] 
               << " ";
      }     
      q_file << "\n"; 
    }
    q_file << "\n";
  }

  q_file.close();

  // Write out fort.t file
  std::ostringstream t_file_path;
  t_file_path << output_path << "/" << file_prefix << "t"
              << std::setfill('0') << std::setw(4) << frame;

  std::fstream t_file;
  t_file.open(t_file_path.str().c_str(), std::ios::out);

  t_file << std::setiosflags(std::ios::fixed) << std::setprecision(16)
              << std::setw(26) << t << "        time\n";
  t_file << std::setiosflags(std::ios::fixed) << std::setprecision(5) 
              << state.num_eqn << "                 meqn\n";
  t_file << std::setiosflags(std::ios::fixed) << std::setprecision(5) 
              << 1 << "                 ngrids\n";
  t_file << std::setiosflags(std::ios::fixed) << std::setprecision(5) 
              << 0 << "                 maux\n";
  t_file << std::setiosflags(std::ios::fixed) << std::setprecision(5) 
              << 2 << "                 ndim\n";

  t_file.close();
}

void State::randomize()
{
  randomize_vector(q);
  randomize_vector(aux);
}

Solver::Solver()
{
  // Default constructor;
}

Solver::Solver(Solution& solution, int num_ghost, int num_wave):
  num_ghost(num_ghost), num_wave(num_wave)
{
  const int nx = solution.grid.num_cells[0];
  const int ny = solution.grid.num_cells[1];
  const int num_eqn = solution.state.num_eqn;

  //const int dim = solution.grid.dim;
  // Outputs on interfaces
  EdgeFieldIndexer efi(nx, ny, num_ghost, num_eqn);
  EdgeFieldIndexer wave_efi(nx, ny, num_ghost, num_eqn, num_wave);
  EdgeFieldIndexer s_efi(nx, ny, num_ghost, num_wave);
  amdq.resize(efi.size());
  apdq.resize(efi.size());
  wave.resize(wave_efi.size());
  wave_speed.resize(s_efi.size());

  // Default values of member data
  cfl_desired = 0.4;
  cfl_max = 0.5;
  dt_variable = true;
  dt_max = 1e99;
  max_steps = 10000;
  dt = 0.1;
}

Solver::Solver(int* num_cells, int num_eqn, int num_ghost, int num_wave):
  num_ghost(num_ghost), num_wave(num_wave)
{
  const int nx = num_cells[0];
  const int ny = num_cells[1];
  //const int dim = solution.grid.dim;
  // Outputs on interfaces
  EdgeFieldIndexer efi(nx, ny, num_ghost, num_eqn);
  EdgeFieldIndexer wave_efi(nx, ny, num_ghost, num_eqn, num_wave);
  EdgeFieldIndexer s_efi(nx, ny, num_ghost, num_wave);
  amdq.resize(efi.size());
  apdq.resize(efi.size());
  wave.resize(wave_efi.size());
  wave_speed.resize(s_efi.size());

  // Default values of member data
  cfl_desired = 0.4;
  cfl_max = 0.5;
  dt_variable = true;
  dt_max = 1e99;
  max_steps = 10000;
  dt = 0.1;
}

void Solver::define(int* num_cells, int num_eqn, int num_ghost, int num_wave)
{
  this->num_wave = num_wave;
  this->num_ghost = num_ghost;
  const int nx = num_cells[0];
  const int ny = num_cells[1];

  // Outputs on interfaces
  EdgeFieldIndexer efi(nx, ny, num_ghost, num_eqn);
  EdgeFieldIndexer wave_efi(nx, ny, num_ghost, num_eqn, num_wave);
  EdgeFieldIndexer s_efi(nx, ny, num_ghost, num_wave);
  amdq.resize(efi.size());
  apdq.resize(efi.size());
  wave.resize(wave_efi.size());
  wave_speed.resize(s_efi.size());

  if (this->dt_variable)
  {
    FieldIndexer fi(nx, ny, num_ghost, num_eqn);
    q_old.resize(fi.size());
  }
}

int Solver::evolve_to_time(Solution &solution, real t_end)
{
  real t_start = solution.t;
  real cfl;
  FieldIndexer fi(solution.grid.num_cells[0], solution.grid.num_cells[1],
                  num_ghost, solution.state.num_eqn);

  // Check to see if we need to size q_old storage
  if (dt_variable && q_old.size() == 0)
    q_old.resize(fi.size());

  // ========================
  //  Main time stepping loop
  for (int n=0; n < max_steps; ++n)
  {
    // Determine time step if bounded by t_end
    if (solution.t + dt > t_end && t_start < t_end)
      dt = t_end - solution.t;
    if (t_end - solution.t - dt < 1e-14 * solution.t)
      dt = t_end - solution.t;

    // Keep backup of q in case we need to retake a step
    // TODO: Probably should save aux here as well
    if (dt_variable)
    {
      for (unsigned index = 0; index < fi.size(); ++index)
        q_old[index] = solution.state.q[index];
    }
      
    // Take one time step using the solver's current dt
    cfl = step(&solution);

    // Check CFL
    if (cfl <= cfl_max)
    {
      // Accept this time step
      if (dt_variable)
        solution.t += dt;
      else
        // Avoid roundoff error if dt_variable == False
        solution.t = t_start + (n + 1) * dt;
    }
    else
    {
      // Reject this step
      if (dt_variable)
      {
        // Go back to saved q array
        for (unsigned index = 0; index < fi.size(); ++index)
          solution.state.q[index] = q_old[index];
      }
      else
      {
        // Give up, we cannot adapt
        return 1;
      }
    }

    // Choose new time step
    if (dt_variable)
    {
      if (cfl > 0.0)
        dt = fmin(dt_max, dt * cfl_desired / cfl);
      else
        dt = dt_max;
    }
    //  End main time stepping loop
    // =============================

    // Check to see if we are done
    if (solution.t >= t_end)
      break;
  }

  return 0;
}

real Solver::step(Solution *solution)
{
  // Note that this all will break if the grid is not uniform!
  real dtdx = dt / solution->grid.dx[0];

  set_bc(&solution->state.q[0], &solution->state.aux[0],
        solution->grid.num_cells[0], solution->grid.num_cells[1],
        num_ghost, solution->state.num_eqn);

  rp_grid_eval(&solution->state.q[0], &solution->state.aux[0], 
                solution->state.aux_global, solution->grid.num_cells[0], 
                solution->grid.num_cells[1],
               &amdq[0], &apdq[0], &wave[0], &wave_speed[0]);

  update(&solution->state.q[0], &solution->state.aux[0], 
          solution->grid.num_cells[0], solution->grid.num_cells[1],
          &amdq[0], &apdq[0], &wave[0], &wave_speed[0],
          num_ghost, solution->state.num_eqn, dtdx);

  // Compute CFL
  return calculate_cfl(solution->grid.num_cells[0], solution->grid.num_cells[1],
                        num_ghost, solution->state.num_eqn, num_wave,
                        dtdx);
}
