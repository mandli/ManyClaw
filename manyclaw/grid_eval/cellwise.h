#ifndef CELLWISE_H
#define CELLWISE_H

struct RiemannCellwiseStepper {
  inline int operator()(const rp_t rp, const rp_grid_params grid_params,
                        const real* q,  const real* aux,
                        const void* aux_global,
                        const int nx, const  int ny,
                        real* amdq, real* apdq, real* wave,
                        real* wave_speeds){
    int col, row, idx_left, idx_center, idx_up, idx_out_x, idx_out_y;
    const int num_ghost = grid_params.num_ghost;
    const int num_eqn = grid_params.num_eqn;
    const int num_waves = grid_params.num_waves;

    for(row = num_ghost; row <= ny + num_ghost; ++row) {
      for(col = num_ghost; col <= nx + num_ghost; ++col) {
        idx_left = col + row*(nx + 2*num_ghost) - 1;
        idx_up = col + (row - 1)*(nx + 2*num_ghost);
        idx_center = idx_left + 1;
        idx_out_x = (col - num_ghost) + (row - num_ghost) * (nx + 1);
        idx_out_y = idx_out_x + ((nx + 1)*(ny + 1));
        
        rp(q + idx_left*num_eqn, q + idx_center*num_eqn,
           aux, aux,  aux_global,
           amdq + idx_out_x*num_eqn, apdq + idx_out_x*num_eqn,
           wave + num_waves*num_eqn*idx_out_x, wave_speeds + num_waves*idx_out_x);
        rp(q + idx_up*num_eqn, q + idx_center*num_eqn,
           aux, aux,  aux_global,
           amdq + idx_out_y*num_eqn, apdq + idx_out_y*num_eqn,
           wave + num_waves*num_eqn*idx_out_y, wave_speeds + num_waves*idx_out_y);
      }
    }
    return 0;
  }
};

#endif // CELLWISE_H
