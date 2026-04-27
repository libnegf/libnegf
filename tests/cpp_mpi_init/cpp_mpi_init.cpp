/*!!--------------------------------------------------------------------------!
 *!! libNEGF: a general library for Non-Equilibrium Greens functions.         !
 *!! Copyright (C) 2012 - 2026                                                !
 *!!                                                                          !
 *!! This file is part of libNEGF: a library for                              !
 *!! Non Equilibrium Green's Function calculation                             !
 *!!                                                                          !
 *!! Developers: Alessandro Pecchia, Daniele Soccodato                        !
 *!! Former Contributors: Gabriele Penazzi, Luca Latessa, Aldo Di Carlo       !
 *!!                                                                          !
 *!! libNEGF is free software: you can redistribute it and/or modify          !
 *!! it under the terms of the GNU Lesse General Public License as published  !
 *!! by the Free Software Foundation, either version 3 of the License, or     !
 *!! (at your option) any later version.                                      !
 *!!                                                                          !
 *!!  You should have received a copy of the GNU Lesser General Public        !
 *!!  License along with libNEGF.  If not, see                                !
 *!!  <http://www.gnu.org/licenses/>.                                         !
 *!!--------------------------------------------------------------------------!
 */
#include "libnegf.hpp"

#include <mpi.h>

#include <cstdio>
#include <cstring>


int main() {
    int ierr = MPI_Init(nullptr, nullptr);
    if(ierr != 0) {
        fprintf(stderr, "Error in mpi_init: %d\n", ierr);
        return 1;
    }

    printf("Initializing libNEGF\n");
    int handler[NEGF_HSIZE];
    negf_init_session(handler);
    negf_init(handler);

    MPI_Fint global_comm_f = MPI_Comm_c2f(MPI_COMM_WORLD);
    negf_set_mpi_fcomm(handler, global_comm_f);
    MPI_Fint cart_comm, k_comm, en_comm;
    negf_cartesian_init(handler, global_comm_f, 1, cart_comm, k_comm, en_comm);

    //Release library
    negf_destruct_libnegf(handler);
    negf_destruct_session(handler);
    printf("Done\n");

    MPI_Finalize();
}
