//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) Mara Strungaru 2022. All rights reserved.
//
//   Email: mara.strungaru@york.ac.uk
//
//   implementation based on the paper Phys. Rev. B 103, 024429, (2021) M.Strungaru, M.O.A. Ellis et al
//------------------------------------------------------------------------------
//

// C++ standard library headers

// Vampire headers
#include "sld.hpp"
#include <iostream>
#include <vector>
#include <math.h>
#include "constants.hpp"
#include "material.hpp"

// sld module headers
#include "internal.hpp"


namespace sld{


   double compute_spin_temperature(const int start_index, // first atom for exchange interactions to be calculated
               const int end_index,
               const std::vector<int>& type_array, // type for atom
               std::vector<double>& x_spin_array, // coord vectors for atoms
               std::vector<double>& y_spin_array,
               std::vector<double>& z_spin_array,
               std::vector<double>& fields_array_x, //  vectors for fields
               std::vector<double>& fields_array_y,
               std::vector<double>& fields_array_z,
               std::vector <double>& mu_s_array){

                double SxH2=0.0;
                double SH=0.0;
                for (int at=start_index;at<end_index;at++){
                    double Sx = x_spin_array[at];
                    double Sy = y_spin_array[at];
                    double Sz = z_spin_array[at];
                    double Hx = fields_array_x[at];
                    double Hy = fields_array_y[at];
                    double Hz = fields_array_z[at];

                     double SxHx = Sy * Hz - Sz * Hy;
                     double SxHy = Sz * Hx - Sx * Hz;
                     double SxHz = Sx * Hy - Sy * Hx;
                     SxH2  = SxH2+ SxHx*SxHx + SxHy*SxHy + SxHz*SxHz;
                     SH  = SH +  Sx * Hx + Sy * Hy + Sz*Hz;

                }
               double T_spin=0.5* mu_s_array[0]*constants::muB /constants::kB * SxH2 / SH;

      return T_spin;

      }//end of spin_temperature


double compute_lattice_temperature(const int start_index, // first atom for exchange interactions to be calculated
            const int end_index,
            const std::vector<int>& type_array, // type for atom
            std::vector<double>& velo_array_x, // coord vectors for atoms
            std::vector<double>& velo_array_y,
            std::vector<double>& velo_array_z){

            double kinetic=0;
             for (int at=start_index;at<end_index;at++){
                 double vx = velo_array_x[at];
                 double vy = velo_array_y[at];
                 double vz = velo_array_z[at];
                 kinetic+=vx*vx+vy*vy+vz*vz;
                 
                 //if (at==2) std::cout<<"SLD at 2 "<<kinetic*sld::internal::mp[0].mass.get()*0.5*( 2.0 /(3.0*constants::kB_eV))<<std::endl;
             }

            kinetic*=sld::internal::mp[0].mass.get()*0.5/(end_index-start_index);
            double T_lat=( 2.0 /(3.0*constants::kB_eV))*kinetic;

   return T_lat;

}//end of lattice_temperature

// angular momentum - Nargiz Ahmadzada
void compute_angular_momentum(
    const int start_index, // first atom for exchange interactions to be calculated
    const int end_index,
    const std::vector<int>& type_array, // type for atom,
    std::vector<double>& x_coord_array, // coord vectors for atoms (r)
    std::vector<double>& y_coord_array,
    std::vector<double>& z_coord_array,
    std::vector<double>& velo_array_x, // velocity vectors for atoms
    std::vector<double>& velo_array_y,
    std::vector<double>& velo_array_z,
    std::array<double, 4>& L){
    
    double Rx, Ry, Rz; // r vector components
    double Vx, Vy, Vz; // velocity vector components
    double Lx = 0.0;
    double Ly = 0.0;
    double Lz = 0.0; // angular momentum vector components
                    
    for (int at=start_index;at<end_index;at++){
      const unsigned int imat = type_array[at];
      double mass = sld::internal::mp[imat].mass.get(); // mass
      Rx = x_coord_array[at];
      Ry = y_coord_array[at];
      Rz = z_coord_array[at];
      Vx = velo_array_x[at];
      Vy = velo_array_y[at];
      Vz = velo_array_z[at];

      Lx += mass * (Ry * Vz - Rz * Vy); // m * R x V (cross product)
      Ly += mass * (Rz * Vx - Rx * Vz);
      Lz += mass * (Rx * Vy - Ry * Vx);
                        
    }
    L[0] = Lx;
    L[1]= Ly;
    L[2]= Lz;
    L[3] = sqrt(Lx*Lx + Ly*Ly + Lz*Lz); // update to module
	return;
  } 
// end of angular momentum

   } // end of sld namespace
     //
