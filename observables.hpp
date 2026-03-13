#ifndef OBSERVABLES_HPP
#define OBSERVABLES_HPP

#include <vector>
#include "verlet.hpp"
//aquí declararemos las variables de los observables a calcular

double kineticEnergy3D(const std::vector<Particle3D>& particles);
double potentialEnergy3D(const std::vector<Particle3D>& particles, const std::vector<double>& k);


#endif