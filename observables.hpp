#ifndef OBSERVABLES_HPP
#define OBSERVABLES_HPP

#include <vector>
#include "particle.hpp"
//aquí declararemos las variables de los observables a calcular

// Energía cinética
double kineticEnergy1D(const std::vector<Particle1D>& particles); // const hace que calcule el observable del sistema pero que no lo cambie, el & es para que no almacene el vector, sino que solo se usen las componentes 
//std:: es el prefijo del namespace estándar de C++.
//Se usa para decirle al compilador de dónde viene un símbolo (función, clase, objeto)
// Energía potencial (forma funcional editable)
double potentialEnergy1D(const std::vector<Particle1D>& particles, double k);

double kineticEnergy2D(const std::vector<Particle2D>& particles);
double potentialEnergy2D(const std::vector<Particle2D>& particles, double k);

double kineticEnergy3D(const std::vector<Particle3D>& particles);
double potentialEnergy3D(const std::vector<Particle3D>& particles, double k);


#endif