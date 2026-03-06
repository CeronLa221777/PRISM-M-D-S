#ifndef VERLET_HPP
#define VERLET_HPP

#include <vector>
#include "particle.hpp"
void computeAccelerations1D(const std::vector<Particle1D>& particles,
                            std::vector<double>& acc);
void computeAccelerations2D(const std::vector<Particle2D>& particles,
                            std::vector<double>& acc_x,
                            std::vector<double>& acc_y,
                            double k);
void velocityVerlet1D(std::vector<Particle1D> &particles, double dt);
void velocityVerlet2D(std::vector<Particle2D> &particles,
                        double dt,
                        double k,
                        double xmin, double xmax,
                        double ymin, double ymax,
                        bool useBoundaries);
void applyReflectiveBC1D(std::vector<Particle1D>& particles,
                         double xmin,
                         double xmax);
void applyReflectiveBC2D(std::vector<Particle2D>& particles,
                         double xmin, double xmax,
                         double ymin, double ymax);                         
//std::vec..bla bla bla es el sistema de N partículas, & se pasa por ref no se copia memoria Estos es la interfaz del integrador
#endif