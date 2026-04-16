#include "verlet.hpp"
#include <cmath>
#include <random>
#include <omp.h> // ¡Librería mágica para paralelizar!

constexpr double RCUT2 = 16.0; // radio de corte definido arriba en el codigo

void computeAccelerations3D(const std::vector<Particle3D>& particles,
                            std::vector<double>& acc_x,
                            std::vector<double>& acc_y,
                            std::vector<double>& acc_z,
                            const std::vector<double>& k,
                            bool usePBounds,
                            double Lx, double Ly, double Lz)
{   
    int N = particles.size();
    double rcut2 = RCUT2;

    double kx = k[0];
    double ky = k[1];
    double kz = k[2];

    // Trampa armónica (Sobrescribe las aceleraciones, eliminando la necesidad de std::fill)
    #pragma omp parallel for
    for(int i = 0; i < N; i++){
        acc_x[i] = -kx * particles[i].x;
        acc_y[i] = -ky * particles[i].y;
        acc_z[i] = -kz * particles[i].z;
    }

    // Interacción por pares (Cuello de botella O(N^2))
    #pragma omp parallel for schedule(dynamic)
    for(int i = 0; i < N; i++){
        for(int j = i + 1; j < N; j++){

            double dx = particles[i].x - particles[j].x;
            double dy = particles[i].y - particles[j].y;
            double dz = particles[i].z - particles[j].z;

            // minimum image convention
            if(usePBounds){
                dx -= Lx * std::round(dx / Lx);
                dy -= Ly * std::round(dy / Ly);
                dz -= Lz * std::round(dz / Lz);
            }

            double r2 = dx*dx + dy*dy + dz*dz;

            if (r2 < rcut2){

                double r4 = r2*r2;
                double r6 = r4*r2;
                double r14 = r6*r6*r2;

                double f_scalar = 12.0 / r14;

                // Sección atómica para evitar "Race Conditions" (colisiones de hilos)
                #pragma omp atomic
                acc_x[i] += f_scalar*dx;
                #pragma omp atomic
                acc_y[i] += f_scalar*dy;
                #pragma omp atomic
                acc_z[i] += f_scalar*dz;

                #pragma omp atomic
                acc_x[j] -= f_scalar*dx;
                #pragma omp atomic
                acc_y[j] -= f_scalar*dy;
                #pragma omp atomic
                acc_z[j] -= f_scalar*dz;
            }
        }
    }
}

// Implementación velocity verlet
void velocityVerlet3D(std::vector<Particle3D>& particles, double dt,
                      const std::vector<double>& k, double xmin, double xmax,
                      double ymin, double ymax, double zmin, double zmax,
                      bool useBoundaries, bool usePBounds, double Lx, double Ly, double Lz)
{
    int N = particles.size();
    std::vector<double> acc_x(N, 0.0);
    std::vector<double> acc_y(N, 0.0);
    std::vector<double> acc_z(N, 0.0);

    computeAccelerations3D(particles, acc_x, acc_y, acc_z, k, usePBounds, Lx, Ly, Lz);

    #pragma omp parallel for
    for(int i = 0; i < N; i++){
        // velocidades medio paso
        particles[i].vx += 0.5 * acc_x[i] * dt;                 
        particles[i].vy += 0.5 * acc_y[i] * dt;
        particles[i].vz += 0.5 * acc_z[i] * dt;
        
        // posiciones paso entero
        particles[i].x += particles[i].vx * dt;                 
        particles[i].y += particles[i].vy * dt;
        particles[i].z += particles[i].vz * dt;
    }

    if (useBoundaries){
        applyReflectiveBC3D(particles, xmin, xmax, ymin, ymax, zmin, zmax);     
    }
    if (usePBounds){
        applyPeriodicBoundary(particles, Lx, Ly, Lz);  
    }
    
    computeAccelerations3D(particles, acc_x, acc_y, acc_z, k, usePBounds, Lx, Ly, Lz);         

    #pragma omp parallel for
    for(int i = 0; i < N; i++){
        particles[i].vx += 0.5 * acc_x[i] * dt;                 
        particles[i].vy += 0.5 * acc_y[i] * dt;
        particles[i].vz += 0.5 * acc_z[i] * dt;
    }
}


void applyReflectiveBC3D(std::vector<Particle3D>& particles,
                         double xmin, double xmax,
                         double ymin, double ymax,
                         double zmin, double zmax)
{
    int N = particles.size();
    #pragma omp parallel for
    for (int i = 0; i < N; i++) {
        auto& p = particles[i]; // Referencia para mantener el código limpio

        // eje x
        if (p.x < xmin) {
            p.x = 2.0 * xmin - p.x;
            p.vx = -p.vx;
        }
        if (p.x > xmax) {
            p.x = 2.0 * xmax - p.x;
            p.vx = -p.vx;
        }

        // eje y
        if (p.y < ymin) {
            p.y = 2.0 * ymin - p.y;
            p.vy = -p.vy;
        }
        if (p.y > ymax) {
            p.y = 2.0 * ymax - p.y;
            p.vy = -p.vy;
        }

        // eje z
        if (p.z < zmin) {
            p.z = 2.0 * zmin - p.z;
            p.vz = -p.vz;
        }
        if (p.z > zmax) {
            p.z = 2.0 * zmax - p.z;
            p.vz = -p.vz;
        }

    }
}




void applyPeriodicBoundary(std::vector<Particle3D>& particles,
                           double Lx, double Ly, double Lz)
{
    int N = particles.size();
    #pragma omp parallel for
    for (int i = 0; i < N; i++) {
        auto& p = particles[i];

        if (p.x >= Lx) p.x -= Lx;
        else if (p.x < 0.0) p.x += Lx;

        if (p.y >= Ly) p.y -= Ly;
        else if (p.y < 0.0) p.y += Ly;

        if (p.z >= Lz) p.z -= Lz;
        else if (p.z < 0.0) p.z += Lz;
    }
}

// Distribución condiciones iniciales
bool tooClose(const std::vector<Particle3D>& particles,
              double x, double y, double z,
              int current,
              double minDist,
              bool usePBounds,
              double Lx, double Ly, double Lz)
{
    double minDist2 = minDist * minDist;

    // Este ciclo es secuencial porque tiene un "return" temprano
    for(int j = 0; j < current; j++){

        double dx = x - particles[j].x;
        double dy = y - particles[j].y;
        double dz = z - particles[j].z;

        if(usePBounds){
            dx -= Lx * std::round(dx / Lx);
            dy -= Ly * std::round(dy / Ly);
            dz -= Lz * std::round(dz / Lz);
        }

        double r2 = dx*dx + dy*dy + dz*dz;

        if(r2 < minDist2)
            return true;
    }

    return false;
}

// Termostato de Andersen
void applyAndersenThermostat(std::vector<Particle3D>& particles,
                             double T_target,
                             double nu,
                             double dt,
                             int dim,
                             std::mt19937& gen)
{
    std::uniform_real_distribution<double> dist_prob(0.0, 1.0);
    std::normal_distribution<double> dist_vel(0.0, std::sqrt(T_target));
    
    double collision_prob = nu * dt;
    
    // Secuencial: El generador aleatorio mt19937 no soporta multihilo (no es thread-safe)
    for(auto& p : particles){
        if(dist_prob(gen) < collision_prob){
            if(dim >= 1) p.vx = dist_vel(gen);
            if(dim >= 2) p.vy = dist_vel(gen);
            if(dim == 3) p.vz = dist_vel(gen);
        }
    }
}