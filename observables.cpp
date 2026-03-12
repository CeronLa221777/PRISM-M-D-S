#include "observables.hpp"
#include <cmath>

//aquí definimos nuestra librería de funciones para calcular observables
// -------1D------------------
// ---- Energía cinética ----

constexpr double RCUT2 = 16;

double kineticEnergy1D(const std::vector<Particle1D>& particles)
{
    double K = 0.0;
    for(const auto& p : particles)// auto le dice al compilador "infiero el tiop p automaticamente a partir de particles", ahorra escritura y ya :|
        K += 0.5 * p.v * p.v; // masa = 1
    return K;
}

// ---- Energía potencial ----
double potentialEnergy1D(const std::vector<Particle1D>& particles, double k)
{
    double U = 0.0;
    int N = particles.size();

    // Energía de la trampa armónica
    for(const auto& p : particles){
        U += 0.5 * k * p.x * p.x;
    }
    // Energía por interacción de pares (soft-core)
    for(int i = 0; i < N; i++){
        for(int j = i+1; j < N; j++){
            double dx = particles[i].x - particles[j].x;
            double r2 = dx*dx;
            double rcut2 = RCUT2; // mismo corte que en aceleraciones
            if(r2 < rcut2){
                double r6 = r2 * r2 * r2;
                U += 1.0 / (r6 * r6); // energía potencial de pares
            }
        }
    }

    return U;
}



//---------2D--------------------------------

double kineticEnergy2D(const std::vector<Particle2D>& particles)
{
    double K = 0.0;
    for(const auto& p : particles)// auto le dice al compilador "infiero el tiop p automaticamente a partir de particles", ahorra escritura y ya :|
        K += 0.5 *( p.vx * p.vx + p.vy * p.vy); // masa = 1
    return K;
}

double potentialEnergy2D(const std::vector<Particle2D>& particles, double k)
{
    double U = 0.0;
    int N = particles.size();

    // Energía de la trampa armónica
    for(const auto& p : particles){
        U += 0.5 * k * (p.x * p.x + p.y * p.y );
    }
    // Energía por interacción de pares (soft-core)
    for(int i = 0; i < N; i++){
        for(int j = i+1; j < N; j++){
            double dx = particles[i].x - particles[j].x;
            double dy = particles[i].y - particles[j].y;
            double r2 = dx*dx + dy*dy;
            double rcut2 = RCUT2; // mismo corte que en aceleraciones
            if(r2 < rcut2){
                double r6 = r2 * r2 * r2;
                U += 1.0 / (r6 * r6); // energía potencial de pares
            }
        }
    }

    return U;
}

//---------3D-------------------------------------------
double kineticEnergy3D(const std::vector<Particle3D>& particles)
{
    double K = 0.0;
    for(const auto& p : particles)// auto le dice al compilador "infiero el tiop p automaticamente a partir de particles", ahorra escritura y ya :|
        K += 0.5 *( p.vx * p.vx + p.vy * p.vy + p.vz * p.vz); // masa = 1
    return K;
}

double potentialEnergy3D(const std::vector<Particle3D>& particles, double k)
{
    double U = 0.0;
    int N = particles.size();

    // Energía de la trampa armónica
    for(const auto& p : particles){
        U += 0.5 * k * (p.x * p.x + p.y * p.y + p.z * p.z);
    }
    // Energía por interacción de pares (soft-core)
    for(int i = 0; i < N; i++){
        for(int j = i+1; j < N; j++){
            double dx = particles[i].x - particles[j].x;
            double dy = particles[i].y - particles[j].y;
            double dz = particles[i].z - particles[j].z;
            double r2 = dx*dx + dy*dy + dz*dz;
            double rcut2 = RCUT2; // mismo corte que en aceleraciones
            if(r2 < rcut2){
                double r6 = r2 * r2 * r2;
                U += 1.0 / (r6 * r6); // energía potencial de pares
            }
        }
    }

    return U;
}






