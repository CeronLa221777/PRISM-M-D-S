#include "observables.hpp"
#include <cmath>
#include <fstream>
#include <iostream>
#include <algorithm> // para std::min
#include <omp.h>     // Librería OpenMP para paralelización

constexpr double RCUT2 = 16.0;
constexpr double PI = 3.14159265358979323846; 

//---------3D-------------------------------------------
double kineticEnergy3D(const std::vector<Particle3D>& particles)
{
    double K = 0.0;
    int N = particles.size();
    
    // Paralelizamos y sumamos de forma segura en la variable K
    #pragma omp parallel for reduction(+:K)
    for(int i = 0; i < N; i++) {
        K += 0.5 * (particles[i].vx * particles[i].vx + 
                    particles[i].vy * particles[i].vy + 
                    particles[i].vz * particles[i].vz);
    }
    return K;
}

double potentialEnergy3D(const std::vector<Particle3D>& particles,
                         const std::vector<double>& k,
                         bool usePBounds,
                         double Lx, double Ly, double Lz)
{
    double U = 0.0;
    int N = particles.size();

    double kx = k[0];
    double ky = k[1];
    double kz = k[2];

    double rcut2 = RCUT2;

    // Energía de la trampa armónica
    #pragma omp parallel for reduction(+:U)
    for(int i = 0; i < N; i++){
        U += 0.5 * (kx*particles[i].x*particles[i].x + 
                    ky*particles[i].y*particles[i].y + 
                    kz*particles[i].z*particles[i].z);
    }

    // Energía de interacción por pares (Cuello de botella)
    #pragma omp parallel for schedule(dynamic) reduction(+:U)
    for(int i = 0; i < N; i++){
        for(int j = i + 1; j < N; j++){

            double dx = particles[i].x - particles[j].x;
            double dy = particles[i].y - particles[j].y;
            double dz = particles[i].z - particles[j].z;

            // minimum image
            if(usePBounds){
                dx -= Lx * std::round(dx / Lx);
                dy -= Ly * std::round(dy / Ly);
                dz -= Lz * std::round(dz / Lz);
            }

            double r2 = dx*dx + dy*dy + dz*dz;

            if(r2 < rcut2){
                double r6 = r2*r2*r2;
                double r12 = r6*r6;
                U += 1.0 / r12;
            }
        }
    }

    return U;
}

// ------------------------------------------------------------------
// FUNCIÓN 1: Recolector del Histograma g(r)
// ------------------------------------------------------------------
void updateRDF3D(const std::vector<Particle3D>& particles,
                 std::vector<double>& rdf_hist,
                 double dr,
                 bool usePBounds,
                 double Lx, double Ly, double Lz)
{
    int N = particles.size();
    
    // Por convención de imagen mínima, no tiene sentido medir más allá de L/2
    double max_dist = std::min({Lx, Ly, Lz}) / 2.0;

    // Paralelizamos la recolección de distancias
    #pragma omp parallel for schedule(dynamic)
    for (int i = 0; i < N; i++) {
        for (int j = i + 1; j < N; j++) {
            
            double dx = particles[i].x - particles[j].x;
            double dy = particles[i].y - particles[j].y;
            double dz = particles[i].z - particles[j].z;

            // minimum image convention
            if (usePBounds) {
                dx -= Lx * std::round(dx / Lx);
                dy -= Ly * std::round(dy / Ly);
                dz -= Lz * std::round(dz / Lz);
            }

            double r = std::sqrt(dx*dx + dy*dy + dz*dz);

            // Solo agregamos si está dentro de la caja mínima medible
            if (r < max_dist) {
                // Cambiado a size_t para solucionar el warning del compilador
                size_t bin = static_cast<size_t>(r / dr); 
                
                if (bin < rdf_hist.size()) {
                    // atomic previene que dos hilos intenten sumar al mismo bin al mismo tiempo
                    #pragma omp atomic
                    rdf_hist[bin] += 2.0; 
                }
            }
        }
    }
}

// ------------------------------------------------------------------
// FUNCIÓN 2: Normalizador y Escritura de g(r)
// ------------------------------------------------------------------
void normalizeAndSaveRDF3D(const std::vector<double>& rdf_hist,
                           const std::string& filename,
                           int N, double Lx, double Ly, double Lz,
                           int num_snapshots, double dr)
{
    // Esta función es muy rápida y solo se llama una vez al final. 
    // No requiere paralelización.
    
    std::ofstream out(filename);
    if (!out.is_open()) {
        std::cerr << "[ERROR] No se pudo crear el archivo " << filename << std::endl;
        return;
    }

    double V = Lx * Ly * Lz;
    double rho = N / V; // Densidad ideal

    out << "# r g(r)\n";

    for (size_t bin = 0; bin < rdf_hist.size(); bin++) {
        double r_inner = bin * dr;
        double r_outer = (bin + 1) * dr;
        double r_center = r_inner + 0.5 * dr;

        // Volumen del cascarón esférico V = 4/3 * PI * (r_outer^3 - r_inner^3)
        double dV = (4.0 / 3.0) * PI * (r_outer * r_outer * r_outer - r_inner * r_inner * r_inner);
        
        // ¿Cuántas partículas habría en este cascarón si fuera un gas ideal sin interacciones?
        double ideal_particles = rho * dV;

        // Normalizamos el conteo usando los snapshots totales tomados, las N particulas y el gas ideal
        double g_r = 0.0;
        if (ideal_particles > 0.0) {
            g_r = rdf_hist[bin] / (N * num_snapshots * ideal_particles);
        }

        out << r_center << " " << g_r << "\n";
    }
    
    out.close();
}