#include <iostream>
#include <sstream>
#include <string>
#include <iomanip>
#include <cmath>
#include <vector>
#include <fstream>
#include <random>
#include "verlet.hpp"
#include "observables.hpp"

int main() {
    constexpr double PI = 3.14159265358979323846;
    int N = 15;                     // número de partículas
    double v_initial = 2.0;         // velocidad inicial
    double radius = 7.0;            // radio esfera 3D

    // Switches para elegir dimensionalidad 

    bool use_1D = false;             // true = 1D
    bool use_2D = true;              // true = 2D (z=0)
    bool use_rotation = false;       // rotación 2D o 3D
    bool perturbation = true;       // perturbación en posiciones
    bool periodicB =true;
    bool enable_walls= false;

    // Condiciones de frontera
    double x_min=-10.0, x_max=10.0;
    double y_min=-10.0, y_max=10.0;
    double z_min=-10.0, z_max=10.0;

    //Condiciones frontera periódica
    double Lx = 10;
    double Ly = 10;
    double Lz = 10;

    // Paso temporal
    double dt = 0.001;
    int steps = 30000;

    // Acoplamientos
    std::vector<double> k_harmonic = {0.0, 0.0, 0.0};  // z=0

    // Vector de partículas
    std::vector<Particle3D> particles(N);

    // Perturbación aleatoria
    double noise_amp = 0.2;
    std::mt19937 gen(42);
    std::uniform_real_distribution<double> dist(-noise_amp, noise_amp);

    // Variables auxiliares 1D
    std::vector<double> pos_init(N);
    double start = -3.5;
    double step  = 1.0;

    for(int i=0; i<N; i++){
        if(use_1D){
            // === Inicialización 1D ===
            if(i==0) pos_init[i] = -8.0;
            else pos_init[i] = start + (i-1)*step;

            particles[i].x = pos_init[i];
            particles[i].y = 0.0;
            particles[i].z = 0.0;

            particles[i].vx = 0.0;
            particles[i].vy = 0.0;
            particles[i].vz = 0.0;
        }
        else if(use_2D){
            // === Inicialización 2D ===
            double phi = 2.0 * PI * i / N; // ángulo para distribuir en círculo en xy
            double base_x = radius * std::cos(phi);
            double base_y = radius * std::sin(phi);
            double base_z = 0.0;           // plano z=0

            if(perturbation){
                particles[i].x = base_x + dist(gen);
                particles[i].y = base_y + dist(gen);
            } else {
                particles[i].x = base_x;
                particles[i].y = base_y;
            }
            particles[i].z = 0.0;           // siempre plano z=0

            if(use_rotation){
                double r_planar = std::sqrt(base_x*base_x + base_y*base_y);
                if(r_planar > 1e-12){
                    particles[i].vx = -v_initial * base_y / r_planar;
                    particles[i].vy =  v_initial * base_x / r_planar;
                 } else {
                    particles[i].vx = particles[i].vy = 0.0;
                }
                particles[i].vz = 0.0;       // velocidad z=0
            } else {
                double mag = std::sqrt(base_x*base_x + base_y*base_y);
                if(mag > 1e-12){
                    particles[i].vx = v_initial * base_x / mag;
                    particles[i].vy = v_initial * base_y / mag;
                } else {
                    particles[i].vx = particles[i].vy = 0.0;
                }
                particles[i].vz = 0.0;       // velocidad z=0
            }
        }
        else {
            // === Inicialización 3D ===
            double phi = std::acos(1.0-2.0 * (i + 0.5) / (double)N);
            double theta = std::sqrt(N * PI) * phi;

            double base_x = radius * std::sin(phi) * std::cos(theta);
            double base_y = radius * std::sin(phi) * std::sin(theta);
            double base_z = radius * std::cos(phi);

            if(perturbation){
                particles[i].x = base_x + dist(gen);
                particles[i].y = base_y + dist(gen);
                particles[i].z = base_z + dist(gen);
            } else {
                particles[i].x = base_x;
                particles[i].y = base_y;
                particles[i].z = base_z;
            }

            if(use_rotation){
                double r_planar = std::sqrt(base_x*base_x + base_y*base_y);
                if(r_planar > 1e-12){
                    particles[i].vx = -v_initial * base_y / r_planar;
                    particles[i].vy =  v_initial * base_x / r_planar;
                    particles[i].vz = 0.0;
                } else {
                    particles[i].vx = particles[i].vy = particles[i].vz = 0.0;
                }
            } else {
                double mag = std::sqrt(base_x*base_x + base_y*base_y + base_z*base_z);
                if(mag > 1e-12){
                    particles[i].vx = v_initial * base_x / mag;
                    particles[i].vy = v_initial * base_y / mag;
                    particles[i].vz = v_initial * base_z / mag;
                }else {
                    particles[i].vx = particles[i].vy = particles[i].vz = 0.0;
                }
            }
        }
    }


    // ------ generar nombres automáticamente para los archivos de datos ---------
    std::stringstream ss;

    // Dimensionalidad
    if(use_1D) {
        ss << "1D_";
    } else if(use_2D) {
        ss << "2D_";
    } else {
        ss << "3D_";
    }

    // Tipo de movimiento
    ss << (use_rotation ? "ROT" : "RAD");

    // Velocidad inicial con 1 decimal
    ss << "_v" << std::fixed << std::setprecision(1) << v_initial;

    // Perturbación
    ss << (perturbation ? "_pert" : "_clean");

    // Fronteras
    ss << (periodicB ? "_periodic" : "_box");

    // Construir nombres de archivo finales
    std::string suffix = ss.str();
    std::string traj_filename = "trayectoria_" + suffix + ".xyz";
    std::string obs_filename  = "observables_" + suffix + ".dat";

    // Ahora traj_filename y obs_filename reflejan correctamente 1D, 2D o 3D

    //guardando los observables
    std::ofstream traj(traj_filename); // Lo mismo que ya teniamos
    std::ofstream obs(obs_filename); // Cambio para obtener y pintar los observables

    // encabezado para el archivo de observables
    obs << "# t K U E\n";

    for(int i = 0; i < steps; i++)
    {
        double t = i * dt;

        velocityVerlet3D(particles, dt, k_harmonic, x_min, x_max, y_min, y_max, z_min, z_max, enable_walls, periodicB, Lx,Ly,Lz);

        // Trayectoria partículas
        traj << particles.size() <<"\n";
        traj << "#t = " << t << "\n";

        for(size_t j = 0; j < particles.size(); j++){
             traj << j << " "
                  << particles[j].x << " "
                  << particles[j].y << " "
                  << particles[j].z << "\n";
        }

        // archivo observables
        double K = kineticEnergy3D(particles);
        double U = potentialEnergy3D(particles, k_harmonic);
        double E = K + U;

        obs << t << " " << K << " " << U << " " << E << "\n";
    }

    traj.close();
    obs.close();
    return 0;
}