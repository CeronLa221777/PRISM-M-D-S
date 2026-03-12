#include <iostream>
#include <sstream>
#include <string>
#include <iomanip>
#include <cmath>
#include <vector>
#include <fstream>
#include <random>
#include "particle.hpp"
#include "verlet.hpp"
#include "observables.hpp"

int main() {
    constexpr double PI = 3.14159265358979323846;
    int N = 50;                     // mas particulas para comportamiento mas interesante
    
    //condiciones experimento
    double k_harmonic = 0.0;       //poner en 0.0 para remover trampa armonica (mesa de pool)
    double radius = 7.0;           //radio de la esfera a la que se ajustaran las particulas 
    double v_initial = 2.0;        //velocidad inicial de las particulas

    //switches
    bool enable_walls = true;      //activar/desactivar condiciones de frontera reflectivas
    bool use_rotation = true;      //true: particulas rotan al rededor del origen/false: particulas son disparadas hacia afuera 
    bool perturbation = false;      //activar/desactivar el cambio en condiciones iniciales


    //condiciones de frontera
    double x_min = -10.0, x_max = 10.0;
    double y_min = -10.0, y_max = 10.0;
    double z_min = -10.0, z_max = 10.0;


    //intervalo de integracion y numero de pasos
    double dt = 0.001;              //paso de tiempo
    int steps = 30000;              //deben ser 30 mil
    std::vector<Particle3D> particles(N);

    //generador pequeña perturbacion de la posicion 
    double noise_amp = 0.2;
    std::mt19937 gen(42);
    std::uniform_real_distribution<double> dist(-noise_amp,noise_amp);

    //poner las particulas en puntos al rededor de una esfera
    for(int i = 0; i < N; i++){
        //calcular posicion angular
        double phi = std::acos(1.0-2.0 * (i + 0.5) / (double)N);
        double theta = std::sqrt(N * PI) * phi;

        double base_x = radius * std::sin(phi) * std::cos(theta);
        double base_y = radius * std::sin(phi) * std::sin(theta);
        double base_z = radius * std::cos(phi);
    
        //aplica perturbacion a la posicion si perturbation = true
        if (perturbation){
            particles[i].x = base_x + dist(gen);
            particles[i].y = base_y + dist(gen);
            particles[i].z = base_z + dist(gen);
        } else {
            particles[i].x = base_x;
            particles[i].y = base_y;
            particles[i].z = base_z;
        }
        
        //definicion de velocidades
        if (use_rotation){
            //rotacion pura alrededor del eje z
            //manteniendo la magnitud v_inicial por medio de normalizacion en el plano xy
            double r_planar = std::sqrt(base_x*base_x + base_y*base_y);
            if (r_planar > 0){
                particles[i].vx = -v_initial * (base_y / r_planar);
                particles[i].vy =  v_initial * (base_x / r_planar);
                particles[i].vz = 0.0; //se mantiene la altura z
            } else {
                particles[i].vx = particles[i].vy = particles[i].vz = 0.0;
            }
        } else {
            //disparo radial (desde el centro hacia afuera)
            double mag = std::sqrt(base_x*base_x + base_y*base_y + base_z*base_z);
            particles[i].vx = v_initial * (base_x) / mag; 
            particles[i].vy = v_initial * (base_y) / mag;
            particles[i].vz = v_initial * (base_z) / mag;
        }
        
        //std::cout<<particles[i].x<<std::endl;
    }

    //------generar nombres automaticamente para los archivos de datos ---------
    std::stringstream ss;
    ss << "3D_" << (use_rotation ? "ROT" : "RAD")
       << "_v" << std::fixed << std::setprecision(1) << v_initial
       << (perturbation ? "_pert" : "_clean");

    std::string suffix = ss.str();
    std::string traj_filename = "trayectoria_" + suffix + ".dat";
    std::string obs_filename = "observables_" + suffix + ".dat";
    //-------------------------------------------------------------------------

    //guardando los observables
    std::ofstream traj(traj_filename); // Lo mismo que ya teniamos
    std::ofstream obs(obs_filename); // Cambio para obtener y pintar los observables

    // encabezado para el archivo de observables
    obs << "# t K U E\n";

    for(int i = 0; i < steps; i++)
    {
        double t = i * dt;

        velocityVerlet3D(particles, dt, k_harmonic, x_min, x_max, y_min, y_max, z_min, z_max, enable_walls);

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