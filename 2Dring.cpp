#include <iostream>
#include <cmath>
#include <vector>
#include <fstream>
#include <random>
#include "particle.hpp"
#include "verlet.hpp"
#include "observables.hpp"

int main() {
    constexpr double PI = 3.14159265358979323846;
    int N = 10;                     // Número de partículas deben ser 10
    
    //condiciones experimento
    double k_harmonic = 0.0;       //poner en 0.0 para remover trampa armonica (mesa de pool)
    double radius = 4.0;           //radio al que se ajustaran las particulas 
    double v_initial = 2.0;        //velocidad inicial de las particulas

    //switches
    bool enable_walls = true;      //activar/desactivar condiciones de frontera reflectivas
    bool use_rotation = true;      //true: particulas rotan al rededor del origen/false: particulas son disparadas hacia afuera 
    bool perturbation = true;      //activar/desactivar el cambio en condiciones iniciales


    //condiciones de frontera
    double x_min = -7.0, x_max = 7.0;
    double y_min = -7.0, y_max = 7.0;

    //intervalo de integracion y numero de pasos
    double dt = 0.001;              //paso de tiempo
    int steps = 30000;              //deben ser 30 mil
    std::vector<Particle2D> particles(N);

    // condiciones iniciales, probablemente no lo mas optimizado posiblemente, pero ahí camella 
    std::vector<double> pos_init_x(N);
    std::vector<double> pos_init_y(N);
    std::vector<double> theta(N);

    //generador pequeña perturbacion de la posicion 
    double noise_amp = 0.15;
    std::mt19937 gen(42);
    std::uniform_real_distribution<double> dist(-noise_amp,noise_amp);

    //partir el circulo en N angulos iguales
    for(int i = 0; i < N; i++){
        //calcular posicion angular
        double theta_i = 2.0 * PI * i / (double)N;

        //R = 4, mapeamos coordendas polares a cartesianas
        double base_x = radius * std::cos(theta_i);
        double base_y = radius * std::sin(theta_i);
    
        //aplica perturbacion a la posicion si perturbation = true
        if (perturbation){
            particles[i].x = base_x + dist(gen);
            particles[i].y = base_y + dist(gen);
        } else {
            particles[i].x = base_x;
            particles[i].y = base_y;
        }
        
        //definicion de velocidades
        if (use_rotation){
            particles[i].vx = -v_initial * std::sin(theta_i); 
            particles[i].vy = v_initial * std::cos(theta_i);
        } else {
            particles[i].vx = v_initial * std::cos(theta_i); 
            particles[i].vy = v_initial * std::sin(theta_i);
        }
        
        std::cout<<particles[i].x<<std::endl;
        
    }

    

    //guardando los observables
    std::ofstream traj("trayectoria2DPTchaos.dat"); // Lo mismo que ya teniamos
    std::ofstream obs("observables2DPTchaos.dat"); // Cambio para obtener y pintar los observables

    // encabezado para el archivo de observables
    obs << "# t K U E\n";

    for(int i = 0; i < steps; i++)
    {
        double t = i * dt;

        velocityVerlet2D(particles, dt, k_harmonic, x_min, x_max, y_min, y_max, enable_walls);

        // Trayectoria partículas
        traj << particles.size() <<"\n";
        traj << "#t = " << t << "\n";

        for(size_t j = 0; j < particles.size(); j++){
             traj << j << " "
                  << particles[j].x << " "
                  << particles[j].y<< " "
                  << 0.0 << "\n";
        }

        // archivo observables
        double K = kineticEnergy2D(particles);
        double U = potentialEnergy2D(particles, k_harmonic);
        double E = K + U;

        obs << t << " " << K << " " << U << " " << E << "\n";
    }

    traj.close();
    obs.close();
    return 0;
}