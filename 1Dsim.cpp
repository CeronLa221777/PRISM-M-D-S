#include <iostream>
#include <vector>
#include <fstream>
#include "particle.hpp"
#include "verlet.hpp"
#include "observables.hpp"

int main() {
    int N = 10;                                 //Número de partículas: deben ser 10
    double dt = 0.001;                          //intervalo de tiempo: debería ser 0.001
    int steps = 30000;                          //numero de pasos de integracion: deben ser 30 mil
    std::vector<Particle1D> particles(N);       //vector que almacena las posiciones y velocidades de cada particula para N particulas 
    std::vector<double> pos_init(N);            //posiciones iniciales de las particulas

    // condiciones iniciales---- Corrección, condiciones iniciales sugeridas por el profesor 
    double start = -3.5;
    double step  = 1.0;


    pos_init[0] = -8.0;  // la partícula aislada

    for(int i = 1; i < N; i++){
        pos_init[i] = start + (i-1)*step;
    }
    for(int i = 0; i < N; i++){
        particles[i].x =pos_init[i];
        particles[i].v = 0.0; 
    }

    std::ofstream traj("trayectoria1Dregen.dat");    //archivo de datos que guarda informacion de  
    std::ofstream obs("observables1Dregen.dat");     //Cambio para obtener y pintar los observables

    // encabezado para el archivo de observables
    obs << "# t K U E\n";

    for(int i = 0; i < steps; i++)
    {
        double t = i * dt;

        velocityVerlet1D(particles, dt);

        // Trayectoria partículas
        traj << "\n";
        traj << "t = " << t << "\n";

        for(size_t j = 0; j < particles.size(); j++){
             traj << j << " "
                  << particles[j].x << " "
                  << 0.0 << " "
                  << 0.0 << "\n";
        }

        // archivo observables
        double K = kineticEnergy1D(particles);
        double U = potentialEnergy1D(particles);
        double E = K + U;

        obs << t << " "
            << K << " "
            << U << " "
            << E << "\n";
    }

    traj.close();
    obs.close();
    return 0;
}