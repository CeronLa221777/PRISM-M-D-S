#include "verlet.hpp"
#include <cmath>

void computeAccelerations1D(const std::vector<Particle1D>& particles,
                          std::vector<double>& acc) //Nos toca calcular por aparte las aceleraciones para no todo meterlo en el verlet
{
    int N = particles.size();

    // Trampa armónica
    for(int i = 0; i < N; i++)
        acc[i] = -0.3*particles[i].x;

    // Interacción por pares
    for(int i = 0; i < N; i++){
        for(int j = i+1; j < N; j++){

            double dx = particles[i].x - particles[j].x; // Vector relativo
            double r2 = dx*dx; // Norma cuadrado del vector relativo
            double rcut2 = 16;
            if (r2 < rcut2){
                double f = 12.0 * dx / std::pow(r2,7); // spft-core luego de la derivada da elevado a la 14
                acc[i] += f;// fuerza que ejerce j sobre i
                acc[j] -= f;// línea que hace que no tengamos cálculos duplicados, fuerza que recibe j de i por tecera ley de newton
            }
            
        }
    }
}
void computeAccelerations2D(const std::vector<Particle2D>& particles,
                            std::vector<double>& acc_x,
                            std::vector<double>& acc_y,
                            double k) //Nos toca calcular por aparte las aceleraciones para no todo meterlo en el verlet
{
    int N = particles.size();
    double rcut2 = 16;
    // Trampa armónica
    for(int i = 0; i < N; i++){
        acc_x[i] = -k*particles[i].x;
        acc_y[i] = -k*particles[i].y;
    }
    // Interacción por pares (soft-core)
    for(int i = 0; i < N; i++){
        for(int j = i+1; j < N; j++){

            double dx = particles[i].x - particles[j].x; // Vector relativo
            double dy = particles[i].y - particles[j].y; 
            double r2 = dx*dx + dy*dy; // Norma cuadrado del vector relativo
            
            if (r2 < rcut2){
                double f_scalar = 12.0  / std::pow(r2,7); // soft-core luego de la derivada da elevado a la 14
                double fx = f_scalar * dx;
                double fy = f_scalar * dy;
                acc_x[i] += fx;// fuerza que ejerce j sobre i
                acc_y[i] += fy;


                acc_x[j] -= fx;// línea que hace que no tengamos cálculos duplicados, fuerza que recibe j de i por tecera ley de newton
                acc_y[j] -= fy;
            }
            
        }
    }
}

void velocityVerlet1D(std::vector<Particle1D>& particles, double dt)
{
    int N = particles.size();
    std::vector<double> acc(N);

    computeAccelerations1D(particles, acc); // Calculamos aceleración en configuración "inicial"
    for(int i = 0; i < N; i++){
        particles[i].v += 0.5 * acc[i] * dt;// afectamos la velocidad de la partícula
        particles[i].x += particles[i].v * dt;// corregimos la posición 
    }

    computeAccelerations1D(particles, acc); // Calculamos la aceleración con la nueva configuración de posiciones

    for(int i = 0; i < N; i++){
        particles[i].v += 0.5 * acc[i] * dt; // Nueva velocidad de la partícula
    }
}

void velocityVerlet2D(std::vector<Particle2D>& particles, double dt,
                        double k, double xmin, double xmax, double ymin, double ymax,
                        bool useBoundaries)
{
    int N = particles.size();
    std::vector<double> acc_x(N);
    std::vector<double> acc_y(N);

    computeAccelerations2D(particles, acc_x, acc_y, k);         //Calculamos aceleración en configuración "inicial"
    for(int i = 0; i < N; i++){
        //Movimiento en x
        particles[i].vx += 0.5 * acc_x[i] * dt;                 //Actualizamos velocidades medio paso
        particles[i].x += particles[i].vx * dt;                 //Corregimos la posición paso entero
        //Movimiento en y 
        particles[i].vy += 0.5 * acc_y[i] * dt;
        particles[i].y += particles[i].vy * dt;
    }

    if (useBoundaries){
    applyReflectiveBC2D(particles, xmin, xmax, ymin, ymax);     //Checkeamos condicion: agregar condiciones de frontera
    }
    
    computeAccelerations2D(particles, acc_x, acc_y, k);         //Calculamos la aceleración con la nueva configuración de posiciones

    for(int i = 0; i < N; i++){
        particles[i].vx += 0.5 * acc_x[i] * dt;                 //Nueva velocidad de la partícula
        particles[i].vy += 0.5 * acc_y[i] * dt;
    }
}

void applyReflectiveBC1D(std::vector<Particle1D>& particles,
                         double xmin,
                         double xmax)
{
    for (auto& p : particles) {

        // pared izquierda
        if (p.x < xmin) {
            p.x = 2.0 * xmin - p.x;  // reflejar posicion
            p.v = -p.v;              // invertir velocidad
        }

        // pared derecha
        if (p.x > xmax) {
            p.x = 2.0 * xmax - p.x;  // reflejar posicion
            p.v = -p.v;              // invertir velocidad
        }
    }
}

//aplica condiciones de frontera reflectivas 
void applyReflectiveBC2D(std::vector<Particle2D>& particles,
                         double xmin, double xmax,
                         double ymin, double ymax)
{
    for (auto& p : particles) {

        // pared izquierda
        if (p.x < xmin) {
            p.x = 2.0 * xmin - p.x;
            p.vx = -p.vx;
        }

        // pared derecha
        if (p.x > xmax) {
            p.x = 2.0 * xmax - p.x;
            p.vx = -p.vx;
        }

        // pared inferior
        if (p.y < ymin) {
            p.y = 2.0 * ymin - p.y;
            p.vy = -p.vy;
        }

        // pared superior
        if (p.y > ymax) {
            p.y = 2.0 * ymax - p.y;
            p.vy = -p.vy;
        }
    }
}