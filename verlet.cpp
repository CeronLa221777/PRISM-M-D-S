#include "verlet.hpp"
#include <cmath>

constexpr double RCUT2 = 16.0; //radio de cort definido arriba en el codigo

void computeAccelerations1D(const std::vector<Particle1D>& particles,
                          std::vector<double>& acc, double k) //Nos toca calcular por aparte las aceleraciones para no todo meterlo en el verlet
{
    int N = particles.size();
    double rcut2 = RCUT2;
    std::fill(acc.begin(), acc.end(), 0.0);

    // Trampa armónica
    for(int i = 0; i < N; i++)
        acc[i] = -k*particles[i].x;

    // Interacción por pares
    for(int i = 0; i < N; i++){
        for(int j = i+1; j < N; j++){

            double dx = particles[i].x - particles[j].x; // Vector relativo
            double r2 = dx*dx; // Norma cuadrado del vector relativo
            
            if (r2 < rcut2){
                double r4 = r2 * r2;
                double r6 = r4 * r2;
                double r14 = r6 * r6 *r2;
                double f = 12.0 / r14;
                double f_scalar = f * dx; // spft-core luego de la derivada da elevado a la 14
                acc[i] += f_scalar;// fuerza que ejerce j sobre i
                acc[j] -= f_scalar;// línea que hace que no tengamos cálculos duplicados, fuerza que recibe j de i por tecera ley de newton
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
    double rcut2 = RCUT2;

    std::fill(acc_x.begin(), acc_x.end(), 0.0);
    std::fill(acc_y.begin(), acc_y.end(), 0.0);

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
                double r4 = r2 * r2;
                double r6 = r4 * r2;
                double r14 = r6 * r6 *r2;
                double f_scalar = 12.0  / r14; // soft-core luego de la derivada da elevado a la 14
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

void computeAccelerations3D(const std::vector<Particle3D>& particles,
                            std::vector<double>& acc_x,
                            std::vector<double>& acc_y,
                            std::vector<double>& acc_z,
                            double k) //Nos toca calcular por aparte las aceleraciones para no todo meterlo en el verlet
{

    int N = particles.size();
    double rcut2 = RCUT2;
    // Trampa armónica

    std::fill(acc_x.begin(), acc_x.end(), 0.0);
    std::fill(acc_y.begin(), acc_y.end(), 0.0);
    std::fill(acc_z.begin(), acc_z.end(), 0.0);

    for(int i = 0; i < N; i++){
        acc_x[i] = -k*particles[i].x;
        acc_y[i] = -k*particles[i].y;
        acc_z[i] = -k*particles[i].z;
    }
    // Interacción por pares (soft-core)
    for(int i = 0; i < N; i++){
        for(int j = i+1; j < N; j++){

            double dx = particles[i].x - particles[j].x; // Vector relativo
            double dy = particles[i].y - particles[j].y; 
            double dz = particles[i].z - particles[j].z;
            double r2 = dx*dx + dy*dy + dz*dz; // Norma cuadrado del vector relativo
            
            if (r2 < rcut2){
                double r4 = r2 * r2;
                double r6 = r4 * r2;
                double r14 = r6 * r6 *r2;
                double f_scalar = 12.0  / r14; // soft-core luego de la derivada da elevado a la 14
                                
                acc_x[i] += f_scalar*dx;// fuerza que ejerce j sobre i
                acc_y[i] += f_scalar*dy;
                acc_z[i] += f_scalar*dz;

                acc_x[j] -= f_scalar*dx;// línea que hace que no tengamos cálculos duplicados, fuerza que recibe j de i por tecera ley de newton
                acc_y[j] -= f_scalar*dy;
                acc_z[j] -= f_scalar*dz;
            }
        }
    }
}

//implementacion velocity verlet para las diferentes dimensiones
void velocityVerlet1D(std::vector<Particle1D>& particles, double dt,
                        double k, double xmin, double xmax,
                        bool useBoundaries)
{
    int N = particles.size();
    std::vector<double> acc(N, 0.0);

    computeAccelerations1D(particles, acc, k); // Calculamos aceleración en configuración "inicial"
    for(int i = 0; i < N; i++){
        particles[i].v += 0.5 * acc[i] * dt;// afectamos la velocidad de la partícula
        particles[i].x += particles[i].v * dt;// corregimos la posición 
    }

    if (useBoundaries){
        applyReflectiveBC1D(particles, xmin, xmax);     //Checkeamos condicion: agregar condiciones de frontera
    }

    computeAccelerations1D(particles, acc, k); // Calculamos la aceleración con la nueva configuración de posiciones

    for(int i = 0; i < N; i++){
        particles[i].v += 0.5 * acc[i] * dt; // Nueva velocidad de la partícula
    }
}

void velocityVerlet2D(std::vector<Particle2D>& particles, double dt,
                        double k, double xmin, double xmax, double ymin, double ymax,
                        bool useBoundaries)
{
    int N = particles.size();
    std::vector<double> acc_x(N, 0.0);
    std::vector<double> acc_y(N, 0.0);

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

void velocityVerlet3D(std::vector<Particle3D>& particles, double dt,
                        double k, double xmin, double xmax,
                                  double ymin, double ymax,
                                  double zmin, double zmax,
                                bool useBoundaries)
{
    int N = particles.size();
    std::vector<double> acc_x(N, 0.0);
    std::vector<double> acc_y(N, 0.0);
    std::vector<double> acc_z(N, 0.0);

    computeAccelerations3D(particles, acc_x, acc_y, acc_z, k);         //Calculamos aceleración en configuración "inicial"
    for(int i = 0; i < N; i++){
        //velocidades en x,y,z
        particles[i].vx += 0.5 * acc_x[i] * dt;                 //Actualizamos velocidades medio paso
        particles[i].vy += 0.5 * acc_y[i] * dt;
        particles[i].vz += 0.5 * acc_z[i] * dt;
        //Movimiento en x,y,z 
        particles[i].x += particles[i].vx * dt;                 //Corregimos la posición paso entero
        particles[i].y += particles[i].vy * dt;
        particles[i].z += particles[i].vz * dt;

    }

    if (useBoundaries){
    applyReflectiveBC3D(particles, xmin, xmax, ymin, ymax, zmin, zmax);     //Checkeamos condicion: agregar condiciones de frontera
    }
    
    computeAccelerations3D(particles, acc_x, acc_y, acc_z, k);         //Calculamos la aceleración con la nueva configuración de posiciones
    for(int i = 0; i < N; i++){
        particles[i].vx += 0.5 * acc_x[i] * dt;                 //Nueva velocidad de la partícula
        particles[i].vy += 0.5 * acc_y[i] * dt;
        particles[i].vz += 0.5 * acc_z[i] * dt;
    }
}

//aplicar las condiciones de frontera
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

void applyReflectiveBC3D(std::vector<Particle3D>& particles,
                         double xmin, double xmax,
                         double ymin, double ymax,
                         double zmin, double zmax)
{
    for (auto& p : particles) {

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