#ifndef PARTICLE_HPP
#define PARTICLE_HPP

struct Particle1D {
    double x; //guarda la posición de la partícula
    double v; //guarda la velocidad de la particula en una dimension
};// El archivo presente no tiene algoritmo porque este es para guardar solo info del sistema, No guarda masa porque m=1

struct Particle2D {
    double x, y; //guarda la posición de la partícula
    double vx, vy;//guardamos posiciónde partícula
};

struct Particle3D {
    double x, y, z;
    double vx, vy, vz;
};


#endif