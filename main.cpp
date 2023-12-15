// // main.cpp
// #include "Simulation.h"
// #include "Visualization.h"
// int main() {
//     double r = 0.002*1.1;
//     double space_r = 12 * r;

//     int numSpheres = 352;

//     double V = numSpheres * (1.3333 * 3.14 * r * r * r);

//     double space_V = (1.3333 * 3.14 * space_r * space_r * space_r);

//     double voulumfraction = 0.4;

//     Simulation simulation(numSpheres,r,space_r,100000000,(1.3333 * 3.14 * r * r * r)*2130,0.9);

    
//     Visualization visualization(simulation);

//     for (int step = 0; step < 1000000; ++step) {
//         simulation.Update();
//         visualization.Update();
//     }

//     return 0;
// }
