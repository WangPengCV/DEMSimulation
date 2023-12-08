#pragma once

#include "Simulation.h"
#include <vtkSmartPointer.h>

class vtkPoints;
class vtkRenderer;
class vtkRenderWindow;
class vtkRenderWindowInteractor;
class vtkPolyData;
class vtkSphereSource;
class vtkActor;

class Visualization {
public:
    explicit Visualization(const Simulation& simulation);
    void Update();

private:
    const Simulation& simulation;
    vtkSmartPointer<vtkPoints> points;
    vtkSmartPointer<vtkRenderer> renderer;
    vtkSmartPointer<vtkRenderWindow> renderWindow;
    vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor;
    vtkSmartPointer<vtkPolyData> polyData;

    // Big sphere
    vtkSmartPointer<vtkSphereSource> bigSphereSource;
    vtkSmartPointer<vtkActor> bigSphereActor;


    void UpdatePoints();
    void UpdateSphere();
};
