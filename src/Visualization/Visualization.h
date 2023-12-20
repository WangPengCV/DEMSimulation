#pragma once

#include <vtkSmartPointer.h>
#include "DEMProperties.h"
class vtkPoints;
class vtkRenderer;
class vtkRenderWindow;
class vtkRenderWindowInteractor;
class vtkPolyData;
class vtkSphereSource;
class vtkActor;
class vtkPlaneSource;

class Visualization {
public:
    explicit Visualization(const DEMProperties& DEMproperties);
    void Update();

private:
    const DEMProperties& DEMproperties;

    //sphere
    vtkSmartPointer<vtkPoints> spherepoints;
    vtkSmartPointer<vtkPolyData> polyData;

    vtkSmartPointer<vtkRenderer> renderer;
    vtkSmartPointer<vtkRenderWindow> renderWindow;
    vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor;

    // Plane wall
    std::vector<vtkSmartPointer<vtkPlaneSource>> planeWallSource;
    std::vector<vtkSmartPointer<vtkActor>> planeWallActor;


    void UpdatePlaneWall();
    void UpdateSphere();
};
