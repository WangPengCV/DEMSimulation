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
<<<<<<< HEAD
class vtkCellArray;
class vtkGlyph3DMapper;

class Visualization
{
public:
    explicit Visualization(const DEMProperties &DEMproperties);
    void Update();

private:
    const DEMProperties &DEMproperties;

    // sphere
    vtkSmartPointer<vtkPoints> spherepoints;
    vtkSmartPointer<vtkPolyData> spherepolyData;

    // cylinder
    vtkSmartPointer<vtkPoints> cylinderpoints;
    vtkSmartPointer<vtkPolyData> cylinderpolydata;
    //fibersphere
    vtkSmartPointer<vtkPoints> fiberspherepoints;
    vtkSmartPointer<vtkPolyData> fiberspherepolyData;
=======

class Visualization {
public:
    explicit Visualization(const DEMProperties& DEMproperties);
    void Update();

private:
    const DEMProperties& DEMproperties;

    //sphere
    vtkSmartPointer<vtkPoints> spherepoints;
    vtkSmartPointer<vtkPolyData> polyData;
>>>>>>> 686cbfa3ebadc1d4aba7bce443978911f7964200

    vtkSmartPointer<vtkRenderer> renderer;
    vtkSmartPointer<vtkRenderWindow> renderWindow;
    vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor;

    // Plane wall
    std::vector<vtkSmartPointer<vtkPlaneSource>> planeWallSource;
    std::vector<vtkSmartPointer<vtkActor>> planeWallActor;

<<<<<<< HEAD
    void UpdatePlaneWall();
    void UpdateSphere();
    void UpdataFibers();

    int count_numer;
=======

    void UpdatePlaneWall();
    void UpdateSphere();
>>>>>>> 686cbfa3ebadc1d4aba7bce443978911f7964200
};
