// Visualization.cpp
#include "Visualization.h"
#include <vtkCylinderSource.h>
#include <vtkSphereSource.h>
#include <vtkPolyDataMapper.h>
#include <vtkActor.h>
#include <vtkGlyph3D.h>
#include <vtkTransform.h>
#include <vtkTransformPolyDataFilter.h>
#include <vtkPolyData.h>
#include <vtkSmartPointer.h>
#include <vtkPoints.h>
#include <vtkRenderer.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkPolyData.h>
#include <vtkProperty.h>

Visualization::Visualization(const Simulation& simulation) : simulation(simulation) {
    points = vtkSmartPointer<vtkPoints>::New();
    UpdatePoints();

    vtkSmartPointer<vtkSphereSource> sphereSource = vtkSmartPointer<vtkSphereSource>::New();
    sphereSource->SetRadius(0.002);
    sphereSource->SetPhiResolution(20);    // Increase vertical resolution
    sphereSource->SetThetaResolution(20);  // Increase horizontal resolution

    vtkSmartPointer<vtkGlyph3D> glyph3D = vtkSmartPointer<vtkGlyph3D>::New();
    glyph3D->SetSourceConnection(sphereSource->GetOutputPort());
    glyph3D->SetInputData(polyData);
    glyph3D->Update();

    vtkSmartPointer<vtkPolyDataMapper> mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    mapper->SetInputConnection(glyph3D->GetOutputPort());

    vtkSmartPointer<vtkActor> actor = vtkSmartPointer<vtkActor>::New();
    actor->SetMapper(mapper);


    // Initialize big sphere
    bigSphereSource = vtkSmartPointer<vtkSphereSource>::New();
    bigSphereSource->SetRadius(simulation.GetSphereWall().radius);
    bigSphereSource->SetCenter(simulation.GetSphereWall().x,simulation.GetSphereWall().y,simulation.GetSphereWall().z);
    bigSphereSource->SetPhiResolution(30);    // Increase vertical resolution
    bigSphereSource->SetThetaResolution(30);  // Increase horizontal resolution
    vtkSmartPointer<vtkPolyDataMapper> bigSphereMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    bigSphereMapper->SetInputConnection(bigSphereSource->GetOutputPort());

    bigSphereActor = vtkSmartPointer<vtkActor>::New();
    bigSphereActor->SetMapper(bigSphereMapper);
    bigSphereActor->GetProperty()->SetOpacity(0.7); // Adjust transparency here



    renderer = vtkSmartPointer<vtkRenderer>::New();
    renderWindow = vtkSmartPointer<vtkRenderWindow>::New();
    renderWindow->AddRenderer(renderer);
    renderer->AddActor(actor);
    renderer->AddActor(bigSphereActor);

    renderer->SetBackground(0.1, 0.2, 0.3);

    renderWindowInteractor = vtkSmartPointer<vtkRenderWindowInteractor>::New();
    renderWindowInteractor->SetRenderWindow(renderWindow);
}

void Visualization::Update() {
    UpdatePoints();
    UpdateSphere();
    renderWindowInteractor->GetRenderWindow()->Render();
    if (renderWindowInteractor) {
        renderWindowInteractor->ProcessEvents();
    }
}

void Visualization::UpdatePoints() {
    points->Reset();
    for (const auto& sphere : simulation.GetSpheres()) {
        points->InsertNextPoint(sphere.x, sphere.y, sphere.z);
    }
    points->Modified();

    if (!polyData) {
        polyData = vtkSmartPointer<vtkPolyData>::New();
        polyData->SetPoints(points);
    }
}

void Visualization::UpdateSphere()
{
    // Update the big sphere (if needed)
    bigSphereSource->SetRadius(simulation.GetSphereWall().radius);
    bigSphereSource->SetCenter(simulation.GetSphereWall().x,simulation.GetSphereWall().y,simulation.GetSphereWall().z);
    bigSphereSource->Update();

}
