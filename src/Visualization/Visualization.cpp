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
#include <vtkCubeSource.h>
#include <vtkFloatArray.h>
#include <vtkUnsignedCharArray.h>
#include <vtkPointData.h>
#include <vtkLight.h>
#include <vtkCamera.h>
#include <vtkPlaneSource.h>

Visualization::Visualization(const DEMProperties& DEMproperties) : DEMproperties(DEMproperties) 
{
    spherepoints = vtkSmartPointer<vtkPoints>::New();
    UpdateSphere();

    vtkSmartPointer<vtkSphereSource> sphereSource = vtkSmartPointer<vtkSphereSource>::New();
    sphereSource->SetPhiResolution(20);    // Increase vertical resolution
    sphereSource->SetThetaResolution(20);  // Increase horizontal resolution

    vtkSmartPointer<vtkGlyph3D> glyph3D = vtkSmartPointer<vtkGlyph3D>::New();
    glyph3D->SetSourceConnection(sphereSource->GetOutputPort());
    glyph3D->SetInputData(polyData);
    glyph3D->SetScaleModeToScaleByScalar();
    glyph3D->SetInputArrayToProcess(0, 0, 0, vtkDataObject::FIELD_ASSOCIATION_POINTS, "scales");
    glyph3D->SetColorModeToColorByScalar();
    glyph3D->SetInputArrayToProcess(1, 0, 0, vtkDataObject::FIELD_ASSOCIATION_POINTS, "colors");
    glyph3D->Update();

    vtkSmartPointer<vtkPolyDataMapper> mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    mapper->SetInputConnection(glyph3D->GetOutputPort());

    vtkSmartPointer<vtkActor> actor = vtkSmartPointer<vtkActor>::New();
    actor->SetMapper(mapper);



   

    renderer = vtkSmartPointer<vtkRenderer>::New();
    renderWindow = vtkSmartPointer<vtkRenderWindow>::New();
    renderWindow->AddRenderer(renderer);
    renderer->AddActor(actor);

     for (const auto& wall : DEMproperties.getPlaneWall()) {
        // Create a plane source
        vtkSmartPointer<vtkPlaneSource> planesource = vtkSmartPointer<vtkPlaneSource>::New();
        planesource->SetOrigin(wall->getCorner2().x(),wall->getCorner2().y(),wall->getCorner2().z()); 
        planesource->SetPoint1(wall->getCorner1().x(),wall->getCorner1().y(),wall->getCorner1().z()); 
        planesource->SetPoint2(wall->getCorner3().x(),wall->getCorner3().y(),wall->getCorner3().z());

       
        // Map the transformed plane to a mapper
        vtkSmartPointer<vtkPolyDataMapper> planemapper = vtkSmartPointer<vtkPolyDataMapper>::New();
        planemapper->SetInputConnection(planesource->GetOutputPort());

        // Create an actor to represent the wall
        vtkSmartPointer<vtkActor> planeactor = vtkSmartPointer<vtkActor>::New();
        planeactor->SetMapper(planemapper);
        // Set properties like color, transparency, etc.
        planeactor->GetProperty()->SetColor(1.0, 1.0, 1.0); // Red color
        planeactor->GetProperty()->EdgeVisibilityOn();
        planeactor->GetProperty()->SetOpacity(0.1);

        // Add the actor to the renderer
        renderer->AddActor(planeactor);

        // Store the sources and actors for later updates
        planeWallSource.push_back(planesource);
        planeWallActor.push_back(planeactor);
    }
   
    renderer->SetBackground(1, 1, 1);

    renderWindowInteractor = vtkSmartPointer<vtkRenderWindowInteractor>::New();
    renderWindowInteractor->SetRenderWindow(renderWindow);

}

void Visualization::Update() {
    UpdateSphere();
    UpdatePlaneWall();
    renderWindowInteractor->GetRenderWindow()->Render();
    if (renderWindowInteractor) {
        renderWindowInteractor->ProcessEvents();
    }
}

void Visualization::UpdateSphere() {
    spherepoints->Reset();
    // Create scale and color arrays
    vtkSmartPointer<vtkFloatArray> scales = vtkSmartPointer<vtkFloatArray>::New();
    scales->SetName("scales");
    vtkSmartPointer<vtkUnsignedCharArray> colors = vtkSmartPointer<vtkUnsignedCharArray>::New();
    colors->SetName("colors");
    colors->SetNumberOfComponents(3);  // 3 components for R, G, B
    for (const auto& particle : DEMproperties.getParticles()) 
    {
        auto manager = particle->getParticlePropertyManager();
        auto typemapping = manager->gettypeMapping();

        if(typemapping[particle->getType().getCategory()] == ParticleType::SPHERE)
        {
            spherepoints->InsertNextPoint(particle->getPosition().x(),particle->getPosition().y(),particle->getPosition().z());
            float scale = 2*manager->getSphereProperties(particle->getType())->getRadius();
            scales->InsertNextValue(scale);
            unsigned char color[3] = {0,0,255};
            colors->InsertNextTypedTuple(color);
        }
    }
    spherepoints->Modified();

    if (!polyData) {
        polyData = vtkSmartPointer<vtkPolyData>::New();
        polyData->SetPoints(spherepoints);
    }
    // Add scale and color arrays to polyData
    polyData->GetPointData()->SetScalars(colors);
    polyData->GetPointData()->AddArray(scales);
}

void Visualization::UpdatePlaneWall()
{
    int plane_index = 0;
     for (const auto& wall : DEMproperties.getPlaneWall()) {
       
        planeWallSource[plane_index]->SetOrigin(wall->getCorner2().x(),wall->getCorner2().y(),wall->getCorner2().z()); 
        planeWallSource[plane_index]->SetPoint1(wall->getCorner1().x(),wall->getCorner1().y(),wall->getCorner1().z()); 
        planeWallSource[plane_index]->SetPoint2(wall->getCorner3().x(),wall->getCorner3().y(),wall->getCorner3().z());
        plane_index++;
     }

}
