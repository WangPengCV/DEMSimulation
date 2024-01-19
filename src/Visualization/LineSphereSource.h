#pragma once
#include <vtkPolyDataAlgorithm.h>
#include <vtkSmartPointer.h>
#include <vtkPoints.h>
#include <vtkCellArray.h>
#include <vtkSphereSource.h>
#include <vtkPolyLine.h>
#include <vtkMath.h>
#include <vtkObjectFactory.h>
#include <vtkInformation.h>
#include <vtkInformationVector.h>
#include <vtkPolyData.h>
#include <vtkPointData.h>
#include <vtkUnsignedCharArray.h>

class LineSphereSource : public vtkPolyDataAlgorithm {
public:
    vtkTypeMacro(LineSphereSource, vtkPolyDataAlgorithm);
    static LineSphereSource* New();

    vtkSetMacro(Radius, double);
    vtkGetMacro(Radius, double);

protected:
    LineSphereSource();
    ~LineSphereSource() override;

    int RequestData(vtkInformation*, vtkInformationVector**, vtkInformationVector*) override;

private:
    double Radius;

    void CreateSphere(vtkPoints* points, vtkCellArray* polys);
    void AddLineOnSphere(vtkPoints* points, vtkCellArray* lines,vtkUnsignedCharArray* colors);
};

