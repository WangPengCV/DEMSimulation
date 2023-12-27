#include "LineSphereSource.h"


vtkStandardNewMacro(LineSphereSource);

LineSphereSource::LineSphereSource() : Radius(1.0) {
    this->SetNumberOfInputPorts(0);
}

int LineSphereSource::RequestData(vtkInformation*, vtkInformationVector**, vtkInformationVector* outputVector) {
    vtkInformation *outInfo = outputVector->GetInformationObject(0);
    vtkPolyData *output = vtkPolyData::SafeDownCast(outInfo->Get(vtkDataObject::DATA_OBJECT()));

    vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
    vtkSmartPointer<vtkCellArray> polys = vtkSmartPointer<vtkCellArray>::New();
    vtkSmartPointer<vtkCellArray> lines = vtkSmartPointer<vtkCellArray>::New();

    

    vtkSmartPointer<vtkUnsignedCharArray> colors = vtkSmartPointer<vtkUnsignedCharArray>::New();
    colors->SetNumberOfComponents(3); // RGB
    colors->SetName("Colors");

    CreateSphere(points, polys);
    AddLineOnSphere(points, lines, colors);

    output->SetPoints(points);
    output->SetPolys(polys);
    output->SetLines(lines);
    output->GetPointData()->SetScalars(colors); // Attach the color array to the output


    return 1;
}

void LineSphereSource::CreateSphere(vtkPoints* points, vtkCellArray* polys) {
   // Basic sphere creation using vtkSphereSource as a reference
    vtkSmartPointer<vtkSphereSource> sphereSrc = vtkSmartPointer<vtkSphereSource>::New();
    sphereSrc->SetRadius(this->Radius);
    sphereSrc->SetThetaResolution(100);
    sphereSrc->SetPhiResolution(100);
    sphereSrc->Update();

    vtkPolyData* sphereData = sphereSrc->GetOutput();
    points->DeepCopy(sphereData->GetPoints());
    polys->DeepCopy(sphereData->GetPolys());
}

void LineSphereSource::AddLineOnSphere(vtkPoints* points, vtkCellArray* lines,vtkUnsignedCharArray* colors) {
    const int numLinePoints = 50;
    std::vector<vtkIdType> linePointIds(numLinePoints);

    for (int i = 0; i < numLinePoints; ++i) {
        double theta = vtkMath::Pi() * i / (numLinePoints - 1);
        double x = this->Radius * sin(theta);
        double y = 0.0;
        double z = this->Radius * cos(theta);

        linePointIds[i] = points->InsertNextPoint(x, y, z);
        unsigned char red[3] = {255, 0, 0};
        colors->InsertNextTypedTuple(red);
    }

    vtkSmartPointer<vtkPolyLine> polyLine = vtkSmartPointer<vtkPolyLine>::New();
    polyLine->GetPointIds()->SetNumberOfIds(numLinePoints);
    for (int i = 0; i < numLinePoints; ++i) {
        polyLine->GetPointIds()->SetId(i, linePointIds[i]);
    }

    lines->InsertNextCell(polyLine);
}
// LineSphereSource.cpp
LineSphereSource::~LineSphereSource() {
    // Destructor code here (if needed)
}

