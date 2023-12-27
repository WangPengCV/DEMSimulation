add_executable(SphereParticleTest SphereParticleTest.cpp)
target_link_libraries(SphereParticleTest Particle_lib DEM_lib)

add_executable( ParticlePropertiesTest ParticlePropertiesTest.cpp)
target_link_libraries(ParticlePropertiesTest ParticleProperties_lib DEM_lib)

add_executable(ParticleSimulationTest ParticleSimulationTest.cpp)
target_link_libraries(ParticleSimulationTest DEM_lib Visualization_lib ${VTK_LIBRARIES})

vtk_module_autoinit(
    TARGETS ParticleSimulationTest
    MODULES ${VTK_LIBRARIES}
  )

add_executable(VisualizeGridBasedContactDetection VisualizeGridBasedContactDetection.cpp)
target_link_libraries(VisualizeGridBasedContactDetection DEM_lib Visualization_lib ${VTK_LIBRARIES})

vtk_module_autoinit(
    TARGETS VisualizeGridBasedContactDetection
    MODULES ${VTK_LIBRARIES}
  )

  add_executable(TwoSphereImpact TwoSphereImpact.cpp)
  target_link_libraries(TwoSphereImpact DEM_lib Visualization_lib ${VTK_LIBRARIES})
  
  vtk_module_autoinit(
      TARGETS TwoSphereImpact
      MODULES ${VTK_LIBRARIES}
    )