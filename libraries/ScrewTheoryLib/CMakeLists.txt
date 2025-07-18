if(NOT orocos_kdl_FOUND AND (NOT DEFINED ENABLE_ScrewTheoryLib OR ENABLE_ScrewTheoryLib))
    message(WARNING "orocos_kdl package not found, disabling ScrewTheoryLib")
endif()

cmake_dependent_option(ENABLE_ScrewTheoryLib "Enable/disable ScrewTheoryLib library" ON
                       orocos_kdl_FOUND OFF)

if(ENABLE_ScrewTheoryLib)

    add_library(ScrewTheoryLib SHARED ScrewTheoryTools.hpp
                                      ScrewTheoryTools.cpp
                                      MatrixExponential.hpp
                                      MatrixExponential.cpp
                                      ProductOfExponentials.hpp
                                      ProductOfExponentials.cpp
                                      ScrewTheoryIkProblem.hpp
                                      ScrewTheoryIkProblem.cpp
                                      ScrewTheoryIkProblemBuilder.cpp
                                      ScrewTheoryIkSubproblems.hpp
                                      PadenKahanSubproblems.cpp
                                      PardosGotorSubproblems.cpp
                                      ConfigurationSelector.hpp
                                      ConfigurationSelector.cpp
                                      ConfigurationSelectorLeastOverallAngularDisplacement.cpp
                                      ConfigurationSelectorHumanoidGait.cpp)

    set_property(TARGET ScrewTheoryLib PROPERTY PUBLIC_HEADER MatrixExponential.hpp
                                                              ProductOfExponentials.hpp
                                                              ScrewTheoryIkProblem.hpp
                                                              ConfigurationSelector.hpp)

    target_link_libraries(ScrewTheoryLib PUBLIC ${orocos_kdl_LIBRARIES})

    target_include_directories(ScrewTheoryLib PUBLIC ${orocos_kdl_INCLUDE_DIRS}
                                                     $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
                                                     $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)

    target_compile_features(ScrewTheoryLib PUBLIC cxx_std_17) # nested namespace definition

    install(TARGETS ScrewTheoryLib
            EXPORT SCREW_THEORY_SOLVERS)

    set_property(GLOBAL APPEND PROPERTY _exported_dependencies orocos_kdl)

    add_library(ROBOTICSLAB::ScrewTheoryLib ALIAS ScrewTheoryLib)

else()

    set(ENABLE_ScrewTheoryLib OFF CACHE BOOL "Enable/disable ScrewTheoryLib library" FORCE)

endif()
