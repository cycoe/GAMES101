set(CMAKE_CXX_FLAGS "-O3")
set(CMAKE_CXX_STANDARD 17)

add_executable(07 main.cpp Object.hpp Vector.cpp Vector.hpp Sphere.hpp global.hpp Triangle.hpp Scene.cpp
        Scene.hpp Light.hpp AreaLight.hpp BVH.cpp BVH.hpp Bounds3.hpp Ray.hpp Material.hpp Intersection.hpp
        Renderer.cpp Renderer.hpp)

target_link_libraries(07 pthread)
