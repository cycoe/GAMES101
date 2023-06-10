#include "Material.hpp"
#include "Renderer.hpp"
#include "Scene.hpp"
#include "Triangle.hpp"
#include "Sphere.hpp"
#include "Vector.hpp"
#include "global.hpp"
#include "Material/Microfacet.hpp"
#include <chrono>

// In the main function of the program, we create the scene (create objects and
// lights) as well as set the options for the render (image width and height,
// maximum recursion depth, field-of-view, etc.). We then call the render
// function().
int main(int argc, char** argv)
{

    // Change the definition here to change resolution
    Scene scene(784, 784);

    Material* red = new Material(DIFFUSE, Vector3f(0.0f));
    red->Kd = Vector3f(0.63f, 0.065f, 0.05f);
    Material* green = new Material(DIFFUSE, Vector3f(0.0f));
    green->Kd = Vector3f(0.14f, 0.45f, 0.091f);
    Material* blue = new Material(DIFFUSE, Vector3f(0.0f));
    blue->Kd = Vector3f(0.091f, 0.32f, 0.65f);
    Material* white = new Material(DIFFUSE, Vector3f(0.0f));
    white->Kd = Vector3f(0.725f, 0.71f, 0.68f);
    Material* light = new Material(DIFFUSE, (12.0f * Vector3f(0.747f+0.058f, 0.747f+0.258f, 0.747f) + 20.6f * Vector3f(0.740f+0.287f,0.740f+0.160f,0.740f) + 25.f *Vector3f(0.737f+0.642f,0.737f+0.159f,0.737f)));
    light->Kd = Vector3f(0.65f);
    Material* mirror = new Material(SPECULAR, Vector3f(0.f));
    mirror->Ks = Vector3f(0.8f, 0.8f, 0.8f);
    Material* glass = new Material(GLASS, Vector3f(0.f));
    glass->Ks = Vector3f(0.98f);
    Material* glaze = new Material(GLAZE, Vector3f(0.f));
    glaze->Ks = Vector3f(0.98f);
    glaze->Kd = Vector3f(0.55f, 0.13f, 0.67f);
    Material* block = new Material(BLOCK, Vector3f(0.0f));
    block->Kd = Vector3f(0.725f, 0.71f, 0.68f);
    Material* mf = new Material(OREN_NAYAR, Vector3f(0.0f));
    mf->Kd = Vector3f(0.725f, 0.71f, 0.68f);
    Vector3f R(1.f, 1.f, 1.f);
    Material* bsdf = new BSDF(R, R, 0.1, 0.1);
    ///Vector3f R2(0.7f, 0.95f, 0.85f);
    ///Vector3f R2(1.f, 0.84f, 0.f);

    MeshTriangle floor("./models/cornellbox/floor.obj", block);
    MeshTriangle wall("./models/cornellbox/wall.obj", white);
    MeshTriangle shortbox("./models/cornellbox/shortbox.obj", white);
    MeshTriangle tallbox("./models/cornellbox/tallbox.obj", mirror);
    MeshTriangle left("./models/cornellbox/left.obj", red);
    MeshTriangle right("./models/cornellbox/right.obj", green);
    MeshTriangle light_("./models/cornellbox/light.obj", light);
    MeshTriangle bunny("./models/bunny/bunny.obj", glass);
    MeshTriangle spot("./models/spot/spot.obj", mirror);
    MeshTriangle pane("./models/cornellbox/glass.obj", glass);
    //MeshTriangle lucy1("./models/lucy1.obj", mfr);
    //MeshTriangle lucy2("./models/lucy2.obj", mf);
    //MeshTriangle lucy3("./models/lucy3.obj", glass);
    //MeshTriangle dragon("./models/xyzrgb_dragon.obj", bsdf);
    MeshTriangle teapot("./models/teapot.obj", bsdf);
    Sphere sphere1(Vector3f(350.f, 60.f, 300.f), 60.f, bsdf);
    Sphere sphere2(Vector3f(150.f, 60.f, 300.f), 60.f, mirror);

    scene.Add(&floor);
    scene.Add(&wall);
    //scene.Add(&shortbox);
    //scene.Add(&tallbox);
    scene.Add(&left);
    scene.Add(&right);
    scene.Add(&light_);
    //scene.Add(&lucy1);
    //scene.Add(&lucy2);
    //scene.Add(&lucy3);
    //scene.Add(&bunny);
    //scene.Add(&pane);
    //scene.Add(&spot);
    scene.Add(&sphere1);
    scene.Add(&sphere2);
    //scene.Add(&dragon);
    //scene.Add(&teapot);

    scene.buildBVH();

    Renderer r;
    int spp = 1024;
    if (argc > 1)
    {
        spp = std::atoi(argv[1]);
        spp = std::max(1, spp);
    }

    auto start = std::chrono::system_clock::now();
    r.Render(scene,spp);
    auto stop = std::chrono::system_clock::now();

    std::cout << "Render complete: \n";
    std::cout << "Time taken: " << std::chrono::duration_cast<std::chrono::hours>(stop - start).count() << " hours\n";
    std::cout << "          : " << std::chrono::duration_cast<std::chrono::minutes>(stop - start).count() << " minutes\n";
    std::cout << "          : " << std::chrono::duration_cast<std::chrono::seconds>(stop - start).count() << " seconds\n";

    return 0;
}
