#include "bounding_volume_hierarchy.h"
#include "disable_all_warnings.h"
#include "draw.h"
#include "image.h"
#include "ray_tracing.h"
#include "screen.h"
#include "trackball.h"
#include "window.h"
// Disable compiler warnings in third-party code (which we cannot change).
DISABLE_WARNINGS_PUSH()
#include <glm/gtc/type_ptr.hpp>
#include <glm/vec2.hpp>
#include <glm/vec4.hpp>
#include <imgui.h>
DISABLE_WARNINGS_POP()
#include <chrono>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <optional>
#include <random>
#include <string>
#include <type_traits>
#ifdef USE_OPENMP
#include <omp.h>
#endif

// This is the main application. The code in here does not need to be modified.
constexpr glm::ivec2 windowResolution{ 1200, 960 };
const std::filesystem::path dataPath{ DATA_DIR };
const std::filesystem::path outputPath{ OUTPUT_DIR };

std::random_device rd;                                  // NOTE: (Jean-Paul) These imports are essential for random point generation
std::mt19937 gen(rd());                                 // https://stackoverflow.com/questions/9878965/rand-between-0-and-1
std::uniform_real_distribution<> dis(-1, 1);
glm::vec3 pixels[windowResolution.y][windowResolution.x];
glm::vec3 bpixels[windowResolution.y][windowResolution.x];

int glossLevel = 8;
int sampleLevel = 8;
int bloomFilterSize = 3;

bool enableAi{ false };
bool enableBloom{ false };
bool renderingToFile{ false };
enum class ViewMode {
    Rasterization = 0,
    RayTracing = 1
};

// Random vector generation happens at the beginning to
//   distribute the spherical lights
std::vector<glm::vec3> randomVectors;
void initialize(int samples) {
    for (int i = 0; i < samples; i++) {
        float xr = dis(gen);
        float yr = dis(gen);
        float zr = dis(gen);

        glm::vec3 randomVector = 0.9f * glm::vec3(xr, yr, zr);
        randomVectors.push_back(randomVector);
    }
}

// Hard Shadows
// tests if any object lies between hitpoint and lightsource
// returns true if position is in a shadow, false if facing light
bool castShadow(glm::vec3 vertexPos, glm::vec3 lightsource, BoundingVolumeHierarchy bvh) {
    float distance = glm::length(vertexPos - lightsource);
    glm::vec3 destination = lightsource;
    glm::vec3 origin = vertexPos;
    glm::vec3 direction = glm::normalize(destination - origin);

    Ray shadowray;
    shadowray.origin = origin+direction*0.00001f;
    shadowray.direction = direction;
    shadowray.t = INFINITY;

    // now here we shoot our ray but nobody knows how ???
    HitInfo h;
    bool intersects = bvh.intersect(shadowray, h);
    if (!(distance < shadowray.t)) {
        drawRay(shadowray, glm::vec3{ 1.0f, 1.0f ,1.0f }); // draw a red debug ray
        return true;
    }
    drawRay(shadowray, glm::vec3{ 0.0f, 1.0f ,0 }); // draw a green debug ray if not in a shadow
    return false;
}

// Soft Shadows
// Creates a spherical light source with center and radius.
// returns the percentage of shadow coverage
float softShadow(glm::vec3 intersectionPos, SphericalLight sphericalLight, const BoundingVolumeHierarchy& bvh) {
    glm::vec3 center = sphericalLight.position;                             // center of our light sphere
    float r = sphericalLight.radius;

    if (r == 0) {                                                       // if pointed to the light
        if (castShadow(intersectionPos, center, bvh)) return 1;              // If in shadow skip further computation for this lightsource
    }
    // we generate random points as point lights
    float contributions = 0;
    for (int i = 0; i < sampleLevel; i++) {
        if (castShadow(intersectionPos, center + randomVectors[i] * r, bvh)) {
            contributions++;
        }
    }

    float coverage = (contributions / sampleLevel);  // calculate percentage of the light coverage
    return coverage;
}

static glm::vec3 getFaceColour(const Scene& scene, const BoundingVolumeHierarchy& bvh, Ray ray, HitInfo hitInfo) {
    glm::vec3 intersectPos = ray.origin + ray.direction * ray.t;
    glm::vec3 kd = hitInfo.material.kd;
    float lightIntensity = 1.0f;

    glm::vec3 diffuse{ 0 };
    glm::vec3 spec{ 0 };
    for (SphericalLight const slight : scene.sphericalLight) {
        lightIntensity -= softShadow(intersectPos, slight, bvh);
        glm::vec3 lightVec = glm::normalize(slight.position - intersectPos);
        glm::vec3 normal = glm::normalize(hitInfo.normal);

        lightIntensity = lightIntensity / scene.sphericalLight.size();

        float dot = glm::dot(normal, lightVec);
        if (dot > 0) {
            diffuse = diffuse + (kd * slight.color * dot) * lightIntensity;
        }

        glm::vec3 viewVec = glm::normalize(intersectPos - ray.origin);
        glm::vec3 reflectionVec = glm::normalize(glm::reflect(lightVec, hitInfo.normal));

        float dotprod = glm::dot(viewVec, reflectionVec);
        if (dotprod > 0) {
            spec = spec + (hitInfo.material.ks * slight.color * pow(dotprod, hitInfo.material.shininess)) * lightIntensity;
        }
    }

    for (PointLight const light : scene.pointLights) {
        if (!castShadow(intersectPos, light.position, bvh)) { lightIntensity = 1; }
        glm::vec3 lightVec = glm::normalize(light.position - intersectPos);
        glm::vec3 normal = glm::normalize(hitInfo.normal);

        lightIntensity = lightIntensity / scene.pointLights.size();
        float dot = glm::dot(normal, lightVec);
        if (dot > 0) {
            diffuse = diffuse + (kd * light.color * dot) * lightIntensity;
        }

        glm::vec3 viewVec = glm::normalize(intersectPos - ray.origin);
        glm::vec3 reflectionVec = glm::normalize(glm::reflect(lightVec, hitInfo.normal));

        float dotprod = glm::dot(viewVec, reflectionVec);
        if (dotprod > 0) {
            spec = spec + (hitInfo.material.ks * light.color * pow(dotprod, hitInfo.material.shininess)) * lightIntensity;
        }

    }
    Ray rey{ intersectPos, hitInfo.normal, 1 };
    glm::vec3 res = diffuse + spec;
    //drawRay(ray, res);
    return res;
}

void drawSquare(double x1, double y1, double sidelength)
{
    double halfside = sidelength / 2;

    glColor3d(0, 0, 0);
    glBegin(GL_POLYGON);

    glVertex2d(x1 + halfside, y1 + halfside);
    glVertex2d(x1 + halfside, y1 - halfside);
    glVertex2d(x1 - halfside, y1 - halfside);
    glVertex2d(x1 - halfside, y1 + halfside);

    glEnd();
}

static glm::vec3 recursiveRayTracing(const Scene& scene, const BoundingVolumeHierarchy& bvh, Ray ray, HitInfo hitInfo, int level, int maxLevel, glm::vec3 hitColor, bool glossy) {
    Ray reflect = calculateReflectionRay(ray, hitInfo);

    if (glossy) {
        float dot = glm::dot(hitInfo.normal, (reflect.direction - reflect.origin));
        float angle = glm::degrees(glm::acos(dot));

        //glm::rotate(, reflect.direction, angle);
        float samples = (float)glossLevel;
        glm::vec3 basedirection = reflect.direction;
        glm::vec3 finalColor = glm::vec3{ 0.0f };
        int hit = 0;
        for (int i = 0; i < samples; i++) {
            reflect.direction = basedirection + (randomVectors[i] * 0.04f);
            reflect.t = 1000.0f;
            if (bvh.intersect(reflect, hitInfo)) {
                hit++;
                drawRay(reflect);
                if (hitInfo.material.ks != glm::vec3{ 0.0f } && hitInfo.material.ks.x < 0.9f)
                {
                    finalColor+= recursiveRayTracing(scene, bvh, reflect, hitInfo, level + 1, maxLevel, hitColor + getFaceColour(scene, bvh, ray, hitInfo), true);
                }
                else if (hitInfo.material.ks != glm::vec3{ 0.0f } && hitInfo.material.ks.x > 0.9f)
                {
                    finalColor += recursiveRayTracing(scene, bvh, reflect, hitInfo, level + 1, maxLevel, hitColor + getFaceColour(scene, bvh, ray, hitInfo), false);
                }
                else {
                    finalColor += getFaceColour(scene, bvh, ray, hitInfo);
                }
            }
         
        }
        finalColor = glm::vec3{ finalColor.x / hit, finalColor.y / hit, finalColor.z / hit };
        return finalColor;
      
    }
    else {
        Ray prev = reflect;
        if (bvh.intersect(reflect, hitInfo) && level < maxLevel) {
            drawRay(reflect);
            //Perfect mirror

            if (hitInfo.material.ks != glm::vec3{ 0 } && hitInfo.material.ks.x > 0.9) {
                return recursiveRayTracing(scene, bvh, reflect, hitInfo, level + 1, maxLevel, hitColor + getFaceColour(scene, bvh, ray, hitInfo), false);
            }

        }
        else {
            drawRay(prev, glm::vec3{ 1.0f, 0 ,0 });
            return hitColor;
        }
    }
    
/*
Ray prev = reflect;
    HitInfo prevHit = hitInfo;
    if (bvh.intersect(reflect, hitInfo) && level < maxLevel) {
        glm::vec3 resColor = hitColor + getFaceColour(scene, bvh, ray, hitInfo);
        drawRay(reflect, resColor);
        if (hitInfo.material.ks != glm::vec3{ 0 }) {
            return recursiveRayTracing(scene, bvh, reflect, hitInfo, level + 1, maxLevel, resColor);



*/
    return hitColor + getFaceColour(scene, bvh, ray, hitInfo);
}

// NOTE(Mathijs): separate function to make recursion easier (could also be done with lambda + std::function).
static glm::vec3 getFinalColor(const Scene& scene, const BoundingVolumeHierarchy& bvh, Ray ray)
{
    HitInfo hitInfo;
    //return recursiveRayTracing(scene, bvh, ray, hitInfo);
    int level = 10;

    if (bvh.intersect(ray, hitInfo)) {
        glm::vec3 res{};
        //Glossy reflection
        if (hitInfo.material.ks != glm::vec3{ 0 } && hitInfo.material.ks.x < 0.9) {
            return recursiveRayTracing(scene, bvh, ray, hitInfo, 1, 5, hitInfo.material.kd, true);
        }
        else if (hitInfo.material.ks != glm::vec3{ 0 } && hitInfo.material.ks.x > 0.9) {
            return recursiveRayTracing(scene, bvh, ray, hitInfo, 1, 5, hitInfo.material.kd, false);
        }
        
        else {
            res = getFaceColour(scene, bvh, ray, hitInfo);
        }
        drawRay(ray, res);
        return res;
    }
    else {
        // Draw a red debug ray if the ray missed.
        drawRay(ray, glm::vec3(1.0f, 0.0f, 0.0f));
        // Set the color of the pixel to black if the ray misses.
        return glm::vec3(0.0f);
    }


}

float boxFilter(glm::vec3 bpix[960][1200], int i, int j, int col, int filterSize) {
    filterSize = std::max(1, filterSize);
    float sum = 0.0f;
    for (int x = -filterSize; x < filterSize + 1; ++x) {
        if ((x + j > 0) && (x + j < windowResolution.x)) {
            for (int y = -filterSize; y < filterSize + 1; ++y) {
                if ((y + i >= 0) && (y + i < windowResolution.y)) {
                    sum = sum + bpixels[i + y][j + x][col];
                }
            }
        }
    }
    int multvalue = (2 * filterSize + 1) * (2 * filterSize + 1);
    sum = sum / multvalue;
    return sum;
}

glm::vec3 getBrightPass(glm::vec3 color) {
    glm::vec3 result = glm::vec3{ 0.0f };
    for (int i = 0; i < 3; i++) {
        if (color[i] > 1.0f) {
            result[i] = color[i];
        }
    }
    return result;
}

static void setOpenGLMatrices(const Trackball& camera);
static void renderOpenGL(const Scene& scene, const Trackball& camera, int selectedLight);

// This is the main rendering function. You are free to change this function in any way (including the function signature).
static void renderRayTracing(const Scene& scene, const Trackball& camera, const BoundingVolumeHierarchy& bvh, Screen& screen)
{

#ifdef USE_OPENMP
#pragma omp parallel for
#endif
    


    for (int y = 0; y < windowResolution.y; y++) {
        for (int x = 0; x != windowResolution.x; x++) {
            // NOTE: (-1, -1) at the bottom left of the screen, (+1, +1) at the top right of the screen.
            const glm::vec2 normalizedPixelPos{
                float(x) / windowResolution.x * 2.0f - 1.0f,
                float(y) / windowResolution.y * 2.0f - 1.0f
            };
            const Ray cameraRay = camera.generateRay(normalizedPixelPos);
            if (enableBloom && renderingToFile) {
                glm::vec3 col = getFinalColor(scene, bvh, cameraRay);
                pixels[y][x] = col;
                glm::vec3 bcol = getBrightPass(col);
                bpixels[y][x] = bcol;
            }
            if (enableAi && renderingToFile) {
                glm::vec3 col = getFinalColor(scene, bvh, cameraRay);
                pixels[y][x] = col;
            }
            else {
                screen.setPixel(x, y, getFinalColor(scene, bvh, cameraRay));
            }
        }
    }

    if (enableBloom && renderingToFile) {
        Screen screenBB = Screen{ windowResolution };

        for (int y = 0; y < windowResolution.y; y++) {
            for (int x = 0; x != windowResolution.x; x++) {
                glm::vec3 color = glm::vec3{ 0.0f };
                for (int i =0; i < 3; i++) {
                    color[i] = boxFilter(bpixels, y, x, i, bloomFilterSize);
                }
                screenBB.setPixel(x, y, color+pixels[y][x]);
            }
        }

        screen.clear(glm::vec3{ 1.0f });
        screen = screenBB;
    }

    if (enableAi && renderingToFile) {
        Screen screenAA = Screen{ windowResolution / 2 };
        int xO = 0;
        int yO = 0;
        yO = 0;

        for (int y = 0; y < windowResolution.y / 2; y++) {
            xO = 0;
            for (int x = 0; x != windowResolution.x / 2; x++) {
                glm::vec3 color = glm::vec3{ 0.0f };
                if ((yO + 1) < windowResolution.x && (xO + 1) < windowResolution.y) {
                    color = pixels[yO][xO] + pixels[yO][xO + 1] + pixels[yO + 1][xO] + pixels[yO + 1][xO + 1];
                    color = glm::vec3(color.x / 4, color.y / 4, color.z / 4);
                    screenAA.setPixel(x, y, color);
                    xO += 2;
                }
            }
            yO += 2;

        }


        screen.clear(glm::vec3{ 1.0f });
        screen = screenAA;
    }

}

int main(int argc, char** argv)
{
    initialize(32);  // initialize random vectors with 8 samples
    Trackball::printHelp();
    std::cout << "\n Press the [R] key on your keyboard to create a ray towards the mouse cursor" << std::endl
        << std::endl;
    Window window{ "Final Project - Part 2", windowResolution, OpenGLVersion::GL2 };
    Screen screen{ windowResolution };
    Trackball camera{ &window, glm::radians(50.0f), 3.0f };
    camera.setCamera(glm::vec3(0.0f, 0.0f, 0.0f), glm::radians(glm::vec3(20.0f, 20.0f, 0.0f)), 3.0f);

    SceneType sceneType{ SceneType::SingleTriangle };
    std::optional<Ray> optDebugRay;
    Scene scene = loadScene(sceneType, dataPath);
    BoundingVolumeHierarchy bvh{ &scene };

    int bvhDebugLevel = 0;
    bool debugBVH{ false };



    ViewMode viewMode{ ViewMode::Rasterization };

    window.registerKeyCallback([&](int key, int /* scancode */, int action, int /* mods */) {
        if (action == GLFW_PRESS) {
            switch (key) {
            case GLFW_KEY_R: {
                // Shoot a ray. Produce a ray from camera to the far plane.
                const auto tmp = window.getNormalizedCursorPos();
                optDebugRay = camera.generateRay(tmp * 2.0f - 1.0f);
                viewMode = ViewMode::Rasterization;
            } break;
            case GLFW_KEY_ESCAPE: {
                window.close();
            } break;
            };
        }
        });

    int selectedLight{ 0 };
    while (!window.shouldClose()) {
        window.updateInput();

        // === Setup the UI ===
        ImGui::Begin("Final Project - Part 2");
        {
            constexpr std::array items{ "SingleTriangle", "Cube", "Cornell Box (with mirror)", "Cornell Box (spherical light and mirror)","Cornell Box Gloss", "Monkey", "Dragon", /* "AABBs",*/ "Spheres", /*"Mixed",*/ "Custom" };
            if (ImGui::Combo("Scenes", reinterpret_cast<int*>(&sceneType), items.data(), int(items.size()))) {
                optDebugRay.reset();
                scene = loadScene(sceneType, dataPath);
                bvh = BoundingVolumeHierarchy(&scene);
                if (optDebugRay) {
                    HitInfo dummy{};
                    bvh.intersect(*optDebugRay, dummy);
                }
            }
        }
        {
            constexpr std::array items{ "Rasterization", "Ray Traced" };
            ImGui::Combo("View mode", reinterpret_cast<int*>(&viewMode), items.data(), int(items.size()));
        }
        if (ImGui::Button("Render to file")) {
            {
                renderingToFile = true;
                using clock = std::chrono::high_resolution_clock;
                const auto start = clock::now();
                Screen prev = screen;
                renderRayTracing(scene, camera, bvh, screen);
                const auto end = clock::now();
                std::cout << "Time to render image: " << std::chrono::duration<float, std::milli>(end - start).count() << " milliseconds" << std::endl;
                screen.writeBitmapToFile(outputPath / "render.bmp");
                renderingToFile = false;
                screen = prev;

            }
            //screen.writeBitmapToFile(outputPath / "render.bmp");
        }
        ImGui::Spacing();
        ImGui::Separator();
        ImGui::Text("Debugging");
        if (viewMode == ViewMode::Rasterization) {
            ImGui::Checkbox("Draw BVH", &debugBVH);
            if (debugBVH)
                ImGui::SliderInt("BVH Level", &bvhDebugLevel, 0, bvh.numLevels() - 1);
        }

        ImGui::Spacing();
        ImGui::Separator();
        ImGui::Text("Lights");  
        if (!scene.pointLights.empty() || !scene.sphericalLight.empty()) {
            {
                std::vector<std::string> options;
                for (size_t i = 0; i < scene.pointLights.size(); i++) {
                    options.push_back("Point Light " + std::to_string(i + 1));
                }
                for (size_t i = 0; i < scene.sphericalLight.size(); i++) {
                    options.push_back("Spherical Light " + std::to_string(i + 1));
                }

                std::vector<const char*> optionsPointers;
                std::transform(std::begin(options), std::end(options), std::back_inserter(optionsPointers),
                    [](const auto& str) { return str.c_str(); });

                ImGui::Combo("Selected light", &selectedLight, optionsPointers.data(), static_cast<int>(optionsPointers.size()));
            }

            {
                const auto showLightOptions = [](auto& light) {
                    ImGui::DragFloat3("Light position", glm::value_ptr(light.position), 0.01f, -3.0f, 3.0f);
                    ImGui::ColorEdit3("Light color", glm::value_ptr(light.color));
                    if constexpr (std::is_same_v<std::decay_t<decltype(light)>, SphericalLight>) {
                        ImGui::DragFloat("Light radius", &light.radius, 0.01f, 0.01f, 0.5f);
                    }
                };
                if (selectedLight < static_cast<int>(scene.pointLights.size())) {
                    // Draw a big yellow sphere and then the small light sphere on top.
                    showLightOptions(scene.pointLights[selectedLight]);
                }
                else {
                    // Draw a big yellow sphere and then the smaller light sphere on top.
                    showLightOptions(scene.sphericalLight[selectedLight - scene.pointLights.size()]);
                }
            }
        }

        if (ImGui::Button("Add point light")) {
            scene.pointLights.push_back(PointLight{ glm::vec3(0.0f), glm::vec3(1.0f) });
            selectedLight = int(scene.pointLights.size() - 1);
        }
        if (ImGui::Button("Add spherical light")) {
            scene.sphericalLight.push_back(SphericalLight{ glm::vec3(0.0f), 0.1f, glm::vec3(1.0f) });
            selectedLight = int(scene.pointLights.size() + scene.sphericalLight.size() - 1);
        }
        if (ImGui::Button("Remove selected light")) {
            if (selectedLight < static_cast<int>(scene.pointLights.size())) {
                scene.pointLights.erase(std::begin(scene.pointLights) + selectedLight);
            }
            else {
                scene.sphericalLight.erase(std::begin(scene.sphericalLight) + (selectedLight - scene.pointLights.size()));
            }
            selectedLight = 0;
        }

        ImGui::Spacing();
        ImGui::Separator();
        ImGui::Text("Soft Shadows");
        if (viewMode == ViewMode::Rasterization) {
            ImGui::SliderInt("Samples", &sampleLevel, 0, 32);
        }
        ImGui::Spacing();
        ImGui::Separator();
        ImGui::Text("Glossy reflections");
        if (viewMode == ViewMode::Rasterization) {
            ImGui::SliderInt("Gloss rays", &glossLevel, 0, 8);
        }
        ImGui::Text("Antialiasing");
        if (viewMode == ViewMode::Rasterization) {
            ImGui::Checkbox("Enable antialiasing", &enableAi);
        }
        ImGui::Text("Bloom Filter");
        if (viewMode == ViewMode::Rasterization) {
            ImGui::Checkbox("Enable bloomfilter", &enableBloom);
        }
        if (viewMode == ViewMode::Rasterization) {
            ImGui::SliderInt("Filter Size", &bloomFilterSize, 0, 10);
        }
        // Clear screen.
        glClearDepth(1.0f);
        glClearColor(0.0, 0.0, 0.0, 0.0);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

        // Draw either using OpenGL (rasterization) or the ray tracing function.
        switch (viewMode) {
        case ViewMode::Rasterization: {
            glPushAttrib(GL_ALL_ATTRIB_BITS);
            renderOpenGL(scene, camera, selectedLight);
            if (optDebugRay) {
                // Call getFinalColor for the debug ray. Ignore the result but tell the function that it should
                // draw the rays instead.
                enableDrawRay = true;
                (void)getFinalColor(scene, bvh, *optDebugRay);
                enableDrawRay = false;
            }
            glPopAttrib();
        } break;
        case ViewMode::RayTracing: {
            screen.clear(glm::vec3(0.0f));
            renderRayTracing(scene, camera, bvh, screen);
            screen.setPixel(0, 0, glm::vec3(1.0f));
            screen.draw(); // Takes the image generated using ray tracing and outputs it to the screen using OpenGL.
        } break;
        default:
            break;
        };

        if (debugBVH) {
            glPushAttrib(GL_ALL_ATTRIB_BITS);
            setOpenGLMatrices(camera);
            glDisable(GL_LIGHTING);
            glEnable(GL_DEPTH_TEST);

            // Enable alpha blending. More info at:
            // https://learnopengl.com/Advanced-OpenGL/Blending
            glEnable(GL_BLEND);
            glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
            bvh.debugDraw(bvhDebugLevel);
            glPopAttrib();
        }

        ImGui::End();
        window.swapBuffers();
    }

    return 0; // execution never reaches this point
}

static void setOpenGLMatrices(const Trackball& camera)
{
    // Load view matrix.
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    const glm::mat4 viewMatrix = camera.viewMatrix();
    glMultMatrixf(glm::value_ptr(viewMatrix));

    // Load projection matrix.
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    const glm::mat4 projectionMatrix = camera.projectionMatrix();
    glMultMatrixf(glm::value_ptr(projectionMatrix));
}

static void renderOpenGL(const Scene& scene, const Trackball& camera, int selectedLight)
{
    // Normals will be normalized in the graphics pipeline.
    glEnable(GL_NORMALIZE);
    // Activate rendering modes.
    glEnable(GL_DEPTH_TEST);
    // Draw front and back facing triangles filled.
    glPolygonMode(GL_FRONT, GL_FILL);
    glPolygonMode(GL_BACK, GL_FILL);
    // Interpolate vertex colors over the triangles.
    glShadeModel(GL_SMOOTH);
    setOpenGLMatrices(camera);

    glDisable(GL_LIGHTING);
    // Render point lights as very small dots
    for (const auto& light : scene.pointLights)
        drawSphere(light.position, 0.01f, light.color);
    for (const auto& light : scene.sphericalLight)
        drawSphere(light.position, light.radius, light.color);

    if (!scene.pointLights.empty() || !scene.sphericalLight.empty()) {
        if (selectedLight < static_cast<int>(scene.pointLights.size())) {
            // Draw     a big yellow sphere and then the small light sphere on top.
            const auto& light = scene.pointLights[selectedLight];
            drawSphere(light.position, 0.05f, glm::vec3(1, 1, 0));
            glDisable(GL_DEPTH_TEST);
            drawSphere(light.position, 0.01f, light.color);
            glEnable(GL_DEPTH_TEST);
        }
        else {
            // Draw a big yellow sphere and then the smaller light sphere on top.
            const auto& light = scene.sphericalLight[selectedLight - scene.pointLights.size()];
            drawSphere(light.position, light.radius + 0.01f, glm::vec3(1, 1, 0));
            glDisable(GL_DEPTH_TEST);
            drawSphere(light.position, light.radius, light.color);
            glEnable(GL_DEPTH_TEST);
        }
    }

    // Activate the light in the legacy OpenGL mode.
    glEnable(GL_LIGHTING);

    int i = 0;
    const auto enableLight = [&](const auto& light) {
        glEnable(GL_LIGHT0 + i);
        const glm::vec4 position4{ light.position, 1 };
        glLightfv(GL_LIGHT0 + i, GL_POSITION, glm::value_ptr(position4));
        const glm::vec4 color4{ glm::clamp(light.color, 0.0f, 1.0f), 1.0f };
        const glm::vec4 zero4{ 0.0f, 0.0f, 0.0f, 1.0f };
        glLightfv(GL_LIGHT0 + i, GL_AMBIENT, glm::value_ptr(zero4));
        glLightfv(GL_LIGHT0 + i, GL_DIFFUSE, glm::value_ptr(color4));
        glLightfv(GL_LIGHT0 + i, GL_SPECULAR, glm::value_ptr(zero4));
        // NOTE: quadratic attenuation doesn't work like you think it would in legacy OpenGL.
        // The distance is not in world space but in NDC space!
        glLightf(GL_LIGHT0 + i, GL_CONSTANT_ATTENUATION, 1.0f);
        glLightf(GL_LIGHT0 + i, GL_LINEAR_ATTENUATION, 0.0f);
        glLightf(GL_LIGHT0 + i, GL_QUADRATIC_ATTENUATION, 0.0f);
        i++;
    };
    for (const auto& light : scene.pointLights)
        enableLight(light);
    for (const auto& light : scene.sphericalLight)
        enableLight(light);

    // Draw the scene and the ray (if any).
    drawScene(scene);

    // Draw a colored sphere at the location at which the trackball is looking/rotating around.
    glDisable(GL_LIGHTING);
    drawSphere(camera.lookAt(), 0.01f, glm::vec3(0.2f, 0.2f, 1.0f));
}
