#include "bounding_volume_hierarchy.h"
#include "draw.h"
#include "CustomCompX.cpp"
#include "CustomCompY.cpp"
#include "CustomCompZ.cpp"
#include <cmath>
#include <glm\common.hpp>
#include <queue>
#include <iostream>

AxisAlignedBox setBox(std::vector<Triangle> triangles, std::vector<Mesh> meshes) {
    glm::vec3 minPos = glm::vec3{ std::numeric_limits<float>::max() };
    glm::vec3 maxPos = glm::vec3{ -std::numeric_limits<float>::min() };
    for (const Mesh mesh : meshes) {
        for (const Triangle tri : triangles) {
            for (int i = 0; i < 3; i++) {
                if (mesh.vertices[tri[i]].p.x < minPos.x) {
                    minPos.x = mesh.vertices[tri[i]].p.x;
                }
                if (mesh.vertices[tri[i]].p.y < minPos.y) {
                    minPos.y = mesh.vertices[tri[i]].p.y;
                }
                if (mesh.vertices[tri[i]].p.z < minPos.z) {
                    minPos.z = mesh.vertices[tri[i]].p.z;
                }
                if (mesh.vertices[tri[i]].p.x > maxPos.x) {
                    maxPos.x = mesh.vertices[tri[i]].p.x;
                }
                if (mesh.vertices[tri[i]].p.y > maxPos.y) {
                    maxPos.y = mesh.vertices[tri[i]].p.y;
                }
                if (mesh.vertices[tri[i]].p.z > maxPos.z) {
                    maxPos.z = mesh.vertices[tri[i]].p.z;
                }
            }
            
        }
    }
    

    AxisAlignedBox box{ minPos, maxPos };
    return box;
}

AxisAlignedBox largestBoxSize(Scene* pScene) {
    glm::vec3 minPos = glm::vec3{ std::numeric_limits<float>::max() };
    glm::vec3 maxPos = glm::vec3{ -std::numeric_limits<float>::min() };
    for (const auto& mesh : pScene->meshes) {
        for (const auto& tri : mesh.vertices) {
            if (tri.p.x < minPos.x) {
                minPos.x = tri.p.x;
            }
            if (tri.p.y < minPos.y) {
                minPos.y = tri.p.y;
            }
            if (tri.p.z < minPos.z) {
                minPos.z = tri.p.z;
            }
            if (tri.p.x > maxPos.x) {
                maxPos.x = tri.p.x;
            }
            if (tri.p.y > maxPos.y) {
                maxPos.y = tri.p.y;
            }
            if (tri.p.z > maxPos.z) {
                maxPos.z = tri.p.z;
            }
        }
    }

    AxisAlignedBox box{ minPos, maxPos };
    return box;
    //return setBox(pScene, glm::vec3{ std::numeric_limits<float>::max() }, glm::vec3{ -std::numeric_limits<float>::max() });
}

std::vector<Triangle> getTriangles(Scene* scene) {
    std::vector<Triangle> res = {};
    /*for (const auto& mesh : m_pScene->meshes) {
        for (const auto& tri : mesh.triangles) {
        }
    }*/
    for (const auto& mesh : scene->meshes) {
        for (const auto& tri : mesh.triangles) {
            res.push_back(tri);
        }
    }
    return res;
}

/* 
*  recursive step that splits the bounding box to be smaller
* box is the axisalignedbox it is in
* scene is the scene of meshes the object is within
* axis is whether or not it will be split among the x or y axis the next time. true = x axis, false = y axis
* direction is the direction that the axis was split in, true = left/up, false = right/down
* level is the level of the bvh that we have split to. max is 5, it will not go beyond that.
*/

void shrinkBox(AxisAlignedBox box, std::vector<Triangle> scene, std::vector<Mesh> mesh, int axis, int level) {
    if (level > 1 ) {
        return;
    }

    /*if (axis) {
        AxisAlignedBox box1 = setBox(scene, box.lower, glm::vec3{ box.upper.x / 2.0f, box.upper.y, box.upper.z });
        AxisAlignedBox box2 = setBox(scene, glm::vec3{ box.lower.x/2.0f, box.lower.y, box.lower.z }, box.upper);
        shrinkBox(box1, scene, !axis, level + 1);
        shrinkBox(box2, scene, !axis, level + 1);
    }

    else {
        AxisAlignedBox box1 = setBox(scene, box.lower, glm::vec3{ box.upper.x , box.upper.y*0.5f, box.upper.z });
        AxisAlignedBox box2 = setBox(scene, glm::vec3{ box.lower.x, box.lower.y*0.5f, box.lower.z }, box.upper);
        shrinkBox(box1, scene, !axis, level + 1);
        shrinkBox(box2, scene, !axis, level + 1);
    }*/
    
    if (axis == 0) {
        CustomCompX compX = CustomCompX(mesh.front());
        std::priority_queue<Triangle, std::vector<Triangle>, CustomCompX> queue (compX);
        //ree queue(CustomCompX);
        //    CustomCompX.setMesh(mesh.front());
        int count = 0;
        for (const auto& tri : scene) {
            if (axis) {
                //float x = tri[0].p.x;
                queue.push(tri);

            }
            else {
                queue.push(tri);
            }
            count++;
        }

        int median = count / 2;

        std::vector<Triangle> firstSet = {};
        for (int i = 0; i < median; i++) {
            firstSet.push_back(queue.top());
            queue.pop();
        }

        std::vector<Triangle> secondSet = {};
        while (!queue.empty()) {
            secondSet.push_back(queue.top());
            queue.pop();
        }

        /*std::vector<glm::vec3> test = {};    //these are only here to debug
        for (const Triangle vert : firstSet) {
            test.push_back((mesh.front().vertices[vert[0]].p + mesh.front().vertices[vert[1]].p + mesh.front().vertices[vert[2]].p)/3.0f);
        }

        std::vector<glm::vec3> test2 = {};
        for (const Triangle vert : secondSet) {
            test2.push_back((mesh.front().vertices[vert[0]].p + mesh.front().vertices[vert[1]].p + mesh.front().vertices[vert[2]].p) / 3.0f);
        }*/

        AxisAlignedBox box1 = setBox(firstSet, mesh);// , box.lower, box.upper);
        AxisAlignedBox box2 = setBox(secondSet, mesh);// , box.lower, box.upper);
        shrinkBox(box1, firstSet, mesh, 0, level++);
        shrinkBox(box2, secondSet, mesh, 0, level++);


    }
    /*else if (axis == 1) {
        CustomCompY compY = CustomCompY(mesh.front());
        std::priority_queue<Triangle, std::vector<Triangle>, CustomCompY> queue(compY);

        int count = 0;
        for (const auto& tri : scene) {
            if (axis) {
                queue.push(tri);

            }
            else {
                queue.push(tri);
            }
            count++;
        }

        int median = count / 2;

        std::vector<Triangle> firstSet = {};
        for (int i = 0; i < median; i++) {
            firstSet.push_back(queue.top());
            queue.pop();
        }

        std::vector<Triangle> secondSet = {};
        while (!queue.empty()) {
            secondSet.push_back(queue.top());
            queue.pop();
        }

        AxisAlignedBox box1 = setBox(firstSet, mesh);
        AxisAlignedBox box2 = setBox(secondSet, mesh); 
        shrinkBox(box1, firstSet, mesh, 2, level++);
        shrinkBox(box2, secondSet, mesh, 2, level++);
    }

    else if (axis == 2) {
        CustomCompZ compZ = CustomCompZ(mesh.front());
        std::priority_queue<Triangle, std::vector<Triangle>, CustomCompZ> queue(compZ);

        int count = 0;
        for (const auto& tri : scene) {
            if (axis) {
                queue.push(tri);

            }
            else {
                queue.push(tri);
            }
            count++;
        }

        int median = count / 2;

        std::vector<Triangle> firstSet = {};
        for (int i = 0; i < median; i++) {
            firstSet.push_back(queue.top());
            queue.pop();
        }

        std::vector<Triangle> secondSet = {};
        while (!queue.empty()) {
            secondSet.push_back(queue.top());
            queue.pop();
        }

        AxisAlignedBox box1 = setBox(firstSet, mesh);
        AxisAlignedBox box2 = setBox(secondSet, mesh);
        shrinkBox(box1, firstSet, mesh, 0, level++);
        shrinkBox(box2, secondSet, mesh, 0, level++);
    }*/
    else {
        std::cout << "There is a problem lmao" << std::endl;
    }
}


BoundingVolumeHierarchy::BoundingVolumeHierarchy(Scene* pScene)
    : m_pScene(pScene)
{
    AxisAlignedBox overallBox = largestBoxSize(pScene);
    std::vector<Triangle> totalTriangles = getTriangles(pScene);
    shrinkBox(overallBox, totalTriangles, pScene->meshes, 0, 0);

    // as an example of how to iterate over all meshes in the scene, look at the intersect method below
}



// Use this function to visualize your BVH. This can be useful for debugging. Use the functions in
// draw.h to draw the various shapes. We have extended the AABB draw functions to support wireframe
// mode, arbitrary colors and transparency.
void BoundingVolumeHierarchy::debugDraw(int level)
{

    // Draw the AABB as a transparent green box.
    //AxisAlignedBox aabb{ glm::vec3(-0.05f), glm::vec3(0.05f, 1.05f, 1.05f) };
    //drawShape(aabb, DrawMode::Filled, glm::vec3(0.0f, 1.0f, 0.0f), 0.2f);
    BoundingVolumeHierarchy::numLevels();
    // Draw the AABB as a (white) wireframe box.
    AxisAlignedBox aabb { glm::vec3(-.85f, -0.65f, -0.33f), glm::vec3(0.85f, 0.57f, 0.725f) };
    //drawAABB(aabb, DrawMode::Wireframe);

    /*AxisAlignedBox layer11{ glm::vec3(-.85f, -0.65f, -0.33f), glm::vec3(0.45f, 0.57f, 0.725f) };
    AxisAlignedBox layer12{ glm::vec3(-0.515f, -0.65f, -0.33f), glm::vec3(0.85f, 0.57f, 0.725f) };*/

    AxisAlignedBox layer11{ glm::vec3(-.85f, -0.65f, -0.23f), glm::vec3(0.85f, 0.57f, 0.725f) };
    AxisAlignedBox layer12{ glm::vec3(-0.534f, -0.65f, -0.33f), glm::vec3(0.534f, 0.433f, 0.0f) };
    AxisAlignedBox check{ glm::vec3(-.75f, -.55f, -0.170370504), glm::vec3(0.75f, 0.47f, 0.0f) };


    //drawAABB(aabb, DrawMode::Filled, glm::vec3(0.05f, 1.0f, 0.05f), 0.1);
    drawAABB(layer11, DrawMode::Wireframe);
    drawAABB(layer12, DrawMode::Wireframe);
    drawAABB(check, DrawMode::Wireframe);

}

int BoundingVolumeHierarchy::numLevels() const
{
    return 5;
}

// Return true if something is hit, returns false otherwise. Only find hits if they are closer than t stored
// in the ray and if the intersection is on the correct side of the origin (the new t >= 0). Replace the code
// by a bounding volume hierarchy acceleration structure as described in the assignment. You can change any
// file you like, including bounding_volume_hierarchy.h .
bool BoundingVolumeHierarchy::intersect(Ray& ray, HitInfo& hitInfo) const
{
    bool hit = false;
    Ray minRay = ray;
    // Intersect with all triangles of all meshes.
    for (const auto& mesh : m_pScene->meshes) {
        for (const auto& tri : mesh.triangles) {
            const auto v0 = mesh.vertices[tri[0]];
            const auto v1 = mesh.vertices[tri[1]];
            const auto v2 = mesh.vertices[tri[2]];
            if (intersectRayWithTriangle(v0.p, v1.p, v2.p, minRay, hitInfo)) {
                if (minRay.t < ray.t) {
                    ray.t = minRay.t;
                    hitInfo.material = mesh.material;
                    hit = true;
                }
            }
        }
    }
    // Intersect with spheres.
    for (const auto& sphere : m_pScene->spheres)
        hit |= intersectRayWithShape(sphere, ray, hitInfo);
    return hit;
}
