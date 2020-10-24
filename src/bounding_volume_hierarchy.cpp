#include "bounding_volume_hierarchy.h"
#include "draw.h"
#include "queue"
#include <cmath>
#include <glm/geometric.hpp>
#include <glm/gtx/component_wise.hpp>
#include <glm/vector_relational.hpp>

class CustomCompX {
    std::vector<Vertex> vertices;

    CustomCompX(std::vector<Vertex> vertices) {
        vertices = vertices;
    }

    bool operator() (Mesh mesh, Triangle t1, Triangle t2) const {
        const auto v1 = glm::min(glm::min(mesh.vertices[t1[0]].p.x, mesh.vertices[t1[1]].p.x), mesh.vertices[t1[2]].p.x);
        const auto v2 = glm::min(glm::min(mesh.vertices[t2[0]].p.x, mesh.vertices[t2[1]].p.x), mesh.vertices[t2[2]].p.x);
        return v1 > v2;
    }
};

AxisAlignedBox setBox(Scene* pScene,glm::vec3 min, glm::vec3 max) {
    glm::vec3 minPos = min;
    glm::vec3 maxPos = max;
    for (const auto& mesh : pScene->meshes) {
        for (const auto& tri : mesh.vertices) {
            if (tri.p.x < minPos.x ) {
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
}

AxisAlignedBox largestBoxSize(Scene* pScene) {
    return setBox(pScene, glm::vec3{ std::numeric_limits<float>::max() }, glm::vec3{ -std::numeric_limits<float>::max() });
}


/* 
*  recursive step that splits the bounding box to be smaller
* box is the axisalignedbox it is in
* scene is the scene of meshes the object is within
* axis is whether or not it will be split among the x or y axis the next time. true = x axis, false = y axis
* direction is the direction that the axis was split in, true = left/up, false = right/down
* level is the level of the bvh that we have split to. max is 5, it will not go beyond that.
*/

void shrinkBox(AxisAlignedBox box, Scene scene, bool axis, bool direction, int level) {
    
    if (axis) {
        std::priority_queue<Triangle, std::vector<float>, CustomCompX> queue;
        int count = 0;
        for (const auto& mesh : scene.meshes) {
            for (const auto& tri : mesh.triangles) {
                if (axis) {
                    float x = mesh.vertices[tri[0]].p.x;
                    queue.push(tri);

                }
                else {
                    queue.push(tri);
                }
                count++;
            }
        }

        int median = count / 2;
        if (count % 2 != 0) {
            if (direction) {
                median++;
            }
        }

        std::vector<Triangle> triangles = {};
        for (int i = 0; i < median; i++) {
            Triangle tri = queue.pop;
            triangles.insert(tri);
        }

        Triangle tri{};
        //mesh.vertices[tri[0]].;

    }
    
}


BoundingVolumeHierarchy::BoundingVolumeHierarchy(Scene* pScene)
    : m_pScene(pScene)
{
    largestBoxSize(pScene);
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

    // Draw the AABB as a (white) wireframe box.
    AxisAlignedBox aabb { glm::vec3(-.85f, -0.65f, -0.33f), glm::vec3(0.85f, 0.57f, 0.725f) };
    //drawAABB(aabb, DrawMode::Wireframe);

    drawAABB(aabb, DrawMode::Filled, glm::vec3(0.05f, 1.0f, 0.05f), 0.1);


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
