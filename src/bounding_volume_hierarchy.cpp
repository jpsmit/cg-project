#include "bounding_volume_hierarchy.h"
#include "draw.h"
#include "CustomCompX.cpp"
#include "CustomCompY.cpp"
#include "CustomCompZ.cpp"
#include <cmath>
#include <glm\common.hpp>
#include <queue>
#include <iostream>

std::vector<BoundingVolumeHierarchy::Node> nodes;

/*Sets the size of the AxisAlignedBox that will be created, which surrounds the given triangles for a specific mesh.
* triangles is a vector with all the triangles that will be considered for this bounding box
* mesh is the mesh with vertex positions that correspond to the set of triangles.
*/
AxisAlignedBox setBox(std::vector<Triangle> triangles, Mesh mesh) {
    glm::vec3 minPos = glm::vec3{ std::numeric_limits<float>::max() };
    glm::vec3 maxPos = glm::vec3{ -std::numeric_limits<float>::min() };
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
    

    AxisAlignedBox box{ minPos, maxPos };
    return box;
}

/*Creates an AxisAlignedBox surrounding the edges of the entire scene.
* the scene that we are looking at
*/
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

/* Gets all the triangles within a specific mesh
*  mesh is the mesh that we will get all triangles from
*/
std::vector<Triangle> getTriangles(Mesh mesh) {
    std::vector<Triangle> res = {};
    for (const auto& tri : mesh.triangles) {
        res.push_back(tri);
    }
    return res;
}

/* 
*  recursive step that splits the bounding box to be smaller, returns the current node, and inserts itself into the vector of nodes. 
*  The vector recieves nodes in depth first order so the root node will always be the last in the vector. Use this for traversal.
* box is the axisalignedbox it is in
* scene is the set of triangles that will be checked in the method. This set halves each time, being split among the x, y or z axes.
* mesh is the mesh that is being used. It must correspond to the set of triangles in scene.
* axis determines which axis it will be split among. 0 = x, 1 = y, 2 = z. It should always start with being split on x.
* level is the number of times the current box has been split. Max number of times is 11, and it will stop there or earlier if each triangle is found.
*/

BoundingVolumeHierarchy::Node shrinkBox(AxisAlignedBox box, std::vector<Triangle> scene, Mesh mesh, int axis, int level) {    
    
    //Stopping the algorithm if the scene becomes empty, or the maximum level is reached.
    if (scene.size() == 0) {
        std::vector<BoundingVolumeHierarchy::Node> children {};
        BoundingVolumeHierarchy::Node node{ AxisAlignedBox{ glm::vec3{0}, glm::vec3{0} }, level, children, true };
        nodes.push_back(node);
        return node;
    }
    else if (scene.size() == 1 || level > 10) {
        AxisAlignedBox box = setBox(scene, mesh);
        std::vector<BoundingVolumeHierarchy::Node> children{};
        BoundingVolumeHierarchy::Node node{ box, level, children, true };
        nodes.push_back(node);
        return node;
    }
    
    // splitting among the X axis
    if (axis == 0) {
        CustomCompX compX = CustomCompX(mesh);
        std::priority_queue<Triangle, std::vector<Triangle>, CustomCompX> queue (compX);;
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

        BoundingVolumeHierarchy::Node node1 = shrinkBox(box1, firstSet, mesh, 1, level + 1);
        BoundingVolumeHierarchy::Node node2 = shrinkBox(box2, secondSet, mesh, 1, level + 1);
        std::vector<BoundingVolumeHierarchy::Node> children{node1, node2};
        BoundingVolumeHierarchy::Node currentNode{ box, level, children, false };
        nodes.push_back(currentNode);
        return currentNode;
            
            
    }

    // Splitting along the Y axis
    else if (axis == 1) {
        CustomCompY compY = CustomCompY(mesh);
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
        BoundingVolumeHierarchy::Node node1 = shrinkBox(box1, firstSet, mesh, 2, level + 1);
        BoundingVolumeHierarchy::Node node2 = shrinkBox(box2, secondSet, mesh, 2, level + 1);
        std::vector<BoundingVolumeHierarchy::Node> children{ node1, node2 };
        BoundingVolumeHierarchy::Node currentNode{ box, level, children, false };
        nodes.push_back(currentNode);
        return currentNode;

    }

    //Splitting along the Z axis - the only difference between the 3 axes are the custom comparator and what the next axis which will be split will be
    else if (axis == 2) {
        CustomCompZ compZ = CustomCompZ(mesh);
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
        BoundingVolumeHierarchy::Node node1 = shrinkBox(box1, firstSet, mesh, 0, level + 1);
        BoundingVolumeHierarchy::Node node2 = shrinkBox(box2, secondSet, mesh, 0, level + 1);
        std::vector<BoundingVolumeHierarchy::Node> children{ node1, node2 };
        BoundingVolumeHierarchy::Node currentNode{ box, level, children, false };
        nodes.push_back(currentNode);
        return currentNode;
    }

    //if it somehow gets here we have a problem
    else {
        std::cout << "There is a problem lmao" << std::endl;
    }
}

/* constructor of the bvh. It initializes the vector 'nodes', and creates a bvh box for each mesh in the scene.
*/
BoundingVolumeHierarchy::BoundingVolumeHierarchy(Scene* pScene)
    : m_pScene(pScene)
{
    nodes = {};
    for (Mesh mesh : pScene->meshes) {
        AxisAlignedBox overallBox = largestBoxSize(pScene);
        std::vector<Triangle> totalTriangles = getTriangles(mesh);
        shrinkBox(overallBox, totalTriangles, mesh, 0, 0);
    }
    

    // as an example of how to iterate over all meshes in the scene, look at the intersect method below
}



// Use this function to visualize your BVH. This can be useful for debugging. Use the functions in
// draw.h to draw the various shapes. We have extended the AABB draw functions to support wireframe
// mode, arbitrary colors and transparency.
void BoundingVolumeHierarchy::debugDraw(int level)
{

    int newLevel = BoundingVolumeHierarchy::numLevels();
    for (Node node : nodes) {
        if (node.level == level) {
            drawAABB(node.box, DrawMode::Wireframe);
        }
    }
}

int BoundingVolumeHierarchy::numLevels() const
{
    return 11;
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
