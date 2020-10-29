#include "bounding_volume_hierarchy.h"
#include "draw.h"
#include "CustomCompX.cpp"
#include "CustomCompY.cpp"
#include "CustomCompZ.cpp"
#include <cmath>
#include <glm\common.hpp>
#include <queue>
#include <iostream>
#include <queue>

std::vector<BoundingVolumeHierarchy::Node> nodes;
std::vector<BoundingVolumeHierarchy::Node> rootNodes;

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

int shrinkBox(AxisAlignedBox box, std::vector<Triangle> scene, Mesh mesh, int axis, int level, int meshNum) {    
    
    //Stopping the algorithm if the scene becomes empty, or the maximum level is reached.
    if (scene.size() == 0) {
        std::vector<int> children {};
        BoundingVolumeHierarchy::Node node{ AxisAlignedBox{ glm::vec3{0}, glm::vec3{0} }, level, children, true, meshNum, scene };
        nodes.push_back(node);
        if (node.level == 0) {
            rootNodes.push_back(node);
        }
        return nodes.size() - 1;
    }
    else if (scene.size() == 1 || level > 6) {
        AxisAlignedBox box = setBox(scene, mesh);
        std::vector<int> children{};
        BoundingVolumeHierarchy::Node node{ box, level, children, true, meshNum, scene };
        nodes.push_back(node);
        if (node.level == 0) {
            rootNodes.push_back(node);
        }
        return nodes.size()-1;
    }
    
    // splitting along the X axis
    if (axis == 0) {
        CustomCompX compX = CustomCompX(mesh);
        std::sort(scene.begin(), scene.end(), compX);

        std::vector<Triangle> firstSet {};
        std::vector<Triangle> secondSet{};

        for (int i = 0; i < scene.size() / 2; i++) {
            firstSet.push_back(scene[i]);
        }

        for (int i = scene.size() / 2; i < scene.size(); i++) {
            secondSet.push_back(scene[i]);
        }

        AxisAlignedBox box1 = setBox(firstSet, mesh);
        AxisAlignedBox box2 = setBox(secondSet, mesh);

        int node1 = shrinkBox(box1, firstSet, mesh, 1, level + 1, meshNum);
        int node2 = shrinkBox(box2, secondSet, mesh, 1, level + 1, meshNum);
        std::vector<int> children{ node1 , node2};
        BoundingVolumeHierarchy::Node currentNode{ box, level, children, false };
        if (level == 0) {
            rootNodes.push_back(currentNode);
        }
        nodes.push_back(currentNode);
        return nodes.size() - 1;
            
            
    }

    // Splitting along the Y axis
    else if (axis == 1) {
        CustomCompY compY = CustomCompY(mesh);
        std::sort(scene.begin(), scene.end(), compY);

        std::vector<Triangle> firstSet{};
        std::vector<Triangle> secondSet{};

        for (int i = 0; i < scene.size() / 2; i++) {
            firstSet.push_back(scene[i]);
        }

        for (int i = scene.size() / 2 ; i < scene.size(); i++) {
            secondSet.push_back(scene[i]);
        }

        AxisAlignedBox box1 = setBox(firstSet, mesh);
        AxisAlignedBox box2 = setBox(secondSet, mesh); 
        int node1 = shrinkBox(box1, firstSet, mesh, 2, level + 1, meshNum);
        int node2 = shrinkBox(box2, secondSet, mesh, 2, level + 1, meshNum);
        std::vector<int> children{ node1, node2 };

        BoundingVolumeHierarchy::Node currentNode{ box, level, children, false };
        nodes.push_back(currentNode);
        return nodes.size()-1;

    }

    //Splitting along the Z axis - the only difference between the 3 axes are the custom comparator and what the next axis which will be split will be
    else if (axis == 2) {
        CustomCompZ compZ = CustomCompZ(mesh);
        std::sort(scene.begin(), scene.end(), compZ);

        std::vector<Triangle> firstSet{};
        std::vector<Triangle> secondSet{};

        for (int i = 0; i < scene.size() / 2; i++) {
            firstSet.push_back(scene[i]);
        }

        for (int i = scene.size() / 2; i < scene.size(); i++) {
            secondSet.push_back(scene[i]);
        }

        AxisAlignedBox box1 = setBox(firstSet, mesh);
        AxisAlignedBox box2 = setBox(secondSet, mesh);
        int node1 = shrinkBox(box1, firstSet, mesh, 0, level + 1, meshNum);
        int node2 = shrinkBox(box2, secondSet, mesh, 0, level + 1, meshNum);
        std::vector<int> children{ node1, node2 };
        BoundingVolumeHierarchy::Node currentNode{ box, level, children, false };
        nodes.push_back(currentNode);
        return nodes.size() - 1;
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
    rootNodes = {};
    int i = 0;
    for (Mesh mesh : pScene->meshes) {
        AxisAlignedBox overallBox = largestBoxSize(pScene);
        std::vector<Triangle> totalTriangles = getTriangles(mesh);
        shrinkBox(overallBox, totalTriangles, mesh, 0, 0, i);
        i++;
    }    

    // as an example of how to iterate over all meshes in the scene, look at the intersect method below
}



// Use this function to visualize your BVH. This can be useful for debugging. Use the functions in
// draw.h to draw the various shapes. We have extended the AABB draw functions to support wireframe
// mode, arbitrary colors and transparency.
void BoundingVolumeHierarchy::debugDraw(int level)
{

    int newLevel = BoundingVolumeHierarchy::numLevels();
    for (int i = 0; i < nodes.size(); i++) {
        if (nodes[i].level == level) {
            drawAABB(nodes[i].box, DrawMode::Wireframe);
        }
    }
}

int BoundingVolumeHierarchy::numLevels() const
{
    return 7;
}

bool intersectRayWithShapeNoChange(const AxisAlignedBox& box, Ray& ray)
{
    float tinX = glm::min(((box.lower.x - ray.origin.x) / ray.direction.x), ((box.upper.x - ray.origin.x) / ray.direction.x));
    float tOutX = glm::max(((box.lower.x - ray.origin.x) / ray.direction.x), ((box.upper.x - ray.origin.x) / ray.direction.x));
    float tinY = glm::min(((box.lower.y - ray.origin.y) / ray.direction.y), ((box.upper.y - ray.origin.y) / ray.direction.y));
    float tOutY = glm::max(((box.lower.y - ray.origin.y) / ray.direction.y), ((box.upper.y - ray.origin.y) / ray.direction.y));
    float tinZ = glm::min(((box.lower.z - ray.origin.z) / ray.direction.z), ((box.upper.z - ray.origin.z) / ray.direction.z));
    float tOutZ = glm::max(((box.lower.z - ray.origin.z) / ray.direction.z), ((box.upper.z - ray.origin.z) / ray.direction.z));

    float tin = glm::max(tinX, glm::max(tinY, tinZ));
    float tout = glm::min(tOutX, glm::min(tOutY, tOutZ));

    if (tin > tout || tout < 0) {
        return false;
    }

    /*if (tin < ray.t) {
        ray.t = tin;
    }*/

    return true;
}

// Return true if something is hit, returns false otherwise. Only find hits if they are closer than t stored
// in the ray and if the intersection is on the correct side of the origin (the new t >= 0). Replace the code
// by a bounding volume hierarchy acceleration structure as described in the assignment. You can change any
bool BoundingVolumeHierarchy::intersect(Ray& ray, HitInfo& hitInfo) const
{
    bool hit = false;
    //Ray minRay = ray;
    Ray minRay = ray;

    for (const auto& node : rootNodes) {

        if (intersectRayWithShapeNoChange(node.box, minRay)) {

            std::queue<Node> queue;
            queue.push(node);

            std::vector<Node> leafNodes{};

            while (!queue.empty()) {
                Node current = queue.front();
                queue.pop();

                if (current.leafNode) {
                    leafNodes.push_back(current);
                    continue;
                }

                if (intersectRayWithShapeNoChange(nodes[current.children[0]].box, minRay)) {
                    queue.push(nodes[current.children[0]]);
                }

                if (intersectRayWithShapeNoChange(nodes[current.children[1]].box, minRay)) {
                    queue.push(nodes[current.children[1]]);
                }
            }
            for (int i = 0; i < leafNodes.size(); i++) {
                for (const auto& tri : leafNodes[i].triangles) {
                    const auto v0 = m_pScene->meshes[leafNodes[i].mesh].vertices[tri[0]];
                    const auto v1 = m_pScene->meshes[leafNodes[i].mesh].vertices[tri[1]];
                    const auto v2 = m_pScene->meshes[leafNodes[i].mesh].vertices[tri[2]];

                    if (intersectRayWithTriangle(v0.p, v1.p, v2.p, minRay, hitInfo)) {

                        if (minRay.t < ray.t) {
                            ray.t = minRay.t;
                            //drawRay(ray, glm::vec3{ 1.0f });
                            hitInfo.material = m_pScene->meshes[leafNodes[i].mesh].material;
                            hit = true;

                        }
                    }
                }
            }
        }
    }

    // Intersect with spheres.
    for (const auto& sphere : m_pScene->spheres)
        hit |= intersectRayWithShape(sphere, ray, hitInfo);
    return hit;
}