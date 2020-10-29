#pragma once
#include "ray_tracing.h"
#include "scene.h"

class BoundingVolumeHierarchy {
public:
    BoundingVolumeHierarchy(Scene* pScene);

    // Use this function to visualize your BVH. This can be useful for debugging.

    struct Node {
        AxisAlignedBox box{ glm::vec3{0.0f}, glm::vec3{1.0f} };
        int level = 0;
        std::vector<int> children {};
        bool leafNode = true;
        int mesh = 0;
        std::vector<Triangle> triangles;
    };

    void debugDraw(int level);
    int numLevels() const;

    // Return true if something is hit, returns false otherwise.
    // Only find hits if they are closer than t stored in the ray and the intersection
    // is on the correct side of the origin (the new t >= 0).
    bool intersect(Ray& ray, HitInfo& hitInfo) const;

private:
    Scene* m_pScene;
};
