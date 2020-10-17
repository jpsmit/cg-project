#include "ray_tracing.h"
#include "disable_all_warnings.h"
// Suppress warnings in third-party code.
DISABLE_WARNINGS_PUSH()
#include <glm/geometric.hpp>
#include <glm/gtx/component_wise.hpp>
#include <glm/vector_relational.hpp>
DISABLE_WARNINGS_POP()
#include <cmath>
#include <iostream>
#include <limits>


Ray calculateReflectionRay(Ray ray, HitInfo hitInfo) {
    glm::vec3 vertexPos = ray.origin + (ray.direction * ray.t);
    glm::vec3 incident = ray.origin - vertexPos;
    incident = glm::normalize(incident);
    glm::vec3 reflectionDirection = glm::reflect(incident, hitInfo.normal);
    reflectionDirection = glm::normalize(-reflectionDirection);
    Ray reflectionRay = Ray{ vertexPos, reflectionDirection};
    return reflectionRay;

}

bool pointInTriangle(const glm::vec3& v0, const glm::vec3& v1, const glm::vec3& v2, const glm::vec3& n, const glm::vec3& p)
{
    float s0 = glm::dot(glm::cross((p - v0), (v2 - v0)), n);
    float s1 = glm::dot(glm::cross((p - v2), (v1 - v2)), n);
    float s2 = glm::dot(glm::cross((p - v1), (v0 - v1)), n);

    if (s0 >= 0 && s1 >= 0 && s2 >= 0) {
        return true;
    }
    else if (s0 <= 0 && s1 <= 0 && s2 <= 0) {
        return true;
    }

    return false;
}

bool intersectRayWithPlane(const Plane& plane, Ray& ray)
{
    float dotprod = glm::dot(plane.normal, ray.direction);
    if (dotprod != 0) {
        float d = (plane.D - glm::dot(plane.normal, ray.origin)) / dotprod;
        if (d < ray.t) {
            ray.t = d;
        }

        return true;
    }

    return false;
}

Plane trianglePlane(const glm::vec3& v0, const glm::vec3& v1, const glm::vec3& v2)
{
    Plane plane;
    plane.normal = glm::normalize(glm::cross((v2 - v0), (v1 - v0)));
    plane.D = glm::dot(v0, glm::normalize(glm::cross((v2 - v0), (v1 - v0))));
    return plane;
}

/// Input: the three vertices of the triangle
/// Output: if intersects then modify the hit parameter ray.t and return true, otherwise return false
bool intersectRayWithTriangle(const glm::vec3& v0, const glm::vec3& v1, const glm::vec3& v2, Ray& ray, HitInfo& hitInfo)
{
    Plane plane = trianglePlane(v0, v1, v2);

    Ray newRay = ray;

    if (intersectRayWithPlane(plane, newRay)) {
        glm::vec3 p = newRay.origin + newRay.direction * newRay.t;

        if (pointInTriangle(v0, v1, v2, plane.normal, p)) {
            if (newRay.t < ray.t) {
                hitInfo.normal = plane.normal;
                ray.t = newRay.t;
            }
            return true;
        }
    }
    return false;
}

/// Input: a sphere with the following attributes: sphere.radius, sphere.center
/// Output: if intersects then modify the hit parameter ray.t and return true, otherwise return false
bool intersectRayWithShape(const Sphere& sphere, Ray& ray, HitInfo& hitInfo)
{
    glm::vec3 rayOrigin = ray.origin;
    ray.origin = ray.origin - sphere.center;

    float A = pow(ray.direction.x, 2) + pow(ray.direction.y, 2) + pow(ray.direction.z, 2);
    float B = 2 * (ray.direction.x * ray.origin.x + ray.direction.y * ray.origin.y + ray.direction.z * ray.origin.z);
    float C = pow(ray.origin.x, 2) + pow(ray.origin.y, 2) + pow(ray.origin.z, 2) - pow(sphere.radius, 2);

    ray.origin = rayOrigin;

    float discriminant = pow(B, 2) - (4 * A * C);
    if (discriminant > 0) {
        float d = glm::min((-B + sqrt(discriminant)) / (2 * A), (-B - sqrt(discriminant)) / (2 * A));
        if (d < ray.t) {
            ray.t = d;
            hitInfo.material = sphere.material;
            glm::vec3 selectedPos = ray.origin + ray.direction * ray.t;
            glm::vec3 normalDirection = glm::normalize(selectedPos - sphere.center);
            hitInfo.normal = selectedPos + normalDirection;
        }
        return true;
    }
    else if (discriminant == 0) {
        float d = (-B / (2 * A));
        if (d < ray.t) {
            ray.t = d;
            hitInfo.material = sphere.material;
            glm::vec3 selectedPos = ray.origin + ray.direction * ray.t;
            glm::vec3 normalDirection = glm::normalize(selectedPos - sphere.center);
            hitInfo.normal = selectedPos + normalDirection;
        }
        return true;
    }

    return false;
}

/// Input: an axis-aligned bounding box with the following parameters: minimum coordinates box.lower and maximum coordinates box.upper
/// Output: if intersects then modify the hit parameter ray.t and return true, otherwise return false
bool intersectRayWithShape(const AxisAlignedBox& box, Ray& ray)
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

    if (tin < ray.t) {
        ray.t = tin;
    }

    return true;
}
