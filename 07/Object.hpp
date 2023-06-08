//
// Created by LEI XU on 5/13/19.
//
#pragma once
#ifndef RAYTRACING_OBJECT_H
#define RAYTRACING_OBJECT_H

#include "Vector.hpp"
#include "global.hpp"
#include "Bounds3.hpp"
#include "Ray.hpp"
#include "Intersection.hpp"

class Object
{
public:
    Object() {}
    virtual ~Object() {}
    virtual bool intersect(const Ray& ray) = 0;
    virtual bool intersect(const Ray& ray, float &, uint32_t &) const = 0;
    virtual Intersection getIntersection(Ray _ray) = 0;
    virtual void getSurfaceProperties(const Vector3f &, const Vector3f &, const uint32_t &, const Vector2f &, Vector3f &, Vector2f &) const = 0;
    virtual Vector3f evalDiffuseColor(const Vector2f &) const =0;
    virtual Bounds3 getBounds()=0;
    virtual float getArea()=0;
    virtual void Sample(Intersection &pos, float &pdf)=0;
    virtual bool hasEmit()=0;
    virtual void get_local_coords(Vector3f const& p, Vector3f const& z, Vector3f& x, Vector3f& y) const {
        if (std::fabs(z.x) > std::fabs(z.y)){
            float invLen = 1.0f / std::sqrt(z.x * z.x + z.z * z.z);
            x = Vector3f(z.z * invLen, 0.0f, -z.x *invLen);
            y = crossProduct(z, x);
        }
        else {
            float invLen = 1.0f / std::sqrt(z.y * z.y + z.z * z.z);
            y = Vector3f(0.0f, z.z * invLen, -z.y *invLen);
            x = crossProduct(y, z);
        }
    }
};



#endif //RAYTRACING_OBJECT_H
