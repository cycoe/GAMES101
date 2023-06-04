//
// Created by Göksu Güvendiren on 2019-05-14.
//

#include "Scene.hpp"
#include <cassert>

void Scene::buildBVH() {
    printf(" - Generating BVH...\n\n");
    this->bvh = new BVHAccel(objects, 1, BVHAccel::SplitMethod::NAIVE);
}

Intersection Scene::intersect(const Ray &ray) const
{
    return this->bvh->Intersect(ray);
}

void Scene::sampleLight(Intersection &pos, float &pdf) const
{
    float emit_area_sum = 0;
    for (uint32_t k = 0; k < objects.size(); ++k) {
        if (objects[k]->hasEmit()){
            emit_area_sum += objects[k]->getArea();
        }
    }
    float p = get_random_float() * emit_area_sum;//[0~1]*13650
    emit_area_sum = 0;
    for (uint32_t k = 0; k < objects.size(); ++k) {
        if (objects[k]->hasEmit()){
            emit_area_sum += objects[k]->getArea();
            if (p <= emit_area_sum){//random get the first area > p light,return
                objects[k]->Sample(pos, pdf);
                break;
            }
        }
    }
}



bool Scene::trace(
        const Ray &ray,
        const std::vector<Object*> &objects,
        float &tNear, uint32_t &index, Object **hitObject)
{
    *hitObject = nullptr;
    for (uint32_t k = 0; k < objects.size(); ++k) {
        float tNearK = kInfinity;
        uint32_t indexK;
        Vector2f uvK;
        if (objects[k]->intersect(ray, tNearK, indexK) && tNearK < tNear) {
            *hitObject = objects[k];
            tNear = tNearK;
            index = indexK;
        }
    }


    return (*hitObject != nullptr);
}


Vector3f Scene::shade(Intersection& hit_obj, Vector3f wo, int depth) const
{
    /// 如果打到的物体是光源，则直接返回光源的强度
    if (hit_obj.m->hasEmission())
    {
        return hit_obj.m->getEmission();
    }

    const float epsilon = 0.0005f;
    // 直接光照贡献
    Vector3f Lo_dir;
    {
		float light_pdf;
		Intersection hit_light;
		sampleLight(hit_light, light_pdf);
		Vector3f obj2Light = hit_light.coords - hit_obj.coords;
        Vector3f obj2LightDir = obj2Light.normalized();

        // 检查光线是否被遮挡
        auto t = intersect(Ray(hit_obj.coords, obj2LightDir));
        if (t.distance - obj2Light.norm() > -epsilon)
        {
          Vector3f f_r = hit_obj.m->eval(hit_obj.coords, obj2LightDir, wo, hit_obj.normal);
          float r2 = dotProduct(obj2Light, obj2Light);
          float cosA = std::abs(dotProduct(hit_obj.normal,obj2LightDir));
          if (hit_obj.m->getType() == SPECULAR || hit_obj.m->getType() == GLASS)
		    cosA = 1.f;
            float cosB = std::max(.0f, dotProduct(hit_light.normal,-obj2LightDir));
			Lo_dir = hit_light.emit * f_r * cosA * cosB / r2 / light_pdf;
        }
    }

    // std::cout << depth << " | lo_dir: " << Lo_dir.x << "," << Lo_dir.y << "," << Lo_dir.z << std::endl;

    // 间接光照贡献
    Vector3f Lo_indir;
    {
		if (get_random_float() < RussianRoulette)
		{
			Vector3f dir2NextObj = hit_obj.m->sample(wo, hit_obj.normal).normalized();
            float pdf = hit_obj.m->pdf(wo, dir2NextObj, hit_obj.normal);
            //std::cout << "pdf in scene " << pdf << std::endl;
            if (pdf > epsilon)
            {
                Intersection nextObj = intersect(Ray(hit_obj.coords, dir2NextObj));
				if (nextObj.happened && !nextObj.m->hasEmission())
				{
                  Vector3f f_r = hit_obj.m->eval(hit_obj.coords, dir2NextObj, wo, hit_obj.normal); //BRDF
                  float cos = std::abs(dotProduct(dir2NextObj, hit_obj.normal));
                  if (hit_obj.m->getType() == SPECULAR || hit_obj.m->getType() == GLASS)
						cos = 1.f;
                    Vector3f decay(1.f);
                    if (dotProduct(hit_obj.normal, dir2NextObj) < 0.f &&
                        dotProduct(nextObj.normal, dir2NextObj) > 0.f) {
                      Vector3f absorb(0.01, 0.002, 0.003);
                      float distance = (nextObj.coords - hit_obj.coords).norm();
                      Vector3f coff = -distance * absorb;
                      decay = Vector3f(std::exp(coff.x), std::exp(coff.y), std::exp(coff.z));
                    }

                    Vector3f mid_shade = shade(nextObj, -dir2NextObj, depth + 1);
                    Lo_indir = mid_shade * f_r * cos * decay / pdf / RussianRoulette;
                }
            }
		}
    }

    return Lo_dir + Lo_indir;
}



// Implementation of Path Tracing
Vector3f Scene::castRay(const Ray &ray, int depth) const
{
	auto hitObj = intersect(ray);
	if (!hitObj.happened) return {};
    return shade(hitObj,-ray.direction, depth);
}
