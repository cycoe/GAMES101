//
// Created by LEI XU on 5/16/19.
//

#ifndef RAYTRACING_MATERIAL_H
#define RAYTRACING_MATERIAL_H

#include "Vector.hpp"
#include "global.hpp"

inline Vector3f local_to_world(Vector3f const& v, Vector3f const& x, Vector3f const& y, Vector3f const& z) {
  return v.x * x + v.y * y + v.z * z;
}

inline Vector3f world_to_local(Vector3f const& v, Vector3f const& x, Vector3f const& y, Vector3f const& z) {
  float c11 = y.y * z.z - y.z * z.y;
  float c12 = x.z * z.y - x.y * z.z;
  float c13 = x.y * y.z - x.z * y.y;
  float c21 = y.z * z.x - y.x * z.z;
  float c22 = x.x * z.z - x.z * z.x;
  float c23 = x.z * y.x - x.x * y.z;
  float c31 = y.x * z.y - y.y * z.x;
  float c32 = x.y * z.x - x.x * z.y;
  float c33 = x.x * y.y - x.y * y.x;
  float det_inv = 1.f / (x.x * c11 - y.x * c12 + z.x * c13);
  float x_ = (v.x * c11 + v.y * c21 + v.z * c31) * det_inv;
  float y_ = (v.x * c12 + v.y * c22 + v.z * c32) * det_inv;
  float z_ = (v.x * c13 + v.y * c23 + v.z * c33) * det_inv;
  return Vector3f(x_, y_, z_);
}

enum MaterialType {
    DIFFUSE,
    SPECULAR,
    GLASS,
    GLAZE,
    BLOCK,
    OREN_NAYAR,
    GXX
};

class Material{
private:

    // Compute reflection direction
    Vector3f reflect(const Vector3f &I, const Vector3f &N) const
    {
        return I - 2 * dotProduct(I, N) * N;
    }

    // Compute refraction direction using Snell's law
    //
    // We need to handle with care the two possible situations:
    //
    //    - When the ray is inside the object
    //
    //    - When the ray is outside.
    //
    // If the ray is outside, you need to make cosi positive cosi = -N.I
    //
    // If the ray is inside, you need to invert the refractive indices and negate the normal N
    Vector3f refract(const Vector3f &I, const Vector3f &N, const float &ior) const
    {
        float cosi = clamp(-1, 1, dotProduct(I, N));
        float etai = 1, etat = ior;
        Vector3f n = N;
        if (cosi < 0) { cosi = -cosi; } else { std::swap(etai, etat); n= -N; }
        float eta = etai / etat;
        float k = 1 - eta * eta * (1 - cosi * cosi);
        return k < 0 ? 0 : eta * I + (eta * cosi - sqrtf(k)) * n;
    }

    // Compute Fresnel equation
    //
    // \param I is the incident view direction
    //
    // \param N is the normal at the intersection point
    //
    // \param ior is the material refractive index
    //
    // \param[out] kr is the amount of light reflected
    void fresnel(const Vector3f &I, const Vector3f &N, const float &ior, float &kr) const
    {
        float cosi = clamp(-1, 1, dotProduct(I, N));
        float etai = 1, etat = ior;
        if (cosi > 0) {  std::swap(etai, etat); }
        // Compute sini using Snell's law
        float sint = etai / etat * sqrtf(std::max(0.f, 1 - cosi * cosi));
        // Total internal reflection
        if (sint >= 1) {
            kr = 1;
        }
        else {
            float cost = sqrtf(std::max(0.f, 1 - sint * sint));
            cosi = fabsf(cosi);
            float Rs = ((etat * cosi) - (etai * cost)) / ((etat * cosi) + (etai * cost));
            float Rp = ((etai * cosi) - (etat * cost)) / ((etai * cosi) + (etat * cost));
            kr = (Rs * Rs + Rp * Rp) / 2;
        }
        // As a consequence of the conservation of energy, transmittance is given by:
        // kt = 1 - kr;
    }

    Vector3f toWorld(const Vector3f &a, const Vector3f &N){
        Vector3f B, C;
        if (std::fabs(N.x) > std::fabs(N.y)){
            float invLen = 1.0f / std::sqrt(N.x * N.x + N.z * N.z);
            C = Vector3f(N.z * invLen, 0.0f, -N.x *invLen);
        }
        else {
            float invLen = 1.0f / std::sqrt(N.y * N.y + N.z * N.z);
            C = Vector3f(0.0f, N.z * invLen, -N.y *invLen);
        }
        B = crossProduct(C, N);
        return a.x * B + a.y * C + a.z * N;
    }

public:
    MaterialType m_type;
    //Vector3f m_color;
    Vector3f m_emission;
    float ior;
    Vector3f Kd, Ks;
    float specularExponent;
    //Texture tex;

    inline Material(MaterialType t=DIFFUSE, Vector3f e=Vector3f(0,0,0));
    inline virtual MaterialType getType();
    //inline Vector3f getColor();
    inline Vector3f getColorAt(double u, double v);
    inline Vector3f getEmission();
    inline bool hasEmission();

    // sample a ray by Material properties
    inline virtual Vector3f sample(const Vector3f &wi);
    // given a ray, calculate the PdF of this ray
    inline virtual float pdf(const Vector3f &wi, const Vector3f &wo);
    // given a ray, calculate the contribution of this ray
    inline virtual Vector3f eval(const Vector3f& coords, const Vector3f &wi, const Vector3f &wo);

};

Material::Material(MaterialType t, Vector3f e){
    m_type = t;
    //m_color = c;
    m_emission = e;
}

MaterialType Material::getType(){return m_type;}
///Vector3f Material::getColor(){return m_color;}
Vector3f Material::getEmission() {return m_emission;}
bool Material::hasEmission() {
    if (m_emission.norm() > EPSILON) return true;
    else return false;
}

Vector3f Material::getColorAt(double u, double v) {
    return Vector3f();
}


Vector3f Material::sample(const Vector3f &wi){
    switch(m_type){
        case DIFFUSE:
        case BLOCK:
        case OREN_NAYAR:
        {
          // uniform sample on the hemisphere
            float x_1 = get_random_float(), x_2 = get_random_float();
            float z = std::fabs(1.0f - 2.0f * x_1);
            float r = std::sqrt(1.0f - z * z), phi = 2 * M_PI * x_2;
            Vector3f localRay(r*std::cos(phi), r*std::sin(phi), z);
            return localRay;

            break;
        }
        case SPECULAR:
        {
            return Vector3f(0.f);
        }
        case GLASS:
        {
          return Vector3f(0.f);
        }
    }
}

float Material::pdf(const Vector3f &wi, const Vector3f &wo){
  switch(m_type){
  case DIFFUSE:
  case BLOCK:
  case OREN_NAYAR: {
    // uniform sample probability 1 / (2 * PI)
    if (wo.z > 0.f)
      return 0.5f / M_PI;
    else
      return 0.0f;
    break;
  }
  case SPECULAR:
    {
      return 0.f;
      }
    }
}

Vector3f Material::eval(const Vector3f& coords, const Vector3f &wi, const Vector3f &wo){
  switch(m_type){
  case DIFFUSE:
    {
      // calculate the contribution of diffuse   model
      if (wo.z > 0.0f) {
        Vector3f diffuse = Kd / M_PI;
              return diffuse;
      }
            else
                return Vector3f(0.0f);
            break;
        }
        // case BLOCK:
        //   {
        //     // calculate the contribution of diffuse   model
        //     float cosalpha = dotProduct(N, wo);
        //     if (cosalpha > 0.0f) {
        //         Vector3f diffuse = Kd / M_PI;
        //         if (coords.y < 0.001 &&
        //             (coords.x - (int)(coords.x / 22) * 22 > 20.f ||
        //              coords.z - (int)(coords.z / 22) * 22 > 20.f))
        //           diffuse = 0.1 * diffuse;
        //         return diffuse;
        //     }
        //     else
        //         return Vector3f(0.0f);
        //     break;
        // }
      //   case SPECULAR:
      //   {
	    // Vector3f wr = reflect(-wo, N);
	    // if (dotProduct(wr, wi) > 0.99) {
      //           return Ks;
	    // } else {
      //           return Vector3f(0.f);
	    // }
      //   }
      //   case GLASS:
      //   {
      //     float ior = 1.5f;
      //     Vector3f wr = reflect(-wo, N);
      //     Vector3f wt = refract(-wo, N, ior);
      //     float kr = 0.f;
      //     fresnel(-wi, N, ior, kr);

      //     if (dotProduct(wi, wr) > 0.99f) {
      //       return kr * Ks;
      //     } else if (dotProduct(wi, wt) > 0.99f) {
      //       return (1 - kr) * Ks;
      //     } else {
      //         return Vector3f(0.f);
      //       }
      //   }
    case OREN_NAYAR:
    {
      return Vector3f(0.f);
      //float sigma = 20.f;
      ////sigma = sigma * M_PI / 180.f;
      //float sigma2 = sigma * sigma;
      //float A = 1.f - (sigma2 / (2.f * (sigma2 + 0.33f)));
      //float B = 0.45f * sigma2 / (sigma2 + 0.09f);

      //float cos_theta_i = dotProduct(wi, N);
      //float cos_theta_o = dotProduct(wo, N);
      //float sin_theta_i = std::sqrt(1 - cos_theta_i * cos_theta_i);
      //float sin_theta_o = std::sqrt(1 - cos_theta_o * cos_theta_o);

      //float max_cos = 0;
      //if (sin_theta_i > 1e-4 && sin_theta_o > 1e-4) {
      //  Vector3f phi_i = (wi - cos_theta_i * N).normalized();
      //  Vector3f phi_o = (wo - cos_theta_o * N).normalized();
      //  max_cos = std::max(0.f, dotProduct(-phi_i, phi_o));
      //}

      //float sin_alpha, tan_beta;
      //if (std::abs(cos_theta_i) > std::abs(cos_theta_o)) {
      //  sin_alpha = sin_theta_o;
      //  tan_beta = sin_theta_i / std::abs(cos_theta_i);
      //} else {
      //  sin_alpha = sin_theta_i;
      //  tan_beta = sin_theta_o / std::abs(cos_theta_o);
      //}

      //return (A + B * max_cos * sin_alpha * tan_beta) / M_PI * Kd;
    }
    }

    return Vector3f(0.0f);
}

#endif //RAYTRACING_MATERIAL_H
