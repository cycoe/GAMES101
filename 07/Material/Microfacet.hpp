#ifndef __RAY_TRAING_MICROFACET_HPP__
#define __RAY_TRAING_MICROFACET_HPP__

#include "../Vector.hpp"
#include "../global.hpp"
#include "../Material.hpp"
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <random>

inline Vector3f local_to_world(Vector3f const& wh, Vector3f const& N, Vector3f const& X, Vector3f const& Y) {
  return wh.x * X + wh.y * Y + wh.z * N;
}

struct TrowbridgeReitzDistribution {
  TrowbridgeReitzDistribution(float alpha_x, float alpha_y)
    : m_alpha_x(alpha_x), m_alpha_y(alpha_y) { }
  TrowbridgeReitzDistribution(TrowbridgeReitzDistribution const&) = default;
  inline float D(Vector3f const& wh, Vector3f const& N) const;
  inline float G1(Vector3f const& w, Vector3f const& N) const;
  inline float Lambda(Vector3f const& w, Vector3f const& N) const;
  inline float G(Vector3f const& wo, Vector3f const& wi, Vector3f const& N) const;

  float m_alpha_x;
  float m_alpha_y;
};

struct MicrofacetReflection : public Material {
  MicrofacetReflection(Vector3f const &R, float alpha_x, float alpha_y,
                       TrowbridgeReitzDistribution const &distribute)
      : m_R(R), m_alpha_x(alpha_x), m_alpha_y(alpha_y),
        m_distribution(distribute) {}
  // sample a ray by Material properties
  virtual Vector3f sample(const Vector3f &wi, const Vector3f &N) override;
  // given a ray, calculate the PdF of this ray
  virtual float pdf(const Vector3f &wi, const Vector3f &wo, const Vector3f &N) override;
  // given a ray, calculate the contribution of this ray
  virtual Vector3f eval(const Vector3f &coords, const Vector3f &wi,
                       const Vector3f &wo, const Vector3f &N) override;
  virtual MaterialType getType() override { return GXX; };

  Vector3f m_R;
  float m_alpha_x;
  float m_alpha_y;
  TrowbridgeReitzDistribution m_distribution;
};

float TrowbridgeReitzDistribution::D(const Vector3f &wh, Vector3f const& N) const {
  float cos_theta = dotProduct(wh, N);
  float cos_theta_2 = cos_theta * cos_theta;
  float sin_theta_2 = 1 - cos_theta_2;
  float tan_theta_2 = sin_theta_2 / cos_theta_2;

  if (std::isinf(tan_theta_2)) return 0.f;

  Vector3f x, y;
  if (std::abs(m_alpha_x - m_alpha_y) < 0.001f) {
    float e = 1 + tan_theta_2 / (m_alpha_x * m_alpha_y);
    return 1.f / (M_PI * m_alpha_x * m_alpha_x * cos_theta_2 * cos_theta_2 * e * e);
  } else {
    return 0.f;
  }
}

float TrowbridgeReitzDistribution::G1(const Vector3f &w, Vector3f const& N) const {
  return 1 / (1 + this->Lambda(w, N));
}

float TrowbridgeReitzDistribution::Lambda(const Vector3f &w, Vector3f const& N) const {
  float cos_theta = dotProduct(w, N);
  float cos_theta_2 = cos_theta * cos_theta;
  float sin_theta_2 = 1 - cos_theta_2;
  float tan_theta_2 = sin_theta_2 / cos_theta_2;

  if (std::isinf(tan_theta_2)) return 0.f;

  if (std::abs(m_alpha_x - m_alpha_y) < 0.001f) {
    float alpha_tan_theta_2 = m_alpha_x * m_alpha_y * tan_theta_2;
    return (-1 + std::sqrt(1.f + alpha_tan_theta_2)) / 2;
  } else {
    return 0.f;
  }
}

float TrowbridgeReitzDistribution::G(const Vector3f &wo, const Vector3f &wi, Vector3f const& N) const {
  return 1 / (1 + Lambda(wo, N) + Lambda(wi, N));
}

Vector3f MicrofacetReflection::sample(const Vector3f &wi, const Vector3f &N) {
  float ksi_1 = get_random_float();
  float ksi_2 = get_random_float();
  float phi = 2 * M_PI * ksi_1;
  float tan_theta_2 = ksi_2 * m_alpha_x * m_alpha_y / (1 - ksi_2);
  float cos_theta = 1 / std::sqrt(1 + tan_theta_2);
  float sin_theta = std::sqrt(1 - cos_theta * cos_theta);
  Vector3f wh(sin_theta * std::cos(phi), sin_theta * std::sin(phi), cos_theta);

  Vector3f X;
  Vector3f Y;
  if (std::fabs(N.x) > std::fabs(N.y)) {
    float invLen = 1.0f / std::sqrt(N.x * N.x + N.z * N.z);
    X = Vector3f(N.z * invLen, 0.0f, -N.x * invLen);
    Y = crossProduct(N, X);
  } else {
    float invLen = 1.0f / std::sqrt(N.y * N.y + N.z * N.z);
    Y = Vector3f(0.0f, N.z * invLen, -N.y * invLen);
    X = crossProduct(Y, N);
  }

  wh = local_to_world(wh, N, X, Y).normalized();
  return 2 * dotProduct(wh, wi) * wh - wi;
}

float MicrofacetReflection::pdf(const Vector3f &wi, const Vector3f &wo, const Vector3f &N) {
  Vector3f wh = (wi + wo).normalized();
  float cos_theta = dotProduct(wh, N);
  if (cos_theta < 0.001f) return 0.f;
  float cos_theta_2 = cos_theta * cos_theta;
  float sin_theta_2 = 1 - cos_theta_2;
  float tan_theta_2 = sin_theta_2 / cos_theta_2;
  float e = 1 + tan_theta_2 / (m_alpha_x * m_alpha_y);
  float pdf = 1 / (M_PI * m_alpha_x * m_alpha_y * cos_theta * cos_theta_2 * e * e);
  return pdf;
}

Vector3f MicrofacetReflection::eval(const Vector3f &coords, const Vector3f &wi, const Vector3f &wo, const Vector3f &N) {
  float cos_theta_i = std::abs(dotProduct(wi, N));
  float cos_theta_o = std::abs(dotProduct(wo, N));
  Vector3f wh = wi + wo;

  if (cos_theta_i < 0.001f || cos_theta_o <  0.001f) return Vector3f(0.f);
  if (wh.norm() < 0.001f) return Vector3f(0.f);

  wh = wh.normalized();
  float d = m_distribution.D(wh, N);
  float g = m_distribution.G(wo, wi, N);
  float f = 0.98f;
  Vector3f R = d * g * f / (4 * cos_theta_i * cos_theta_o) * m_R;

  return R;
}

#endif  // __RAY_TRAING_MICROFACET_HPP__