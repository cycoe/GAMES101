#ifndef __RAY_TRAING_MICROFACET_HPP__
#define __RAY_TRAING_MICROFACET_HPP__

#include "../Vector.hpp"
#include "../global.hpp"
#include "../Material.hpp"
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <random>


inline Vector3f reflect(const Vector3f &I, const Vector3f &N)
{
    return I - 2 * dotProduct(I, N) * N;
}

inline bool refract(const Vector3f &wi, const Vector3f &n, float eta, Vector3f& wt) {
    // Compute $\cos \theta_\roman{t}$ using Snell's law
    float cosThetaI = dotProduct(n, wi);
    float sin2ThetaI = std::max(0.f, 1 - cosThetaI * cosThetaI);
    float sin2ThetaT = eta * eta * sin2ThetaI;

    // Handle total internal reflection for transmission
    if (sin2ThetaT >= 1) return false;
    float cosThetaT = std::sqrt(1 - sin2ThetaT);
    wt = eta * -wi + (eta * cosThetaI - cosThetaT) * n;
    return true;
}

inline float FrDielectric(float cosThetaI, float etaI, float etaT) {
    cosThetaI = clamp(-1, 1, cosThetaI);
    // Potentially swap indices of refraction
    bool entering = cosThetaI > 0.f;
    if (!entering) {
        std::swap(etaI, etaT);
        cosThetaI = std::abs(cosThetaI);
    }

    // Compute _cosThetaT_ using Snell's law
    float sinThetaI = std::sqrt(std::max(0.f, 1 - cosThetaI * cosThetaI));
    float sinThetaT = etaI / etaT * sinThetaI;

    // Handle total internal reflection
    if (sinThetaT >= 1) return 1;
    float cosThetaT = std::sqrt(std::max(0.f, 1 - sinThetaT * sinThetaT));
    float Rparl = ((etaT * cosThetaI) - (etaI * cosThetaT)) /
                  ((etaT * cosThetaI) + (etaI * cosThetaT));
    float Rperp = ((etaI * cosThetaI) - (etaT * cosThetaT)) /
                  ((etaI * cosThetaI) + (etaT * cosThetaT));
    return (Rparl * Rparl + Rperp * Rperp) / 2;
}

inline float get_cos_theta(Vector3f const& w) {
    return w.z;
}

inline float get_cos_theta_2(Vector3f const& w) {
    return w.z * w.z;
}

inline float get_sin_theta_2(Vector3f const& w) {
    return std::max(0.f, 1.f - get_cos_theta_2(w));
}

inline float get_sin_theta(Vector3f const& w) {
    return std::sqrt(get_sin_theta_2(w));
}

inline float get_cos_phi(Vector3f const& w) {
    float sin_theta = get_sin_theta(w);
    return (sin_theta == 0) ? 1.f : clamp(-1.f, 1.f, w.x / sin_theta);
}

inline float get_sin_phi(Vector3f const& w) {
    float sin_theta = get_sin_theta(w);
    return (sin_theta == 0) ? 0.f : clamp(-1.f, 1.f, w.y / sin_theta);
}

struct TrowbridgeReitzDistribution {
  TrowbridgeReitzDistribution(float alpha_x, float alpha_y, bool sample_visible_area=true)
    : m_alpha_x(alpha_x), m_alpha_y(alpha_y), m_sample_visible_area(sample_visible_area) { }
  TrowbridgeReitzDistribution(TrowbridgeReitzDistribution const&) = default;
  inline float D(Vector3f const& wh) const;
  inline float G1(Vector3f const& w) const;
  inline float Lambda(Vector3f const& w) const;
  inline float G(Vector3f const& wo, Vector3f const& wi) const;
  inline Vector3f Sample_wh(Vector3f const& wo, float r1, float r2) const;
  inline float Pdf(Vector3f const& wo, Vector3f const& wh) const;

  float m_alpha_x;
  float m_alpha_y;
  bool m_sample_visible_area;
};

float TrowbridgeReitzDistribution::D(const Vector3f &wh) const {
  float cos_theta = wh.z;
  float cos_theta_2 = cos_theta * cos_theta;
  float sin_theta_2 = 1 - cos_theta_2;
  float tan_theta_2 = sin_theta_2 / cos_theta_2;

  if (std::isinf(tan_theta_2)) return 0.f;

  if (std::abs(m_alpha_x - m_alpha_y) < 0.001f) {
    float e = 1 + tan_theta_2 / (m_alpha_x * m_alpha_y);
    return 1.f / (M_PI * m_alpha_x * m_alpha_y * cos_theta_2 * cos_theta_2 * e * e);
  } else {
    float cos_phi = get_cos_phi(wh);
    float cos_phi_2 = cos_phi * cos_phi;
    float sin_phi_2 = 1 - cos_phi_2;
    float e = 1 + tan_theta_2 * (cos_phi_2 / (m_alpha_x * m_alpha_x) +
                                 sin_phi_2 / (m_alpha_y * m_alpha_y));
    return 1.f / (M_PI * m_alpha_x * m_alpha_y * cos_theta_2 * cos_theta_2 * e * e);
  }
}

float TrowbridgeReitzDistribution::G1(const Vector3f &w) const {
  return 1 / (1 + this->Lambda(w));
}

float TrowbridgeReitzDistribution::Lambda(const Vector3f &w) const {
  float cos_theta = w.z;
  float cos_theta_2 = cos_theta * cos_theta;
  float sin_theta_2 = 1 - cos_theta_2;
  float tan_theta_2 = sin_theta_2 / cos_theta_2;

  if (std::isinf(tan_theta_2)) return 0.f;

  if (std::abs(m_alpha_x - m_alpha_y) < 0.001f) {
    float alpha_tan_theta_2 = m_alpha_x * m_alpha_y * tan_theta_2;
    return (-1 + std::sqrt(1.f + alpha_tan_theta_2)) / 2;
  } else {
    float cos_phi = get_cos_phi(w);
    float cos_phi_2 = cos_phi * cos_phi;
    float sin_phi_2 = 1 - cos_phi_2;
    float alpha_2 = cos_phi_2 * m_alpha_x * m_alpha_x + sin_phi_2 * m_alpha_y * m_alpha_y;
    return (-1 + std::sqrt(1.f + alpha_2 * tan_theta_2)) / 2;
  }
}

float TrowbridgeReitzDistribution::G(const Vector3f &wo, const Vector3f &wi) const {
  return 1 / (1 + Lambda(wo) + Lambda(wi));
  //return G1(wo, N) * G1(wi, N);
}

inline void TrowbridgeReitzSample11(float cos_theta, float r1, float r2, float* slope_x, float* slope_y) {
  /// special case (normal incidence)
  if (cos_theta > .9999) {
    float r = std::sqrt(r1 / (1 - r1));
    float phi = 2 * M_PI * r2;
    *slope_x = r * std::cos(phi);
    *slope_y = r * std::sin(phi);
    return;
  }

  float sin_theta = std::sqrt(std::max(0.f, 1 - cos_theta * cos_theta));
  float tan_theta = sin_theta / cos_theta;
  float a = 1.f / tan_theta;
  float G1 = 2 / (1.f + std::sqrt(1.f + 1.f / (a * a)));

  /// sample slope_x
  float A = 2 * r1 / G1 - 1;
  float tmp = 1.f / (A * A - 1.f);
  if (tmp > 1e10) tmp = 1e10;
  float B = tan_theta;
  float D = std::sqrt(std::max(0.f, B * B * tmp * tmp - (A * A - B * B) * tmp));
  float slope_x_1 = B * tmp - D;
  float slope_x_2 = B * tmp + D;
  *slope_x = (A < 0 || slope_x_2 > 1.f / tan_theta) ? slope_x_1 : slope_x_2;

  /// sample slope_y
  float S;
  if (r2 > 0.5f) {
    S = 1.f;
    r2 = 2.f * (r2 - 0.5f);
  } else {
    S = -1.f;
    r2 = 2.f * (0.5f - r2);
  }
  float z = (r2 * (r2 * (r2 * 0.27385f - 0.73369f) + 0.46341f)) /
            (r2 * (r2 * (r2 * 0.093073f + 0.309420f) - 1.000000f) + 0.597999f);
  *slope_y = S * z * std::sqrt(1.f + *slope_x * * slope_y);
}

inline Vector3f TrowbridgeReitzSample(Vector3f const& wi, float alpha_x, float alpha_y, float r1, float r2) {
  /// 1. stretch wi
  Vector3f wi_stretched = Vector3f(alpha_x * wi.x, alpha_y * wi.y, wi.z).normalized();

  /// 2. simulate P22_{wi}(x_slope, y_slope, 1, 1)
  float slope_x, slope_y;
  TrowbridgeReitzSample11(wi_stretched.z, r1, r2, &slope_x, &slope_y);

  /// 3. rotate
  float cos_theta = wi_stretched.z;
  float r = std::sqrt(std::max(0.f, 1 - cos_theta * cos_theta));
  float tmp = wi_stretched.x / r * slope_x - wi_stretched.y / r * slope_y;
  slope_y = wi_stretched.y / r * slope_x + wi_stretched.x / r * slope_y;
  slope_x = tmp;

  /// 4. unstretch
  slope_x = alpha_x * slope_x;
  slope_y = alpha_y * slope_y;

  /// 3. compute normal
  return Vector3f(-slope_x, -slope_y, 1.f).normalized();
}

Vector3f TrowbridgeReitzDistribution::Sample_wh(const Vector3f &wo, float r1, float r2) const {
  Vector3f wh;
  if (!m_sample_visible_area) {
    float cos_theta = 0.f;
    float phi = 2 * M_PI * r2;
    if (std::abs(m_alpha_x - m_alpha_y) < 0.001f) {
      float tan_theta_2 = m_alpha_x * m_alpha_y * r1 / (1.f - r1);
      cos_theta = 1 / std::sqrt(1 + tan_theta_2);
    } else {
      /// TODO
    }
    float sin_theta = std::sqrt(std::max(0.f, 1 - cos_theta * cos_theta));
    wh = Vector3f(sin_theta * std::cos(phi), sin_theta * std::sin(phi), cos_theta);
    if (wh.z * wo.z < 0.f) wh = -wh;
  } else {
    bool flip = wo.z < 0;
    wh = TrowbridgeReitzSample(flip ? -wo : wo, m_alpha_x, m_alpha_y, r1, r2);
    if (flip) wh = -wh;
  }
  return wh;
}

float TrowbridgeReitzDistribution::Pdf(const Vector3f &wo, const Vector3f &wh) const {
  if (m_sample_visible_area) {
    return this->D(wh) * this->G1(wo) * std::abs(dotProduct(wo, wh)) / std::abs(wo.z);
  } else {
    return this->D(wh) * std::abs(wh.z);
  }
}

class BxDF {
public:
    virtual Vector3f Sample_wh(Vector3f const& wo, float r1, float r2) const = 0;
    virtual Vector3f f(Vector3f const& wo, Vector3f const& wi) const = 0;
    virtual Vector3f Sample_f(Vector3f const& wo, Vector3f const& wh, Vector3f& wi, float& pdf) const = 0;
    virtual float Pdf(Vector3f const& wo, Vector3f const& wi) const = 0;
};

class MicrofacetReflection : public BxDF {
public:
    MicrofacetReflection(Vector3f const& R, TrowbridgeReitzDistribution const& distribution)
    : m_R(R), m_distribution(distribution) { }
    inline virtual Vector3f Sample_wh(Vector3f const& wo, float r1, float r2) const override;
    inline virtual Vector3f f(Vector3f const& wo, Vector3f const& wi) const override;
    inline virtual Vector3f Sample_f(Vector3f const& wo, Vector3f const& wh, Vector3f& wi, float& pdf) const override;
    inline virtual float Pdf(Vector3f const& wo, Vector3f const& wi) const override;

protected:
    Vector3f m_R;
    TrowbridgeReitzDistribution m_distribution;
};

Vector3f MicrofacetReflection::Sample_wh(Vector3f const& wo, float r1, float r2) const {
  return m_distribution.Sample_wh(wo, r1, r2);
}

Vector3f MicrofacetReflection::f(const Vector3f &wo, const Vector3f &wi) const {
  float cos_theta_i = std::abs(wi.z);
  float cos_theta_o = std::abs(wo.z);
  Vector3f wh = wi + wo;

  if (cos_theta_i < 0.001f || cos_theta_o <  0.001f) return Vector3f(0.f);
  if (wh.norm() < 0.001f) return Vector3f(0.f);

  wh = wh.normalized();
  float d = m_distribution.D(wh);
  float g = m_distribution.G(wi, wo);
  if (wh.z < 0) wh = -wh;
  float m_eta_a = 1.f;
  float m_eta_b = 1.5f;
  float f = FrDielectric(dotProduct(wo, wh), m_eta_a, m_eta_b);
  Vector3f R = d * g * f / (4 * cos_theta_i * cos_theta_o) * m_R;

  return R;
}

Vector3f MicrofacetReflection::Sample_f(Vector3f const& wo, Vector3f const& wh, Vector3f& wi, float& pdf) const {
    /// sample microfacet orientation wh and reflected direction wi
    if (wo.z == 0) {
        pdf = 0.f;
        return Vector3f(0.f);
    }

    if (dotProduct(wo, wh) < 0) {
        pdf = 0.f;
        return Vector3f(0.f);
    }

    wi = reflect(-wo, wh);
    if (wo.z * wi.z < 0) {
        pdf = 0.f;
        return Vector3f(0.f);
    }

    pdf = m_distribution.Pdf(wo, wh) / (4 * dotProduct(wo, wh));
    return this->f(wo, wi);
}

float MicrofacetReflection::Pdf(const Vector3f &wo, const Vector3f &wi) const {
    if (wo.z * wi.z < 0) return 0.f;
    Vector3f wh = (wo + wi).normalized();
    return m_distribution.Pdf(wo, wh) / (4 * dotProduct(wo, wh));
}

class MicrofacetTransmission : public BxDF {
public:
    MicrofacetTransmission(Vector3f const& T, TrowbridgeReitzDistribution const& distribution)
    : m_eta_a(1), m_eta_b(1.5), m_T(T), m_distribution(distribution) { }
    inline virtual Vector3f Sample_wh(Vector3f const& wo, float r1, float r2) const override;
    inline virtual Vector3f f(Vector3f const& wo, Vector3f const& wi) const override;
    inline virtual Vector3f Sample_f(Vector3f const& wo, Vector3f const& wh, Vector3f& wi, float& pdf) const override;
    inline virtual float Pdf(Vector3f const& wo, Vector3f const& wi) const override;

protected:
    float m_eta_a;
    float m_eta_b;
    Vector3f m_T;
    TrowbridgeReitzDistribution m_distribution;
};

Vector3f MicrofacetTransmission::Sample_wh(Vector3f const& wo, float r1, float r2) const {
  return m_distribution.Sample_wh(wo, r1, r2);
}

Vector3f MicrofacetTransmission::f(const Vector3f &wo, const Vector3f &wi) const {
  if (wo.z * wi.z > 0) return Vector3f(0.f);

  float cos_theta_i = wi.z;
  float cos_theta_o = wo.z;
  if (std::abs(cos_theta_i) < 0.001f || std::abs(cos_theta_o) <  0.001f) return Vector3f(0.f);

  /// compute wh from wo and wi
  float eta = cos_theta_o > 0 ? m_eta_b / m_eta_a : m_eta_a / m_eta_b;
  Vector3f wh = (wo + eta * wi).normalized();
  if (wh.z < 0) wh = -wh;

  /// same side
  if (dotProduct(wo, wh) * dotProduct(wi, wh) > 0) return Vector3f(0.f);

  float sqrt_denom = dotProduct(wo, wh) + eta * dotProduct(wi, wh);
  float factor = true ? 1 / eta : 1;
  float d = m_distribution.D(wh);
  float g = m_distribution.G(wi, wo);
  float f = FrDielectric(dotProduct(wo, wh), m_eta_a, m_eta_b);
  Vector3f T = (1 - f) * std::abs(d * g * eta * eta * std::abs(dotProduct(wi, wh)) *
                std::abs(dotProduct(wo, wh)) * factor * factor /
                (cos_theta_i * cos_theta_o * sqrt_denom * sqrt_denom)) * m_T;

  return T;
}

Vector3f MicrofacetTransmission::Sample_f(Vector3f const& wo, Vector3f const& wh, Vector3f& wi, float& pdf) const {
    /// sample microfacet orientation wh and reflected direction wi
    if (wo.z == 0) {
        pdf = 0.f;
        return Vector3f(0.f);
    }

    if (dotProduct(wo, wh) < 0) {
        pdf = 0.f;
        return Vector3f(0.f);
    }

    float eta = wo.z > 0 ? m_eta_a / m_eta_b : m_eta_b / m_eta_a;
    if (!refract(wo, wh, eta, wi)) {
        pdf = 0.f;
        return Vector3f(0.f);
    }
    pdf = this->Pdf(wo, wi);
    return this->f(wo, wi);
}

float MicrofacetTransmission::Pdf(const Vector3f &wo, const Vector3f &wi) const {
    if (wo.z * wi.z > 0) return 0.f;
    float eta = wo.z > 0 ? m_eta_b / m_eta_a : m_eta_a / m_eta_b;
    Vector3f wh = (wo + eta * wi).normalized();

    if (dotProduct(wo, wh) * dotProduct(wi, wh) > 0) return 0.f;

    float sqrt_denom = dotProduct(wo, wh) + eta * dotProduct(wi, wh);
    float dwh_dwi = std::abs(eta * eta * dotProduct(wi, wh) / (sqrt_denom * sqrt_denom));
    return m_distribution.Pdf(wo, wh) * dwh_dwi;
}

class BSDF : public Material {
public:
    inline BSDF(Vector3f const& R, Vector3f const& T, float alpha_x, float alpha_y);
    inline virtual Vector3f sample(const Vector3f &wi) override;
    inline virtual float pdf(const Vector3f &wi, const Vector3f &wo) override;
    inline virtual Vector3f eval(const Vector3f& coords, const Vector3f &wi, const Vector3f &wo) override;

protected:
    Vector3f m_R;
    Vector3f m_T;
    float m_alpha_x;
    float m_alpha_y;
    BxDF* m_reflection;
    BxDF* m_transmission;
};

BSDF::BSDF(Vector3f const &R, Vector3f const &T, float alpha_x, float alpha_y)
    : m_R(R), m_T(T), m_alpha_x(alpha_x), m_alpha_y(alpha_y)
{
    TrowbridgeReitzDistribution td(m_alpha_x, m_alpha_y);
    m_reflection = new MicrofacetReflection(m_R, td);
    m_transmission = new MicrofacetTransmission(m_T, td);
}

Vector3f BSDF::sample(const Vector3f &wi) {
    float r1 = get_random_float();
    float r2 = get_random_float();
    float rt = get_random_float();

    Vector3f wh = m_reflection->Sample_wh(wi, r1, r2);
    float m_eta_a = 1.f;
    float m_eta_b = 1.5f;
    float f = FrDielectric(dotProduct(wi, wh), m_eta_a, m_eta_b);
    Vector3f wo;
    float pdf;
    if (rt < f) {
      (void)m_reflection->Sample_f(wi, wh, wo, pdf);
    } else {
      (void)m_transmission->Sample_f(wi, wh, wo, pdf);
    }
    return wo;
}

float BSDF::pdf(const Vector3f &wi, const Vector3f &wo) {
    float m_eta_a = 1.f;
    float m_eta_b = 1.5f;
    Vector3f wh;
    if (wi.z * wo.z > 0) {
      wh = (wi + wo).normalized();
    } else {
      float eta = wo.z > 0 ? m_eta_b / m_eta_a : m_eta_a / m_eta_b;
      wh = (wo + eta * wi).normalized();
    }
    if (wh.z < 0) wh = -wh;
    float f = FrDielectric(dotProduct(wi, wh), m_eta_a, m_eta_b);
    if (wi.z * wo.z > 0) {
      return m_reflection->Pdf(wi, wo) * f;
    } else {
      return m_transmission->Pdf(wi, wo) * (1 - f);
    }
}

Vector3f BSDF::eval(const Vector3f &coords, const Vector3f &wi, const Vector3f &wo) {
    if (wi.z * wo.z > 0) {
      return m_reflection->f(wo, wi);
    } else {
      return m_transmission->f(wo, wi);
    }
}

#endif  // __RAY_TRAING_MICROFACET_HPP__