/*
    This file is part of TinyRender, an educative rendering system.

    Designed for ECSE 446/546 Realistic/Advanced Image Synthesis.
    Derek Nowrouzezahrai, McGill University.
*/

#pragma once

#include <tiny_obj_loader.h>
#define RAY_EPS_CV 1e-5 // Use when setting min and max dist for ray in control variates code
TR_NAMESPACE_BEGIN

/**
 * Direct illumination integrator for polygonal light sources
 * Follows Arvo '94.
 */
struct PolygonalIntegrator : Integrator {

	float m_alpha;             // Control variates "strength"
	size_t m_visSamples;       // # of samples to estimate h - alpha*g
	bool m_traceShadows;       // Trace shadows or not
	EPolygonalMethod m_method; // Method to use (Arvo, or control variates)

	std::vector<std::vector<v3f>> m_triangles; // Data structure to store triangles

    explicit PolygonalIntegrator(const Scene& scene) : Integrator(scene) {
        m_alpha = scene.config.integratorSettings.poly.alpha;
        m_visSamples = scene.config.integratorSettings.poly.visSamples;
        m_traceShadows = scene.config.integratorSettings.poly.traceShadows;
        m_method = scene.config.integratorSettings.poly.method;

		/**
		 * 1) Get # of triangles on emitter
		 * 2) Store vertices in m_triangles
		 */
// TODO(A4): Implement this
    }

    /// Reflect
    inline v3f reflect(const v3f& d) const {
        return v3f(-d.x, -d.y, d.z);
    }

    /**
     * === PHONG BONUS ONLY ===
     * Compute the following integral:
     *    T(a, b, n, x) = \int_0^x [a \cos(\theta) + b \sin(\theta)]ˆn d\theta
     * Uses a recurrent relation (see Snyder's note, 1996)
     *
     * Series function:
     *    T_sum(a, b, n, x) = \sum_{i=0}ˆ{(n-1)/2} T(a, b, 2i+1, x)
     * assuming n is _odd_
     */
    float cosineSinePowerIntegralSum(float a, float b, int exp, float theta) const {
        if (exp % 2 == 0) exp += 1; // Make exponent odd
        float Tsum = 0.f;

		// Implementing this function may be useful if you attempt the bonus

        // TODO(A4): Implement this

        return Tsum;
    }

    /**
     * Compute edge (v1--v2) contribution
	 * The exp term is only needed if you attempt the bonus, otherwise, you can ignore it
     */
    float getEdgeContrib(const v3f& v1, const v3f& v2, const SurfaceInteraction& i, int exp = 0) const {
        float contrib = 0.f;

        // TODO(A4): Implement this

        return contrib;
    }
	   

    /// Direct illumination using Arvo '94 analytic solution for polygonal lights
    v3f renderAnalytic(const Ray& ray, Sampler& sampler) const {
        v3f Lr(0.f);

        // TODO(A4): Implement this

        return Lr;
    }

    /**
     * Stand-alone estimator for h - alpha*g (with primary ray)
     * Trace a primary ray, check for emitter hit, and then call `estimateVisDiff()`
     * Used by polygonal render pass
     */
    v3f estimateVisDiffRealTime(const Ray& ray, Sampler& sampler, const Emitter& em) {
        v3f D(0.f);

        SurfaceInteraction hit;
        if (!scene.bvh->intersect(ray, hit)) return D;

        const BSDF* bsdf = getBSDF(hit);
        if (bsdf->isEmissive()) return D;

        hit.wi = v3f(0, 0, 1); // Trick to get 1/pi * albedo without cosine term
        D = estimateVisDiff(sampler, hit, em);

        return D;
    }

    /// Stand-alone estimator for h - alpha*g (without primary ray)
	/// Use RAY_EPS_CV when setting min and max dist for shadow ray
    v3f estimateVisDiff(Sampler& sampler, SurfaceInteraction& i, const Emitter& em) const {
        v3f sum(0.f);

        // TODO(A4): Implement this

        return sum;
    }

    /// Control variates using Arvo '94 for direct illumination; ray trace shadows
	
    v3f renderControlVariates(const Ray& ray, Sampler& sampler) const {
        v3f Lr(0.f);

        // TODO(A4): Implement this

        return Lr;
    }

    /// Direct illumination using surface area sampling
    v3f renderArea(const Ray& ray, Sampler& sampler) const {
        v3f Lr(0.f);

        // TODO(A4): Implement this

        return Lr;
    }

    /// Branch to corresponding method
    v3f render(const Ray& ray, Sampler& sampler) const override {
        switch (m_method) {
            case EPolygonalMethod::ESurfaceArea:
                return PolygonalIntegrator::renderArea(ray, sampler);
                break;
            case EPolygonalMethod::EControlVariates:
                return PolygonalIntegrator::renderControlVariates(ray, sampler);
                break;
            default:
                return PolygonalIntegrator::renderAnalytic(ray, sampler);
                break;
        }
    }

};

TR_NAMESPACE_END