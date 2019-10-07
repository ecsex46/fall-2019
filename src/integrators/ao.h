/*
    This file is part of TinyRender, an educative rendering system.

    Designed for ECSE 446/546 Realistic/Advanced Image Synthesis.
    Derek Nowrouzezahrai, McGill University.
*/

#pragma once

TR_NAMESPACE_BEGIN

/**
 * Ambient occlusion integrator
 */
struct AOIntegrator : Integrator {

	// Use this in your switch statement to select the sampling type 
	ESamplingType m_samplingStrategy;

    explicit AOIntegrator(const Scene& scene) : Integrator(scene) { 
		m_samplingStrategy = scene.config.integratorSettings.ao.sampling_type;
	}

    v3f render(const Ray& ray, Sampler& sampler) const override {
        v3f Li(0.f);
		
		/*
		Use the m_sampling_type variable to set wi and the corresponding pdf 
		appropriately for sphere, hemisphere, or cosine sampling.

		You can use a switch statement or an if/else block.

		The m_sampling_type variable is an enum. The different values of the enum 
		can be accessed through:
		ESamplingType::ESpherical
		ESamplingType::EHemispherical
		ESamplingType::ECosineHemispherical
		*/
		
        // TODO(A3): Implement this

        return Li;
    }
};

TR_NAMESPACE_END