/*
    This file is part of TinyRender, an educative rendering system.

    Designed for ECSE 446/546 Realistic/Advanced Image Synthesis.
    Derek Nowrouzezahrai, McGill University.
*/

#pragma once

#include <core/core.h>
#include <kdtree.h>
#include "direct.h"

TR_NAMESPACE_BEGIN

/* This is the photon
 * The power is not compressed so the
 * size is 28 bytes
*/
typedef struct Photon {
    v3f pos;                  // photon position
    v3f dir;      // incoming direction
    v3f power;   // photon power (uncompressed)
    v3f n;
} Photon;


/**
 * Photon mapping integrator
 */
struct PPMIntegrator : Integrator {

    std::vector<Photon> m_photonMap;

    typedef unsigned int PhotonMapIdx;
    typedef GenericKDTreeNode<v3f, PhotonMapIdx> PhotonKDTreeNode;
    PointKDTree<PhotonKDTreeNode> m_KDTree;

    //1st pass
    int m_photonCount;
    float m_photonRrProb;
    int m_photonRrDepth;

    //2nd pass
    int m_emittedPhotonCount;
    bool m_usePhotonsForDirect;
    float m_radiusSearch;
    int m_nbPhotonsSearch;
    bool m_useFinalGather;
    int m_nbFinalGather;

    std::unique_ptr<DirectIntegrator> m_directIntegrator;

    explicit PPMIntegrator(const Scene& scene) : Integrator(scene)
    {
        m_directIntegrator = std::unique_ptr<DirectIntegrator>(new DirectIntegrator(scene));

        m_directIntegrator->m_emitterSamples = scene.config.integratorSettings.pm.emitterSamplesCount;
        m_directIntegrator->m_bsdfSamples = scene.config.integratorSettings.pm.emitterSamplesCount;
        m_directIntegrator->m_samplingStrategy = ESamplingStrategy::EMIS;

        //1st pass
        m_photonCount = scene.config.integratorSettings.pm.photonCount;
        m_photonRrDepth = scene.config.integratorSettings.pm.photonRrDepth;
        m_photonRrProb = scene.config.integratorSettings.pm.photonRrProb;

        //2nd pass
        m_radiusSearch = scene.config.integratorSettings.pm.searchRadius;
        m_nbPhotonsSearch = scene.config.integratorSettings.pm.photonsSearchCount;
        m_useFinalGather = scene.config.integratorSettings.pm.useFinalGather;
        m_nbFinalGather = scene.config.integratorSettings.pm.finalGatherSamplesCount;
        m_usePhotonsForDirect = scene.config.integratorSettings.pm.usePhotonsForDirect;
    }

    bool init() override {
        Integrator::init();

        std::cout << "Start emitting photons. " << std::endl;
        generatePhotonMap();

        return true;
    }

    void generatePhotonMap() {
        // TODO(A6): Implement this
    }

    v3f render(const Ray& ray, Sampler& sampler) const override {
        v3f throughput(0.f);
        // TODO(A6): Implement this

        return throughput;
    }
};

TR_NAMESPACE_END