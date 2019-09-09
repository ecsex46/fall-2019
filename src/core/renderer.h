/*
    This file is part of TinyRender, an educative rendering system.

    Designed for ECSE 446/546 Realistic/Advanced Image Synthesis.
    Derek Nowrouzezahrai, McGill University.
*/

#pragma once

#include <core/platform.h>
#include <core/core.h>
#include <core/integrator.h>
#include <core/renderpass.h>

TR_NAMESPACE_BEGIN

/**
 * Renderer structure (offline and real-time).
 */
struct Renderer {
    std::unique_ptr<Integrator> integrator;
    std::unique_ptr<RenderPass> renderpass;
    Scene scene;
    bool realTime;
    bool nogui;
    bool realTimeCameraFree;
    unsigned int previousTime = 0, currentTime = 0;
    const int frameDuration = 30;

    explicit Renderer(const Config& config);
    bool init(bool isRealTime, bool nogui);
    void render();
    void cleanUp();
};

TR_NAMESPACE_END