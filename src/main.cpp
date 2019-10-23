/*
    This file is part of TinyRender, an educative rendering system.

    Designed for ECSE 446/546 Realistic/Advanced Image Synthesis.
    Derek Nowrouzezahrai, McGill University.
*/

#include <core/core.h>
#include <core/platform.h>
#include <core/renderer.h>
#define TINYEXR_IMPLEMENTATION
#include "tinyexr.h"
#define TINYOBJLOADER_IMPLEMENTATION
#include "tiny_obj_loader.h"


/**
 * Load TOML scene file and create scene objects.
 */
bool loadTOML(TinyRender::Config& config, const std::string& inputFile) {
    // Scene and Wavefront OBJ files
    const auto data = cpptoml::parse_file(inputFile);
    config.tomlFile = inputFile;
    const auto input = data->get_table("input");
    config.objFile = *input->get_as<std::string>("objfile");

    // Camera settings
    const auto camera = data->get_table("camera");
    config.camera.fov = camera->get_as<double>("fov").value_or(30.);
    auto eye = camera->get_array_of<double>("eye").value_or({1., 1., 0.});
    config.camera.o = v3f(eye[0], eye[1], eye[2]);
    auto at = camera->get_array_of<double>("at").value_or({0., 0., 0.});
    config.camera.at = v3f(at[0], at[1], at[2]);
    auto up = camera->get_array_of<double>("up").value_or({0., 1., 0.});
    config.camera.up = v3f(up[0], up[1], up[2]);

    // Film settings
    const auto film = data->get_table("film");
    config.width = film->get_as<int>("width").value_or(768);
    config.height = film->get_as<int>("height").value_or(576);

    // Renderer settings
    const auto renderer = data->get_table("renderer");
    auto realTime = renderer->get_as<bool>("realtime").value_or(false);
    auto type = renderer->get_as<std::string>("type").value_or("normal");

	// Bonus 
	auto bonus = renderer->get_as<bool>("bonus").value_or(false);
	config.bonus = bonus;

	// Test for automated scripting
	auto test = renderer->get_as<bool>("test").value_or(false);
	config.test = test;
		
    // Real-time renderpass
    if (realTime) {
        if (type == "normal") {
            config.renderpass = TinyRender::ENormalRenderPass;
        }
        else if (type == "direct") {
            config.renderpass = TinyRender::EDirectRenderPass;
        }
        else if (type == "ssao") {
            config.renderpass = TinyRender::ESSAORenderPass;
        }		

		else if (type == "polygonal") {
			config.renderpass = TinyRender::EPolygonalRenderPass;
            config.integratorSettings.poly.alpha = renderer->get_as<double>("alpha").value_or(1);
            config.integratorSettings.poly.visSamples = renderer->get_as<size_t>("visSamples").value_or(1);
		}

        else if (type == "gi") {
            config.renderpass = TinyRender::EGIRenderPass;
            config.integratorSettings.gi.maxDepth = renderer->get_as<int>("maxDepth").value_or(5);
            config.integratorSettings.gi.rrDepth = renderer->get_as<int>("rrDepth").value_or(5);
            config.integratorSettings.gi.rrProb = renderer->get_as<double>("rrProb").value_or(9.95f);
            config.integratorSettings.gi.samplesByVertex = renderer->get_as<int>("samplesByVertex").value_or(100);
        }
        else {
            throw std::runtime_error("Invalid renderpass type");
        }
    }

    // Offline integrator
    else {
        if (type == "normal") {
            config.integrator = TinyRender::ENormalIntegrator;
        }
        else if (type == "simple") {
            config.integrator = TinyRender::ESimpleIntegrator;
        }
        else if (type == "ao") {
            config.integrator = TinyRender::EAOIntegrator;
			auto sample_type_str = renderer->get_as<std::string>("sampling_type").value_or("cosine");
			if (sample_type_str == "sphere")
				config.integratorSettings.ao.sampling_type = TinyRender::ESamplingType::ESpherical;
			else if (sample_type_str == "hemisphere")
				config.integratorSettings.ao.sampling_type = TinyRender::ESamplingType::EHemispherical;
			else
				config.integratorSettings.ao.sampling_type = TinyRender::ESamplingType::ECosineHemispherical;

        }
        else if (type == "ro") {
            config.integrator = TinyRender::EROIntegrator;
            config.integratorSettings.ro.exponent = renderer->get_as<double>("exponent").value_or(30);
        }
        else if (type == "direct") {
            config.integrator = TinyRender::EDirectIntegrator;
            config.integratorSettings.di.emitterSamples = renderer->get_as<size_t>("emitterSamples").value_or(1);
            config.integratorSettings.di.bsdfSamples = renderer->get_as<size_t>("bsdfSamples").value_or(1);
            string samplingStrategy = renderer->get_as<string>("samplingStrategy").value_or("emitter");
            if (samplingStrategy == "mis")
                config.integratorSettings.di.samplingStrategy = TinyRender::ESamplingStrategy::EMIS;
            else if (samplingStrategy == "area")
                config.integratorSettings.di.samplingStrategy = TinyRender::ESamplingStrategy::EArea;
            else if (samplingStrategy == "solidAngle")
                config.integratorSettings.di.samplingStrategy = TinyRender::ESamplingStrategy::ESolidAngle;
            else if (samplingStrategy == "cosineHemisphere")
                config.integratorSettings.di.samplingStrategy = TinyRender::ESamplingStrategy::ECosineHemisphere;
            else
                config.integratorSettings.di.samplingStrategy = TinyRender::ESamplingStrategy::EBSDF;
        }
        else if (type == "polygonal") {
            config.integrator = TinyRender::EPolygonalIntegrator;
            config.integratorSettings.poly.alpha = renderer->get_as<double>("alpha").value_or(1);
            config.integratorSettings.poly.visSamples = renderer->get_as<size_t>("visSamples").value_or(1);
            config.integratorSettings.poly.traceShadows = renderer->get_as<bool>("traceShadows").value_or(false);
            string method = renderer->get_as<string>("method").value_or("area");
            if (method == "area")
                config.integratorSettings.poly.method = TinyRender::EPolygonalMethod::ESurfaceArea;
            else if (method == "controlVariates")
                config.integratorSettings.poly.method = TinyRender::EPolygonalMethod::EControlVariates;
            else
                config.integratorSettings.poly.method = TinyRender::EPolygonalMethod::EArvoAnalytic;
        }
        else if (type == "path") {
            config.integrator = TinyRender::EPathTracerIntegrator;
            config.integratorSettings.pt.isExplicit = renderer->get_as<bool>("isExplicit").value_or(true);
            config.integratorSettings.pt.maxDepth = renderer->get_as<int>("maxDepth").value_or(-1);
            config.integratorSettings.pt.rrDepth = renderer->get_as<int>("rrDepth").value_or(5);
            config.integratorSettings.pt.rrProb = renderer->get_as<double>("rrProb").value_or(0.95f);
        }
        else {
            throw std::runtime_error("Invalid integrator type");
        }

        config.spp = renderer->get_as<int>("spp").value_or(1);
    }

    return realTime;
}

/**
 * Launch rendering job.
 */
void run(std::string& inputTOMLFile, bool nogui) {
    TinyRender::Config config;
    bool isRealTime;

    try {
        isRealTime = loadTOML(config, inputTOMLFile);
    } catch (std::exception const& e) {
        std::cerr << "Error while parsing scene file: " << e.what() << std::endl;
        exit(EXIT_FAILURE);
    }

    TinyRender::Renderer renderer(config);
    renderer.init(isRealTime, nogui);
    renderer.render();
    renderer.cleanUp();
}

/**
 * Main TinyRender program.
 */
int main(int argc, char* argv[]) {
    if (argc != 2 && argc !=3) {
        cerr << "Syntax: " << argv[0] << " <scene.toml>" << endl;
        exit(EXIT_FAILURE);
    }

    bool nogui = false;
    if(argc == 3) {
        if(std::string(argv[2]) == "nogui") {
            nogui = true;
        }
    }

    auto inputTOMLFile = std::string(argv[1]);
    run(inputTOMLFile, nogui);

#ifdef _WIN32
    if(!nogui) system("pause");
#endif

    return EXIT_SUCCESS;
}
