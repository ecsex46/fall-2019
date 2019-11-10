/*
    This file is part of TinyRender, an educative rendering system.

    Designed for ECSE 446/546 Realistic/Advanced Image Synthesis.
    Derek Nowrouzezahrai, McGill University.
*/

#pragma once

#include <GL/glew.h>
#include <functional>
#include <thread>
#include "platform.h"
#include "math.h"
#include "utils.h"
#include "cpptoml.h"
#include "tiny_obj_loader.h"
#include "camera.h"

TR_NAMESPACE_BEGIN

enum ESamplingMeasure{
    ESolidAngleMeasure = 0,
    EAreaMeasure
};

/**
 * Integrator enumeration.
 * A new item needs to be added when creating a new integrator.
 */
enum EIntegrator {
    ENormalIntegrator = 0,
    ESimpleIntegrator,
    EAOIntegrator,
    EROIntegrator,
    EDirectIntegrator,
    EPolygonalIntegrator,
    EPathTracerIntegrator,
    EPhotonMapperIntegrator,
    EIntegrators
};

/**
 * RenderPass enumeration.
 * A new item needs to be added when creating a new renderpass.
 */
enum ERenderPass {
	ENormalRenderPass = 0,
	EDirectRenderPass,
	ESSAORenderPass,
	ERORenderPass,
	EGIRenderPass,
	EPolygonalRenderPass,
    ERenderPasses
};

/**
 * BSDF enumeration.
 */
enum EBSDF {
    EDiffuseBSDF = 0,
    EMirrorBSDF,
    EPhongBSDF,
    EBSDFs
};

/**
 * SamplingType enumeration for Ambient Occlusion integrator
 */

enum ESamplingType { 
	ESpherical = 0, 
	EHemispherical, 
	ECosineHemispherical,
	ESamplingTypes
};

/**
 * Sampling strategy for direct integrator
 */
enum ESamplingStrategy {
    EMIS = 0,
    EArea,
    ESolidAngle,
    ECosineHemisphere,
    EBSDF,
    ESamplingStrategies
};

/**
 * Sampling strategy for polygonal integrator
 */
enum EPolygonalMethod {
    EArvoAnalytic,
    ESurfaceArea,
    EControlVariates,
    EPolygonalMethods
};

// Forward declarations
struct Scene;
struct WorldData;

/**
 * Bounding sphere structure.
 */
struct BSphere {
    inline BSphere() : center(0.0f), radius(0.0f) { }
    inline BSphere(const v3f& center, float radius) : center(center), radius(radius) { }
    inline BSphere(const BSphere& boundingSphere) : center(boundingSphere.center), radius(boundingSphere.radius) { }
    inline bool isEmpty() const { return radius <= 0.0f; }
    inline bool contains(const v3f p) const { return glm::length(p - center) <= radius; }
    v3f center;
    float radius;
};

/**
 * Axis-aligned bounding box.
 */
struct AABB {
    v3f min, max;
    inline AABB() { reset(); }
    inline AABB(const v3f& p) : min(p), max(p) { }
    inline void reset() {
        min = v3f(std::numeric_limits<float>::infinity());
        max = v3f(-std::numeric_limits<float>::infinity());
    }
    inline v3f getCenter() const { return (max + min) * 0.5f; }
    inline BSphere getBSphere() const {
        v3f center = getCenter();
        return BSphere(center, float((center - max).length()));
    }
    inline void expandBy(const v3f& p) {
        for (int i = 0; i < v3f::length(); ++i) {
            min[i] = std::min(min[i], p[i]);
            max[i] = std::max(max[i], p[i]);
        }
    }
    inline void expandBy(const AABB& aabb) {
        for (int i = 0; i < v3f::length(); ++i) {
            min[i] = std::min(min[i], aabb.min[i]);
            max[i] = std::max(max[i], aabb.max[i]);
        }
    }
    inline int getLargestAxis() const {
        v3f d = max - min;
        int largest = 0;

        for (int i=1; i<v3f::length(); ++i)
            if (d[i] > d[largest])
                largest = i;
        return largest;
    }
};

/**
 * Ray structure.
 * Stores origin, direction, min, and max t-values.
 */
struct Ray {
    /* Origin point */
    v3f o;
    /* Ray direction */
    v3f d;
    /* Minimum distance to walk along the ray before allowing intersections */
    float min_t;
    /* Maximum distance to walk along the ray while allowing intersections */
    float max_t;
    Ray(const v3f& co, const v3f& cd, float min_t = Epsilon, float max_t = std::numeric_limits<float>::max())
        : o(co), d(cd), min_t(min_t), max_t(max_t) { }
};

/**
 * Render buffer.
 * Where pixels are stored.
 */
struct RenderBuffer {
    int width, height;
    std::unique_ptr<v3f[]> data;
    RenderBuffer(int w, int h) : width(w), height(h) {
        data = std::unique_ptr<v3f[]>(new v3f[width * height]);
    }
    void scale(float s) {
        for (int i = 0; i < height * width; data[i++] *= s);
    }
    void add(const RenderBuffer& s) {
        for (int i = 0; i < height * width; data[i] += s.data[i], i++);
    }
    void set(const RenderBuffer& s) {
        for (int i = 0; i < height * width; data[i] = s.data[i], i++);
    }
    void clear() {
        for (int i = 0; i < height * width; data[i++] = v3f(0.f));
    }
};

/**
 * Coordinate frame structure.
 * Stores canonical frame and transforms.
 */
struct Frame {
    /* Normal vector of this frame (normal at the surface of the hit) */
    v3f n;
    /* Vector perpendicular to this frame's normal vector and to vector t*/
    v3f s;
    /* Vector perpendicular to this frame's normal vector and to vector t*/
    v3f t;
    explicit Frame() = default;
    explicit Frame(const v3f& n) : n(n) {
        coordinateSystem(n, s, t);
    }
    /* Convert a given vector from World coordinates to this frame's coordinate system */
    v3f toLocal(const v3f& v) const {
        return v3f(glm::dot(v, s), glm::dot(v, t), glm::dot(v, n));
    }
    /* Convert a given vector from this frame's coordinate system to World coordinates */
    v3f toWorld(const v3f& v) const {
        return s * v.x + t * v.y + n * v.z;
    }
    /* Compute the cosine of the angle between the given vector and the xy plane. Assumes the vector has length 1 */
    static float cosTheta(const v3f& v) {
        return v.z;
    }
};

/**
 * Intersection hit structure.
 * Stores hit point incoming/outgoing directions, normal frame, geometry info, etc.
 */
struct SurfaceInteraction {
    /* Hit point position */
    v3f p;
    /* Outgoing direction (towards eye) */
    v3f wo;
    /* Incident direction (towards light source) */
    v3f wi;
    float t, u, v;
    /* ID of the shape that was hit */
    size_t shapeID;
    /* ID of the primitive that was hit */
    size_t primID;
    /* Geometric normal at the hit point on this surface */
    Frame frameNg;
    /* Shading normal at the hit point on this surface */
    Frame frameNs;
    /* ID of the material on this surface */
    int matID;
    unsigned int sampledComponent, sampledType;
};

/**
 * Camera structure.
 * Stores origin, lookat, and up vector.
 */
struct Camera {
    /* Origin point (eye position) */
    v3f o;
    /* Point the camera is looking at (center point) */
    v3f at;
    /* Vector indicating the upwards direction of the camera */
    v3f up;
    /* Full field of view of the camera (vertical) */
    float fov;
};

/**
 * Configuration structure to render a scene.
 * Stores integrator, camera setup, image plane dimensions, sample count, etc.
 */
struct Config {
    /* Which integrator method to use to produce the render */
    EIntegrator integrator;
    /* Which renderpass method to use to produce the render */
    ERenderPass renderpass;
    /* Camera config properties */
    Camera camera;
    fs::path objFile, tomlFile;
    /* width (in pixels) of the image to render */
    int width;
    /* height (in pixels) of the image to render */
    int height;
    /* Number of samples per pixel to take */
    int spp;
	/* Enable code for bonus*/
	bool bonus;
	/* Integer used to specify an automated test */
	bool test;

    struct IntegratorConfig {
        IntegratorConfig() : di{}{};
        ~IntegratorConfig() {}

        /* Config options for the 'direct' integrator */
        struct direct_s {
            /* Number of emitter samples to take */
            size_t emitterSamples{};
            /* Number of bsdf samples to take */
            size_t bsdfSamples{};
            ESamplingStrategy samplingStrategy;
        } di;

        /* Config options for polygonal lights integrator */
        struct poly_s {
            float alpha;
            size_t visSamples;
            bool traceShadows;
            EPolygonalMethod method;
        } poly;

        /* Config options for the photon mapping integrator */
        struct pm_s{
            int photonRrDepth;
            float photonRrProb;

            int photonCount;
            float searchRadius;
            int photonsSearchCount;
            int finalGatherSamplesCount;
            bool useFinalGather;
            size_t emitterSamplesCount;
            bool usePhotonsForDirect;
        } pm{};
        /* Config options for the ambient occlusion integrator */
        struct ao_s{
			ESamplingType sampling_type;
        } ao;
        /* Config options for the reflection occlusion integrator */
        struct ro_s{
            float exponent;
        } ro;
        /* Config options for the path tracer integrator */
        struct pt_s{
            /* Whether to explicitly sample the emitters in the scene */
            bool isExplicit;
            /* Maximum recursion depth */
            int maxDepth;
            /* Depth at which to start Russian Roulette */
            int rrDepth;
            /* Russian Roulette probability (to keep bouncing a ray) */
            float rrProb;
        } pt;
        struct gi_s{
            int maxDepth;
            int rrDepth;
            float rrProb;
            int samplesByVertex;
        } gi;
    } integratorSettings;
};

struct Integrator;

/**
 * Bidirectional scattering distribution function (BSDF) structure.
 * Stores info on scene, methods to evaluate/sample & get PDF, etc.
 */
struct BSDF {

/**
 * BSDF lobe types
 */
    enum EBSDFType {
        // 'null' scattering event, i.e. particles do not undergo deflection
            ENull = 0x00001,
        // Ideally diffuse reflection
            EDiffuseReflection = 0x00002,
        // Ideally diffuse transmission
            EDiffuseTransmission = 0x00004,
        // Glossy reflection
            EGlossyReflection = 0x00008,
        // Glossy transmission
            EGlossyTransmission = 0x00010,
        // Reflection into a discrete set of directions
            EDeltaReflection = 0x00020,
        // Transmission into a discrete set of directions
            EDeltaTransmission = 0x00040,
        // Reflection into a 1D space of directions
            EDelta1DReflection = 0x00080,
        // Transmission into a 1D space of directions
            EDelta1DTransmission = 0x00100
    };

    enum ETypeCombinations {
        // Any reflection component (scattering into discrete, 1D, or 2D set of directions)
            EReflection = EDiffuseReflection | EDeltaReflection | EDelta1DReflection | EGlossyReflection,
        // Any transmission component (scattering into discrete, 1D, or 2D set of directions)
            ETransmission =
            EDiffuseTransmission | EDeltaTransmission | EDelta1DTransmission | EGlossyTransmission | ENull,
        // Diffuse scattering into a 2D set of directions
            EDiffuse = EDiffuseReflection | EDiffuseTransmission,
        // Non-diffuse scattering into a 2D set of directions
            EGlossy = EGlossyReflection | EGlossyTransmission,
        // Scattering into a 2D set of directions
            ESmooth = EDiffuse | EGlossy,
        // Scattering into a discrete set of directions
            EDelta = ENull | EDeltaReflection | EDeltaTransmission,
        // Scattering into a 1D space of directions
            EDelta1D = EDelta1DReflection | EDelta1DTransmission,
        // Any kind of scattering
            EAll = EDiffuse | EGlossy | EDelta | EDelta1D
    };

    const WorldData& worldData;
    const Config& config;
    /* Emission of the material. Has length 0 if not emissive */
    v3f emission;
    std::vector<unsigned int> components;
    unsigned int combinedType;
    BSDF(const WorldData& d, const Config& c, const size_t matID);
    virtual float getExponent(const SurfaceInteraction&) const = 0;
    virtual v3f eval(const SurfaceInteraction&) const = 0;
    virtual float pdf(const SurfaceInteraction&) const = 0;
    virtual v3f sample(SurfaceInteraction&, Sampler&, float* pdf = nullptr) const = 0;
    bool isEmissive() const {
        return glm::length2(emission) > 0.f;
    }
    unsigned int getType() const {
        return combinedType;
    }
    virtual std::string toString() const = 0;
};

/**
 * Emitter/light structure.
 * Stores ID of shape attached to it, and radiance.
 */
struct Emitter {
    /* ID of the shape the emitter is associated with */
    size_t shapeID;
    /* Total surface area of the area light */
    float area;
    /* Emitted radiance (uniform distribution) */
    v3f radiance;
    /* Discrete probability density function over the mesh triangles */
    Distribution1D faceAreaDistribution;
    v3f getRadiance() const { return radiance; }
    v3f getPower() const { return area * M_PI * radiance; }
    bool operator==(const Emitter& other) const { return shapeID == other.shapeID; }
};

/**
 * World data structure.
 * Stores all shapes and BSDFs with their attributes.
 */
struct WorldData {
    tinyobj::attrib_t attrib;
    std::vector<tinyobj::shape_t> shapes;
    std::vector<tinyobj::material_t> materials;
    std::vector<v3f> shapesCenter;
    std::vector<AABB> shapesAABOX;
};

struct AcceleratorBVH;

/**
 * Scene structure.
 * Stores all objects, BVH, list of emitters, list of BSDFs, etc.
 */
struct Scene {
    /* Config properties of the scene and camera, read from the given TOML file */
    const Config& config;
    /* The world data. Contains the vertices, faces and materials of all the objects in the scene */
    WorldData worldData;
    /* Bounding Volume Hierarchy acceleration structure. Provides functions to intersect a ray with the scene geometry. */
    std::unique_ptr<AcceleratorBVH> bvh;
    /* List of all the emitters in the scene */
    std::vector<Emitter> emitters;
    std::vector<std::unique_ptr<BSDF>> bsdfs;
    AABB aabb;

    explicit Scene(const Config& config);
    bool load(bool isRealTime);
    float getShapeArea(size_t shapeID, Distribution1D& faceAreaDistribution);
    float getShapeRadius(const size_t shapeID) const;
    v3f getShapeCenter(const size_t shapeID) const;
    size_t getFirstLight() const;
    v3f getFirstLightPosition() const;
    v3f getFirstLightIntensity() const;

    v3f getObjectVertexPosition(size_t objectIdx, size_t vertexIdx) const;
    v3f getObjectVertexNormal(size_t objectIdx, size_t vertexIdx) const;
    size_t getObjectNbVertices(size_t objectIdx) const;
    int getPrimitiveID(size_t vertexIdx) const;
    int getMaterialID(size_t objectIdx, int primID) const;
};

/**
 * Quick and robust ray-triangle intersection.
 * For a fast method of intersecting a ray with the entire scene, see the BVH object.
 */
inline bool rayTriangleIntersect(const Ray& r,
                                 const v3f& v0,
                                 const v3f& v1,
                                 const v3f& v2,
                                 float& t,
                                 float& u,
                                 float& v) {
    v3f v0v1 = v1 - v0;
    v3f v0v2 = v2 - v0;
    v3f pvec = glm::cross(r.d, v0v2);
    float det = glm::dot(v0v1, pvec);
    if (std::fabs(det) < Epsilon) return false;
    float invDet = 1 / det;
    v3f tvec = r.o - v0;
    u = glm::dot(tvec, pvec) * invDet;
    if (u < 0 || u > 1) return false;
    v3f qvec = glm::cross(tvec, v0v1);
    v = glm::dot(r.d, qvec) * invDet;
    if (v < 0 || u + v > 1) return false;
    t = glm::dot(v0v2, qvec) * invDet;
    return true;
}

/**
 * Texture (templated) structure.
 */
template<class T>
struct Texture {
    virtual T eval(const WorldData&, const SurfaceInteraction&) const = 0;
    virtual T getAverage() const = 0;
    virtual T getMin() const = 0;
    virtual T getMax() const = 0;
};

/**
 * Main texture structure.
 * Stores width, height, and post-processing & loading methods.
 */
struct Tex {
    int w;
    int h;
    std::vector<float> cs;

    void pink() {
        w = 1;
        h = 1;
        cs.emplace_back(1.f);
        cs.emplace_back(0.f);
        cs.emplace_back(1.f);
    }

    // Calculate pixel coordinate of the vertically-flipped image
    inline int fl(int i) {
        const int j = i / 3;
        const int x = j % w;
        const int y = j / w;
        return 3 * ((h - y - 1) * w + x) + i % 3;
    }

    // Post procses a pixel for pmp textures
    inline float pf(int i, float e, std::vector<uint8_t>& ct) {
        // Gamma expansion
        return std::pow(float(ct[fl(i)]) / e, 2.2f);
    }

    // Post procses a pixel for pmf textures
    inline float pf(int i, float e, std::vector<float>& ct) {
        if (e < 0) {
            return ct[fl(i)];
        }
        int32_t m = bswap(*(int32_t*) &ct[fl(i)]);
        float* ptr = reinterpret_cast<float*>(&m);
        return *ptr;
    }

    // Load a ppm or a pfm texture
    template<class T>
    inline void loadpxm(std::vector<float>& c, std::string p) {
        puts(p.c_str());
        static vector<T> ct;
        FILE* f = fopen(p.c_str(), "rb");
        if (!f) {
            std::cout << "Err loading texture : " << p << std::endl;
            pink();
            return;
        }
        int e;
        size_t read = fscanf(f, "%*s %d %d %d%*c", &w, &h, &e);
        if (read == 0) {
            std::cout << "Err loading texture : " << p << std::endl;
            pink();
            return;
        }
        const int sz = w * h * 3;
        ct.assign(sz, 0);
        c.assign(sz, 0);
        read = fread(ct.data(), sizeof(T), sz, f);
        if (read == 0) {
            std::cout << "Err loading texture : " << p << std::endl;
            pink();
            return;
        }
        for (int i = 0; i < sz; i++) {
            c[i] = pf(i, float(e), ct);
        }
        fclose(f);
        std::cout << "Success loading texture : " << p << " w: " << w << " h: " << h << std::endl;
    }

    // Load pfm texture
    inline void loadpfm(string p) {
        loadpxm<float>(cs, p);
    }

    // Load ppm texture
    inline void load(string p) {
        auto b = fs::path(p);
        auto pc = b.replace_extension(".ppm").string();
        auto pa = (b.parent_path() /
            fs::path(b.stem().string() + "_alpha.ppm")).string();
        loadpxm<uint8_t>(cs, pc);
    }
};

struct ConstantTexture3f : Texture<v3f> {
    v3f value;
    explicit ConstantTexture3f(const v3f& v) : value(v) { }
    v3f eval(const WorldData&, const SurfaceInteraction& i) const override {
        return value;
    }
    v3f getAverage() const override { return value; }
    v3f getMin() const override { return value; }
    v3f getMax() const override { return value; }
};

struct ConstantTexture1f : Texture<float> {
    float value;
    explicit ConstantTexture1f(const float& v) : value(v) { }
    float eval(const WorldData& s, const SurfaceInteraction& i) const override {
        return value;
    }
    float getAverage() const override { return value; }
    float getMin() const override { return value; }
    float getMax() const override { return value; }
};

struct BitmapTexture3f : Texture<v3f> {
    std::unique_ptr<Tex> texturePtr;

    explicit BitmapTexture3f(const Config& config, const std::string& filename) {
        texturePtr = std::unique_ptr<Tex>(new Tex());

        fs::path fullpath(config.objFile);
        if (!fullpath.is_absolute())
            fullpath = config.tomlFile.parent_path() / fullpath;

        fs::path file(filename);
        if (!file.is_absolute())
            fullpath = fullpath.parent_path() / file;
        else
            fullpath = file;

        texturePtr->load(fullpath.make_preferred().string());
    }

    v3f getAverage() const override {
        v3f s(0);
        for (size_t i = 0; i < texturePtr->cs.size() / 3; i++) {
            s += v3f(texturePtr->cs[i * 3 + 0], texturePtr->cs[i * 3 + 1], texturePtr->cs[i * 3 + 2]);
        }
        float scale = 1.0 / (texturePtr->cs.size() / 3.0);
        return s * scale;
    }

    v3f getMin() const override {
        v3f s(std::numeric_limits<float>::max(), std::numeric_limits<float>::max(), std::numeric_limits<float>::max());
        for (size_t i = 0; i < texturePtr->cs.size() / 3; i++)
            s = glm::min(s, v3f(texturePtr->cs[i * 3 + 0], texturePtr->cs[i * 3 + 1], texturePtr->cs[i * 3 + 2]));
        return s;
    }

    v3f getMax() const override {
        v3f s(-std::numeric_limits<float>::min(),
              -std::numeric_limits<float>::min(),
              -std::numeric_limits<float>::min());
        for (size_t i = 0; i < texturePtr->cs.size() / 3; i++)
            s = glm::max(s, v3f(texturePtr->cs[i * 3 + 0], texturePtr->cs[i * 3 + 1], texturePtr->cs[i * 3 + 2]));
        return s;
    }

    v3f eval(const WorldData& s, const SurfaceInteraction& hit) const override {
        const tinyobj::shape_t& shape = s.shapes[hit.shapeID];

        int i0 = shape.mesh.indices[hit.primID * 3 + 0].texcoord_index;
        int i1 = shape.mesh.indices[hit.primID * 3 + 1].texcoord_index;
        int i2 = shape.mesh.indices[hit.primID * 3 + 2].texcoord_index;
        v2f st0{s.attrib.texcoords[2 * i0 + 0], s.attrib.texcoords[2 * i0 + 1]};
        v2f st1{s.attrib.texcoords[2 * i1 + 0], s.attrib.texcoords[2 * i1 + 1]};
        v2f st2{s.attrib.texcoords[2 * i2 + 0], s.attrib.texcoords[2 * i2 + 1]};

        v2f st = barycentric(st0, st1, st2, hit.u, hit.v) + v2f(1.0, 1.0);
        st = st - glm::floor(st);

        const int x = clamp(int(st.x * texturePtr->w), 0, texturePtr->w - 1);
        const int y = clamp(int(st.y * texturePtr->h), 0, texturePtr->h - 1);
        const int i = texturePtr->w * y + x;

        return {texturePtr->cs[i * 3 + 0], texturePtr->cs[i * 3 + 1], texturePtr->cs[i * 3 + 2]};
    };
};

struct BitmapTexture1f : Texture<float> {
    std::unique_ptr<Tex> texturePtr;

    explicit BitmapTexture1f(const std::string& filename) {
        texturePtr = std::unique_ptr<Tex>(new Tex());
        texturePtr->load(filename);
    }

    float getAverage() const override {
        float s(0);
        for (size_t i = 0; i < texturePtr->cs.size(); i++) {
            s += texturePtr->cs[i];
        }
        float scale = 1.0f / texturePtr->cs.size();
        return s * scale;
    }

    float getMin() const override {
        float s = std::numeric_limits<float>::max();
        for (size_t i = 0; i < texturePtr->cs.size() / 3; i++)
            s = min(texturePtr->cs[i], s);
        return s;
    }

    float getMax() const override {
        float s = std::numeric_limits<float>::min();
        for (size_t i = 0; i < texturePtr->cs.size() / 3; i++)
            s = max(texturePtr->cs[i], s);
        return s;
    }

    float eval(const WorldData& s, const SurfaceInteraction& hit) const override {
        const tinyobj::shape_t& shape = s.shapes[hit.shapeID];

        int i0 = shape.mesh.indices[hit.primID * 3 + 0].texcoord_index;
        int i1 = shape.mesh.indices[hit.primID * 3 + 1].texcoord_index;
        int i2 = shape.mesh.indices[hit.primID * 3 + 2].texcoord_index;
        v2f st0{s.attrib.texcoords[2 * i0 + 0], s.attrib.texcoords[2 * i0 + 1]};
        v2f st1{s.attrib.texcoords[2 * i1 + 0], s.attrib.texcoords[2 * i1 + 1]};
        v2f st2{s.attrib.texcoords[2 * i2 + 0], s.attrib.texcoords[2 * i2 + 1]};

        v2f st = barycentric(st0, st1, st2, hit.u, hit.v) + v2f(1.0, 1.0);
        st = st - glm::floor(st);

        const int x = clamp(int(st.x * texturePtr->w), 0, texturePtr->w - 1);
        const int y = clamp(int(st.y * texturePtr->h), 0, texturePtr->h - 1);
        const int i = texturePtr->w * y + x;

        return texturePtr->cs[i];
    };
};

/**
 * Utilities for managing multi-threading
 * Provides static methods for running a `for` loop in parallel
 */
struct ThreadPool {

    /* Run a `for` loop in parallel */
	template<typename Index, typename Callable>
	static void ParallelFor(Index start, Index end, Callable func) {
		// Estimate number of threads in the pool
		const static unsigned nb_threads_hint = std::thread::hardware_concurrency();
		const static unsigned nb_threads = (nb_threads_hint == 0u ? 8u : nb_threads_hint);

		// Size of a slice for the range functions
		Index n = end - start + 1;
		Index slice = (Index)std::round(n / static_cast<double> (nb_threads));
		slice = std::max(slice, Index(1));

		// [Helper] Inner loop
		auto launchRange = [&func](int k1, int k2) {
			for (Index k = k1; k < k2; k++) {
				func(k);
			}
		};

		// Create pool and launch jobs
		std::vector<std::thread> pool;
		pool.reserve(nb_threads);
		Index i1 = start;
		Index i2 = std::min(start + slice, end);
		for (unsigned i = 0; i + 1 < nb_threads && i1 < end; ++i) {
			pool.emplace_back(launchRange, i1, i2);
			i1 = i2;
			i2 = std::min(i2 + slice, end);
		}
		if (i1 < end) {
			pool.emplace_back(launchRange, i1, end);
		}

		// Wait for jobs to finish
		for (std::thread &t : pool) {
			if (t.joinable()) {
				t.join();
			}
		}
	}

    /* Run a `for` loop sequentially for easy comparison. Equivalent to a regular C++ `for` loop */
	template<typename Index, typename Callable>
	static void SequentialFor(Index start, Index end, Callable func) {
		for (Index i = start; i < end; i++) {
			func(i);
		}
	}
};

/**
 * Converts 3-dimensional vector to string.
 * 2 digits after decimal precision.
 */
inline std::string toString(const v3f v) {
    std::ostringstream s;
    s << std::fixed << std::setprecision(2);
    std::vector<float> vec{v.x, v.y, v.z};
    std::copy(vec.begin(), vec.end() - 1, std::ostream_iterator<float>(s, ", "));
    s << vec.back();
    return "[" + s.str() + "]";
}

typedef std::function<v3f(const Ray&, Sampler&)> renderer_t;

TR_NAMESPACE_END
