#include <nori/warp.h>
#include <nori/integrator.h>
#include <nori/scene.h>
#include <nori/emitter.h>
#include <nori/bsdf.h>
#include <iostream>
#include <vector>
#include <cstdlib>
#include <random>

NORI_NAMESPACE_BEGIN

class DirectBSSRDF : public Integrator {
public:
    DirectBSSRDF(const PropertyList &props) {
        // Constructor code, no parameters for now
    }
    std::map<std::string, std::vector<Point3f>> m_sampledPointsPerMesh;
    std::map<std::string, std::vector<Color3f>> m_radianceAtPointsPerMesh;
    void preprocess(const Scene* scene) override {
        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_real_distribution<> dis(0.0, 1.0);
        
        // Iterate through all meshes in the scene
        for (const auto& mesh : scene->getMeshes()) {
            // Check if the mesh has a BSSRDF
            auto bssrdf = mesh->getBSDF();
            if (bssrdf && std::string(bssrdf->type()) == "BSSRDF") {
                std::vector<Point3f> sampledPoints;
                std::vector<Color3f> radiances;

                for (int i = 0; i < m_numSamplesPerMesh; ++i) {
                    Point2f sample = Point2f(dis(gen), dis(gen));
                    Point3f p;
                    Normal3f n;
                    Point2f uv;
                    mesh->samplePosition(sample, p, n, uv);  // Sample a point on the mesh

                    // Create a random incoming direction for the evaluation
                    Vector3f wi = Warp::squareToCosineHemisphere(Point2f(dis(gen), dis(gen)));

                    // Construct a BSDFQueryRecord
                    BSDFQueryRecord bRec(wi, Vector3f(0.f, 0.f, 0.f), uv, ESolidAngle);
                    Color3f radiance = bssrdf->eval(bRec);
                
                    sampledPoints.push_back(p);
                    radiances.push_back(radiance);
                }

                // Store the sampled points and radiances for this mesh
                m_sampledPointsPerMesh[mesh->getName()] = sampledPoints;
                m_radianceAtPointsPerMesh[mesh->getName()] = radiances;

                std::cout << "Sampled " << sampledPoints.size() << " points on BSSRDF mesh " << mesh->getName() << "." << std::endl;
            }
        }
    }

    Color3f Li(const Scene* scene, Sampler* sampler, const Ray3f& ray) const override {
        Color3f Lo(0.0f);

        // Find the surface that is visible in the requested direction
        Intersection its;
        if (!scene->rayIntersect(ray, its))
            return scene->getBackground(ray); // Return background if no hit

        // Step 1: Direct emission if intersection is on an emitter
        if (its.mesh->isEmitter()) {
            const Emitter* emitter = its.mesh->getEmitter();
            EmitterQueryRecord emitterRecord(emitter, ray.o, its.p, its.shFrame.n, its.uv);
            Lo += emitter->eval(emitterRecord);
        }

        // Retrieve all emitters in the scene
        const std::vector<Emitter*>& lights = scene->getLights();
        if (lights.empty())
            return Lo;

        // Step 2: Emitter sampling
        int lightIdx = std::min(int(sampler->next1D() * lights.size()), (int)lights.size() - 1);
        const Emitter* emitter = lights[lightIdx];
        EmitterQueryRecord emitterRecord(its.p);
        Color3f Le = emitter->sample(emitterRecord, sampler->next2D(), sampler->next1D());
        emitterRecord.pdf = emitter->pdf(emitterRecord);

        if (!Le.isZero()) {
            // Perform a visibility test
            Ray3f shadowRay(its.p, emitterRecord.wi, Epsilon, emitterRecord.dist - Epsilon);
            if (!scene->rayIntersect(shadowRay)) {
                // Evaluate BSDF for emitter sampling
                BSDFQueryRecord bsdfRecord(its.toLocal(-ray.d), its.toLocal(emitterRecord.wi), its.uv, ESolidAngle);
                Color3f bsdfVal = its.mesh->getBSDF()->eval(bsdfRecord);

                // Calculate MIS weight for emitter sampling
                float pdf_emitter = emitterRecord.pdf;
                float pdf_bsdf = its.mesh->getBSDF()->pdf(bsdfRecord);
                float w_emitter = pdf_emitter / (pdf_emitter + pdf_bsdf + Epsilon);

                float cosTheta = std::max(0.0f, its.shFrame.n.dot(emitterRecord.wi));
                Lo += w_emitter * Le * bsdfVal * cosTheta / std::max(pdf_emitter, Epsilon);
            }
        }

        // Step 3: BSSRDF (Subsurface scattering) sampling
        if (its.mesh->getBSDF()->type() == "BSSRDF") {
            BSDFQueryRecord bsdfRecord(its.p, its.shFrame.n, its.uv, ESolidAngle);
            Color3f bssrdfSample = its.mesh->getBSDF()->sample(bsdfRecord, sampler->next2D());

            if (!bssrdfSample.isZero()) {
                Vector3f wo = its.toWorld(bsdfRecord.wo);
                float pdf_bssrdf = its.mesh->getBSDF()->pdf(bsdfRecord);

                // Setup a new shadow ray to test visibility
                Ray3f shadowRay(its.p, wo, Epsilon, emitterRecord.dist - Epsilon);
                Intersection lightIts;
                if (scene->rayIntersect(shadowRay, lightIts)) {
                    if (lightIts.mesh->isEmitter()) {
                        const Emitter* hitEmitter = lightIts.mesh->getEmitter();
                        EmitterQueryRecord lightRecord(hitEmitter, its.p, lightIts.p, lightIts.shFrame.n, its.uv);
                        Color3f Le_bssrdf = hitEmitter->eval(lightRecord);

                        // Calculate MIS weight for BSSRDF sampling
                        float pdf_emitter = hitEmitter->pdf(lightRecord);
                        float w_bssrdf = pdf_bssrdf / (pdf_emitter + pdf_bssrdf + Epsilon);

                        float cosTheta = std::max(0.0f, its.shFrame.n.dot(wo));
                        Lo += w_bssrdf * Le_bssrdf * bssrdfSample;
                    }
                }
            }

            // Gradually incorporate the precomputed radiance for subsurface scattering
            for (size_t i = 0; i < m_sampledPointsPerMesh.at(its.mesh->getName()).size(); ++i) {
                const Point3f& samplePoint = m_sampledPointsPerMesh.at(its.mesh->getName())[i];
                const Color3f& sampleRadiance = m_radianceAtPointsPerMesh.at(its.mesh->getName())[i];

                float dist = (its.p - samplePoint).norm();
                float threshold = 1.25;
                if (dist < threshold) {
                    float attenuation = exp(-dist * dist / (2 * threshold * threshold))/3500;  // Gaussian falloff
                    Lo += sampleRadiance * attenuation;
                }
            }
        }

        // Step 4: BSDF (material) sampling
        BSDFQueryRecord bsdfRecord(its.toLocal(-ray.d), its.uv);
        Color3f bsdfSample = its.mesh->getBSDF()->sample(bsdfRecord, sampler->next2D());

        if (!bsdfSample.isZero()) {
            Vector3f wo = its.toWorld(bsdfRecord.wo);
            float pdf_bsdf = its.mesh->getBSDF()->pdf(bsdfRecord);

            // Setup a new shadow ray to test visibility
            Ray3f shadowRay(its.p, wo, Epsilon, emitterRecord.dist - Epsilon);
            Intersection lightIts;
            if (scene->rayIntersect(shadowRay, lightIts)) {
                if (lightIts.mesh->isEmitter()) {
                    const Emitter* hitEmitter = lightIts.mesh->getEmitter();
                    EmitterQueryRecord lightRecord(hitEmitter, its.p, lightIts.p, lightIts.shFrame.n, its.uv);
                    Color3f Le_bsdf = hitEmitter->eval(lightRecord);

                    // Calculate MIS weight for BSDF sampling
                    float pdf_emitter = hitEmitter->pdf(lightRecord);
                    float w_bsdf = pdf_bsdf / (pdf_emitter + pdf_bsdf + Epsilon);

                    float cosTheta = std::max(0.0f, its.shFrame.n.dot(wo));
                    Lo += w_bsdf * Le_bsdf * bsdfSample;
                }
            }
        }

        return Lo;
    }

    std::string toString() const override {
        return "Direct Light with BSSRDF Integrator []";
    }

private:
    int m_numSamplesPerMesh = 800;  // Number of points to sample per mesh
    std::vector<Point3f> m_sampledPoints;  // Store sampled points on BSSRDF meshes
    std::vector<Color3f> m_radianceAtPoints;  // Store radiance at each sampled point
};

// Register the integrator in the Nori framework
NORI_REGISTER_CLASS(DirectBSSRDF, "direct_bssrdf");

NORI_NAMESPACE_END