#include <nori/warp.h>
#include <nori/integrator.h>
#include <nori/scene.h>
#include <nori/emitter.h>
#include <nori/bsdf.h>
#include <iostream>
#include <vector>
#include <cstdlib>
#include <random>
#include <map>
#include <cmath>

NORI_NAMESPACE_BEGIN

// Optimized KD-Tree Node structure with better splitting strategy
struct KDTreeNode {
    Point3f point;           // The point for this node
    Color3f radiance;        // The average radiance of the points in this node
    KDTreeNode* left;        // Left child
    KDTreeNode* right;       // Right child

    // Constructor that calculates average radiance for this node
    KDTreeNode(const Point3f& p, const Color3f& r) : point(p), radiance(r), left(nullptr), right(nullptr) {}
};

class KDTree {
public:
    KDTree(int maxPointsPerLeaf = 10) : root(nullptr), m_maxPointsPerLeaf(maxPointsPerLeaf) {}

    void build(std::vector<std::pair<Point3f, Color3f>>& points) {
        root = buildRecursive(points, 0);
    }

    void nearestNeighbor(const Point3f& query, Point3f& nearestPoint, Color3f& nearestRadiance) const {
        nearestNeighborRecursive(root, query, 0, nearestPoint, nearestRadiance);
    }

private:
    KDTreeNode* root;
    int m_maxPointsPerLeaf;

    KDTreeNode* buildRecursive(std::vector<std::pair<Point3f, Color3f>>& points, int depth) {
        if (points.size() <= m_maxPointsPerLeaf) {
            Color3f avgRadiance(0.0f);
            for (const auto& p : points) {
                avgRadiance += p.second;
            }
            avgRadiance /= points.size();
            return new KDTreeNode(points[0].first, avgRadiance);
        }

        int axis = depth % 3;
        size_t medianIndex = points.size() / 2;

        // Sort by the current axis
        std::sort(points.begin(), points.end(), [axis](const auto& a, const auto& b) {
            return a.first[axis] < b.first[axis];
        });

        Color3f avgRadiance(0.0f);
        for (const auto& p : points) {
            avgRadiance += p.second;
        }
        avgRadiance /= points.size();
        
        KDTreeNode* node = new KDTreeNode(points[medianIndex].first, avgRadiance);

        // Split the points
        std::vector<std::pair<Point3f, Color3f>> leftPoints(points.begin(), points.begin() + medianIndex);
        std::vector<std::pair<Point3f, Color3f>> rightPoints(points.begin() + medianIndex + 1, points.end());

        node->left = buildRecursive(leftPoints, depth + 1);
        node->right = buildRecursive(rightPoints, depth + 1);

        return node;
    }

    void nearestNeighborRecursive(KDTreeNode* node, const Point3f& query, int depth, Point3f& nearestPoint, Color3f& nearestRadiance) const {
        if (!node) return;

        int axis = depth % 3;
        float dist = (node->point - query).norm();

        if ((nearestPoint - query).norm() > dist) {
            nearestPoint = node->point;
            nearestRadiance = node->radiance;
        }

        if (query[axis] < node->point[axis]) {
            nearestNeighborRecursive(node->left, query, depth + 1, nearestPoint, nearestRadiance);
            if (std::fabs(query[axis] - node->point[axis]) < (nearestPoint - query).norm()) {
                nearestNeighborRecursive(node->right, query, depth + 1, nearestPoint, nearestRadiance);
            }
        } else {
            nearestNeighborRecursive(node->right, query, depth + 1, nearestPoint, nearestRadiance);
            if (std::fabs(query[axis] - node->point[axis]) < (nearestPoint - query).norm()) {
                nearestNeighborRecursive(node->left, query, depth + 1, nearestPoint, nearestRadiance);
            }
        }
    }
};

class KdBSSRDF : public Integrator {
public:
    KdBSSRDF(const PropertyList &props) {
        // Constructor code
    }

    void preprocess(const Scene* scene) override {
        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_real_distribution<> dis(0.0, 1.0);

        // Create a map to store KD-Trees for each mesh
        m_kdTrees.clear(); // Clear previous data

        // Iterate through all meshes in the scene
        for (const auto& mesh : scene->getMeshes()) {
            auto bssrdf = mesh->getBSDF();
            if (bssrdf && std::string(bssrdf->type()) == "BSSRDF") {
                std::vector<std::pair<Point3f, Color3f>> allSampledPoints;

                // Sample points from this mesh
                for (int i = 0; i < m_numSamplesPerMesh; ++i) {
                    Point2f sample = Point2f(dis(gen), dis(gen));
                    Point3f p;
                    Normal3f n;
                    Point2f uv;
                    mesh->samplePosition(sample, p, n, uv);

                    Vector3f wi = Warp::squareToCosineHemisphere(Point2f(dis(gen), dis(gen)));

                    BSDFQueryRecord bRec(wi, Vector3f(0.f, 0.f, 0.f), uv, ESolidAngle);
                    Color3f radiance = bssrdf->eval(bRec);

                    // Store points and radiances for the KD-Tree
                    allSampledPoints.emplace_back(p, radiance);
                }

                // Build the KD-Tree for this specific mesh
                KDTree kdTree;
                kdTree.build(allSampledPoints);

                // Store the KD-Tree in the map, using mesh pointer or a unique identifier as the key
                m_kdTrees[mesh] = kdTree;
            }
        }

        std::cout << "Preprocessed and built KD-Trees for " << m_kdTrees.size() << " BSSRDF meshes." << std::endl;
    }

    Color3f Li(const Scene* scene, Sampler* sampler, const Ray3f& ray) const override {
        Color3f Lo(0.0f);
        Intersection its;
        
        if (!scene->rayIntersect(ray, its)) {
            return scene->getBackground(ray);
        }

        // Lighting for emitters
        if (its.mesh->isEmitter()) {
            const Emitter* emitter = its.mesh->getEmitter();
            EmitterQueryRecord emitterRecord(emitter, ray.o, its.p, its.shFrame.n, its.uv);
            Lo += emitter->eval(emitterRecord);
        }

        const std::vector<Emitter*>& lights = scene->getLights();
        if (lights.empty()) return Lo;

        // Sample a light
        int lightIdx = std::min(int(sampler->next1D() * lights.size()), (int)lights.size() - 1);
        const Emitter* emitter = lights[lightIdx];
        EmitterQueryRecord emitterRecord(its.p);
        Color3f Le = emitter->sample(emitterRecord, sampler->next2D(), sampler->next1D());
        emitterRecord.pdf = emitter->pdf(emitterRecord);

        if (!Le.isZero()) {
            Ray3f shadowRay(its.p, emitterRecord.wi, Epsilon, emitterRecord.dist - Epsilon);
            if (!scene->rayIntersect(shadowRay)) {
                BSDFQueryRecord bsdfRecord(its.toLocal(-ray.d), its.toLocal(emitterRecord.wi), its.uv, ESolidAngle);
                Color3f bsdfVal = its.mesh->getBSDF()->eval(bsdfRecord);

                float pdf_emitter = emitterRecord.pdf;
                float pdf_bsdf = its.mesh->getBSDF()->pdf(bsdfRecord);
                float w_emitter = pdf_emitter / (pdf_emitter + pdf_bsdf + Epsilon);

                float cosTheta = std::max(0.0f, its.shFrame.n.dot(emitterRecord.wi));
                Lo += w_emitter * Le * bsdfVal * cosTheta / std::max(pdf_emitter, Epsilon);
            }
        }

        // Handle subsurface scattering with BSSRDF
        auto bssrdf = its.mesh->getBSDF();
        if (bssrdf && std::string(bssrdf->type()) == "BSSRDF") {
            // Query the KD-Tree for this mesh (access by its mesh pointer or ID)
            Point3f nearestPoint;
            Color3f nearestRadiance;
            auto it = m_kdTrees.find(its.mesh);  // Find the KD-Tree for this mesh
            float cosTheta = std::max(0.0f, its.shFrame.n.dot(emitterRecord.wi));
            if (it != m_kdTrees.end()) {
                it->second.nearestNeighbor(its.p, nearestPoint, nearestRadiance);
                float dist = (its.p - nearestPoint).norm();
                float threshold = 0.5f;
                if (dist < threshold) {
                    float attenuation = exp(-dist * dist / (2 * threshold * threshold));
                    Lo += nearestRadiance * attenuation * cosTheta;
                }
            }
        }

        // BSDF evaluation for direct illumination
        BSDFQueryRecord bsdfRecord(its.toLocal(-ray.d), its.uv);
        Color3f bsdfSample = its.mesh->getBSDF()->sample(bsdfRecord, sampler->next2D());

        if (!bsdfSample.isZero()) {
            Vector3f wo = its.toWorld(bsdfRecord.wo);
            float pdf_bsdf = its.mesh->getBSDF()->pdf(bsdfRecord);

            Ray3f shadowRay(its.p, wo, Epsilon, emitterRecord.dist - Epsilon);
            Intersection lightIts;
            if (scene->rayIntersect(shadowRay, lightIts)) {
                if (lightIts.mesh->isEmitter()) {
                    const Emitter* hitEmitter = lightIts.mesh->getEmitter();
                    EmitterQueryRecord lightRecord(hitEmitter, its.p, lightIts.p, lightIts.shFrame.n, its.uv);
                    Color3f Le_bsdf = hitEmitter->eval(lightRecord);

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
        return "Direct Light with BSSRDF Integrator using Optimized KD-Tree []";
    }

private:
    int m_numSamplesPerMesh = 100000;  // Number of points to sample per mesh
    KDTree m_kdTree;  // KD-Tree for storing sampled points and radiances
    std::map<const Mesh*, KDTree> m_kdTrees;
};

NORI_REGISTER_CLASS(KdBSSRDF, "Kd_bssrdf");

NORI_NAMESPACE_END

