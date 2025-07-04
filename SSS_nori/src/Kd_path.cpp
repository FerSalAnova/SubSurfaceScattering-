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

// Simple KD-Tree Node structure
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

    // Build the KD-Tree recursively with a better partitioning strategy
    void build(std::vector<std::pair<Point3f, Color3f>>& points) {
        root = buildRecursive(points, 0);
    }

    // Find the closest point to a given query point
    void nearestNeighbor(const Point3f& query, Point3f& nearestPoint, Color3f& nearestRadiance) const {
        nearestNeighborRecursive(root, query, 0, nearestPoint, nearestRadiance);
    }

private:
    KDTreeNode* root;
    int m_maxPointsPerLeaf;  // Limit the number of points per leaf node

    // Recursive function to build the KD-Tree
    KDTreeNode* buildRecursive(std::vector<std::pair<Point3f, Color3f>>& points, int depth) {
        if (points.size() <= m_maxPointsPerLeaf) {
            // Base case: if the number of points is small enough, make it a leaf node
            Color3f avgRadiance(0.0f);
            for (const auto& p : points) {
                avgRadiance += p.second;
            }
            avgRadiance /= points.size();  // Compute average radiance
            return new KDTreeNode(points[0].first, avgRadiance);
        }

        int axis = depth % 3;
        size_t medianIndex = points.size() / 2;

        // Sort the points by the current axis using the variance-based partitioning strategy
        std::sort(points.begin(), points.end(), [axis](const auto& a, const auto& b) {
            return a.first[axis] < b.first[axis];
        });

        // Create the current node
        Color3f avgRadiance(0.0f);
        for (const auto& p : points) {
            avgRadiance += p.second;
        }
        avgRadiance /= points.size();  // Compute average radiance for this node
        KDTreeNode* node = new KDTreeNode(points[medianIndex].first, avgRadiance);

        // Split the points into left and right subtrees
        std::vector<std::pair<Point3f, Color3f>> leftPoints(points.begin(), points.begin() + medianIndex);
        std::vector<std::pair<Point3f, Color3f>> rightPoints(points.begin() + medianIndex + 1, points.end());

        // Recursively build the left and right subtrees
        node->left = buildRecursive(leftPoints, depth + 1);
        node->right = buildRecursive(rightPoints, depth + 1);

        return node;
    }

    // Recursive function to find the nearest neighbor
    void nearestNeighborRecursive(KDTreeNode* node, const Point3f& query, int depth, Point3f& nearestPoint, Color3f& nearestRadiance) const {
        if (!node) return;

        int axis = depth % 3;
        float dist = (node->point - query).norm();

        // Check if this node is closer and update the nearest point if necessary
        if ((nearestPoint - query).norm() > dist) {
            nearestPoint = node->point;
            nearestRadiance = node->radiance;
        }

        // Determine whether to go left or right based on query's position
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

        // If the query is far from the node, use the average radiance of the node
        if (dist > 0.5f) {  // Some threshold distance where we use the average radiance
            nearestRadiance = node->radiance;
        }
    }
};

class KdPathBSSRDF : public Integrator {
public:
    KdPathBSSRDF(const PropertyList& props) {
        // Constructor code
    }

    void preprocess(const Scene* scene) override {
        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_real_distribution<> dis(0.0, 1.0);

        // Create a map to store KD-Trees for each mesh
        m_kdTrees.clear();  // Clear previous data

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
        Color3f Lo(0.0f);  // Accumulated radiance for this bounce
        Ray3f currentRay = ray;  // Current ray to trace
        Color3f throughput(1.0f);  // Initial throughput (starting with 1.0)
        const int maxDepth = 10;  // Maximum recursion depth
        int depth = 0;  // Start from depth 0

        // Iterative loop for path tracing
        while (depth < maxDepth) {
            Intersection its;
            if (!scene->rayIntersect(currentRay, its))
                return Lo + throughput * scene->getBackground(currentRay);  // Return background if no hit

            // Direct emission contribution (if the intersected object is an emitter)
            if (its.mesh->isEmitter()) {
                const Emitter* emitter = its.mesh->getEmitter();
                EmitterQueryRecord emitterRecord(its.p);  // Initialize emitter query record
                Lo += emitter->eval(emitterRecord);  // Add emitter radiance to Lo
            }

            // Retrieve all emitters in the scene
            const std::vector<Emitter*>& lights = scene->getLights();
            if (lights.empty())
                return Lo;  // No lights, return accumulated radiance

            // Step 2: Emitter sampling
            int lightIdx = std::min(int(sampler->next1D() * lights.size()), (int)lights.size() - 1);
            const Emitter* emitter = lights[lightIdx];
            EmitterQueryRecord emitterRecord(its.p);  // Initialize emitter record at the intersection point
            Color3f Le = emitter->sample(emitterRecord, sampler->next2D(), sampler->next1D());
            emitterRecord.pdf = emitter->pdf(emitterRecord);

            if (!Le.isZero()) {
                // Perform a visibility test between the intersection point and the emitter
                Ray3f shadowRay(its.p, emitterRecord.wi, Epsilon, emitterRecord.dist - Epsilon);
                Intersection lightIts;
                if (!scene->rayIntersect(shadowRay, lightIts)) {
                    // The emitter is visible, compute the BSDF evaluation for this direction
                    BSDFQueryRecord bsdfRecord(its.toLocal(-currentRay.d), its.toLocal(emitterRecord.wi), its.uv, ESolidAngle);
                    Color3f bsdfVal = its.mesh->getBSDF()->eval(bsdfRecord);

                    // Calculate MIS weight for emitter sampling
                    float pdf_emitter = emitterRecord.pdf;
                    float pdf_bsdf = its.mesh->getBSDF()->pdf(bsdfRecord);
                    float w_emitter = pdf_emitter / (pdf_emitter + pdf_bsdf + Epsilon);

                    // Clamp the radiance to avoid fireflies
                    Color3f clampedLe = Le.cwiseMin(10.0f);  // Example threshold for clamping
                    Lo += throughput * w_emitter * clampedLe * bsdfVal / std::max(pdf_emitter, Epsilon);
                }
            }

            // Handle BSSRDF scattering using KD-Tree
            auto bssrdf = its.mesh->getBSDF();
            if (bssrdf && std::string(bssrdf->type()) == "BSSRDF") {
                // Query the KD-Tree for this mesh (access by its mesh pointer or ID)
                Point3f nearestPoint;
                Color3f nearestRadiance;
                auto it = m_kdTrees.find(its.mesh);  // Find the KD-Tree for this mesh
                if (it != m_kdTrees.end()) {
                    it->second.nearestNeighbor(its.p, nearestPoint, nearestRadiance);

                    // Clamp the radiance to avoid fireflies
                    Color3f clampedRadiance = nearestRadiance.cwiseMin(10.0f);
                    float dist = (its.p - nearestPoint).norm();
                    float threshold = 0.5f;
                    if (dist < threshold) {
                        float attenuation = exp(-dist * dist / (2 * threshold * threshold));
                        Lo += clampedRadiance * attenuation * std::max(0.0f, its.shFrame.n.dot(nearestPoint - its.p));
                    }
                }
            }

            // Step 4: BSDF (material) sampling
            BSDFQueryRecord bsdfRecord(its.toLocal(-currentRay.d), its.uv);
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

                        // Clamp the radiance to avoid fireflies
                        Color3f clampedLe = Le_bsdf.cwiseMin(10.0f);  // Example threshold for clamping
                        Lo += throughput * w_bsdf * clampedLe * bsdfSample * std::max(0.0f, its.shFrame.n.dot(wo));
                    }
                }
            }

            // Convert the sampled direction back to world space
            Vector3f wiWorld = its.toWorld(bsdfRecord.wo);

            // Update throughput with the sampled BSDF color
            float cosTheta = std::max(0.0f, Frame::cosTheta(bsdfRecord.wo));
            throughput *= bsdfSample * cosTheta;

            // Russian roulette termination
            if (throughput.maxCoeff() < 0.001f) {  // If throughput is too low
                return Lo;  // Skip adding this path's contribution
            }

            float rrProbability = std::min(throughput.maxCoeff(), 0.90f);  // Adjust the termination probability based on throughput
            if (sampler->next1D() > rrProbability)
                return Lo;  // Terminate path

            throughput /= rrProbability;  // Compensate for the terminated paths

            // Trace the next ray (iteratively)
            currentRay = Ray3f(its.p, wiWorld);
            depth++;  // Increment depth for the next iteration
        }

        return Lo;  // Return accumulated radiance
    }

    std::string toString() const override {
        return "Direct Light with BSSRDF Integrator using KD-Tree []";
    }

private:
    int m_numSamplesPerMesh = 100000;  // Number of points to sample per mesh
    KDTree m_kdTree;  // KD-Tree for storing sampled points and radiances
    std::map<const Mesh*, KDTree> m_kdTrees;
};

NORI_REGISTER_CLASS(KdPathBSSRDF, "Kd_path_bssrdf");

NORI_NAMESPACE_END