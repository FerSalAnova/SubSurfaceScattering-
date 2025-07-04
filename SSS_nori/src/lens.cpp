#include <nori/camera.h>
#include <nori/rfilter.h>
#include <nori/warp.h>
#include <Eigen/Geometry>

NORI_NAMESPACE_BEGIN

class ThinLensCamera : public Camera {
public:
    ThinLensCamera(const PropertyList &propList) {
        /* Width and height in pixels. Default: 720p */
        m_outputSize.x() = propList.getInteger("width", 1280);
        m_outputSize.y() = propList.getInteger("height", 720);
        m_invOutputSize = m_outputSize.cast<float>().cwiseInverse();

        /* Specifies an optional camera-to-world transformation. Default: none */
        m_cameraToWorld = propList.getTransform("toWorld", Transform());

        /* Horizontal field of view in degrees */
        m_fov = propList.getFloat("fov", 30.0f);

        /* Near and far clipping planes in world-space units */
        m_nearClip = propList.getFloat("nearClip", 1e-4f);
        m_farClip = propList.getFloat("farClip", 1e4f);

        /* Depth of field parameters */
        m_lensRadius = propList.getFloat("lensRadius", 0.0f);
        m_focalDistance = propList.getFloat("focalDistance", 10.0f);

        m_rfilter = NULL;
    }

    void activate() {
        float aspect = m_outputSize.x() / (float) m_outputSize.y();

        /* Project vectors in camera space onto a plane at z=1 */
        float recip = 1.0f / (m_farClip - m_nearClip),
              cot = 1.0f / std::tan(degToRad(m_fov / 2.0f));

        Eigen::Matrix4f perspective;
        perspective <<
            cot, 0,   0,   0,
            0, cot,   0,   0,
            0,   0,   m_farClip * recip, -m_nearClip * m_farClip * recip,
            0,   0,   1,   0;

        m_sampleToCamera = Transform(
            Eigen::DiagonalMatrix<float, 3>(Vector3f(0.5f, -0.5f * aspect, 1.0f)) *
            Eigen::Translation<float, 3>(1.0f, -1.0f/aspect, 0.0f) * perspective).inverse();

        if (!m_rfilter) {
            m_rfilter = static_cast<ReconstructionFilter *>(
                NoriObjectFactory::createInstance("gaussian", PropertyList()));
            m_rfilter->activate();
        }
    }

    Color3f sampleRay(Ray3f &ray, const Point2f &samplePosition, const Point2f &apertureSample) const {
        /* Compute the corresponding position on the near plane (in local camera space) */
        Point3f nearP = m_sampleToCamera * Point3f(
            samplePosition.x() * m_invOutputSize.x(),
            samplePosition.y() * m_invOutputSize.y(), 0.0f);

        Vector3f d = nearP.normalized();
        float invZ = 1.0f / d.z();
        ray.o = m_cameraToWorld * Point3f(0, 0, 0);
        ray.d = m_cameraToWorld * d;

        if (m_lensRadius > 0) {
            /* Compute point on plane of focus */
            Point2f pLens = m_lensRadius * Warp::squareToUniformDisk(apertureSample);
            float tf = m_focalDistance / d.z();
            Point3f pFocus = d * tf;

            /* Update ray for effect of lens */
            Point3f o(pLens.x(), pLens.y(), 0);
            ray.o = m_cameraToWorld * o;
            d = (pFocus - o).normalized();
            ray.d = m_cameraToWorld * d;
        }

        ray.mint = m_nearClip * invZ;
        ray.maxt = m_farClip * invZ;
        ray.update();

        return Color3f(1.0f);
    }

    void addChild(NoriObject *obj) {
        switch (obj->getClassType()) {
            case EReconstructionFilter:
                if (m_rfilter)
                    throw NoriException("Camera: tried to register multiple reconstruction filters!");
                m_rfilter = static_cast<ReconstructionFilter *>(obj);
                break;

            default:
                throw NoriException("Camera::addChild(<%s>) is not supported!",
                    classTypeName(obj->getClassType()));
        }
    }

    std::string toString() const {
        return tfm::format(
            "ThinLensCamera[\n"
            "  cameraToWorld = %s,\n"
            "  outputSize = %s,\n"
            "  fov = %f,\n"
            "  clip = [%f, %f],\n"
            "  lensRadius = %f,\n"
            "  focalDistance = %f,\n"
            "  rfilter = %s\n"
            "]",
            indent(m_cameraToWorld.toString(), 18),
            m_outputSize.toString(),
            m_fov,
            m_nearClip,
            m_farClip,
            m_lensRadius,
            m_focalDistance,
            indent(m_rfilter->toString())
        );
    }

private:
    Vector2f m_invOutputSize;
    Transform m_sampleToCamera;
    Transform m_cameraToWorld;
    float m_fov;
    float m_nearClip;
    float m_farClip;
    float m_lensRadius;
    float m_focalDistance;
};

NORI_REGISTER_CLASS(ThinLensCamera, "thinlens");
NORI_NAMESPACE_END