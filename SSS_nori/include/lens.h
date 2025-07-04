//#pragma once
//
//#include <nori/camera.h>
//#include <nori/warp.h>
//
//NORI_NAMESPACE_BEGIN
//
//class ThinLensCamera : public Camera {
//public:
//    ThinLensCamera(const PropertyList &props) {
//        // Near plane size
//        m_size = props.getFloat("size", 35.0f);
//        
//        // Aspect ratio (width/height)
//        m_aspect = props.getFloat("aspect", 1.0f);
//        
//        // Aperture radius
//        m_lensRadius = props.getFloat("lensRadius", 0.0f);
//        
//        // Focal distance
//        m_focalDistance = props.getFloat("focalDistance", 1.0f);
//    }
//    Color3f sampleRay(Ray3f& ray,
//        const Point2f& samplePosition,
//        const Point2f& apertureSample) const override {
//        return Color3f(1);
//        }
//    Ray3f sampleRay(const Point2f &samplePosition, const Point2f &lensSample) {
//        // Convert raster position to normalized coordinates [-1, 1]
//        Point2f screenPosition = 2.0f * samplePosition - Point2f(1.0f);
//
//        // Map to near plane coordinates
//        Point3f nearPlanePosition(screenPosition.x() * m_size, screenPosition.y() * m_size / m_aspect, 0.0f);
//
//        // Ray direction (through near plane)
//        Vector3f rayDirection = (nearPlanePosition - m_cameraPosition).normalized();
//
//        // Thin lens sampling (only applies if lensRadius > 0)
//        if (m_lensRadius > 0.0f) {
//            // Sample a point on the lens
//            Point2f lensPoint = m_lensRadius * Warp::squareToUniformDiskConcentric(lensSample);
//            Point3f lensOrigin = m_cameraPosition + Vector3f(lensPoint.x(), lensPoint.y(), 0.0f);
//
//            // Compute intersection with the focal plane
//            float t_focal = m_focalDistance / rayDirection.z();
//            Point3f focalPoint = m_cameraPosition + t_focal * rayDirection;
//
//            // Update ray direction to pass through the focal point
//            rayDirection = (focalPoint - lensOrigin).normalized();
//
//            return Ray3f(lensOrigin, rayDirection);
//        }
//
//        // Default pinhole camera behavior
//        return Ray3f(m_cameraPosition, rayDirection);
//    }
//
//    std::string toString() const override {
//        return tfm::format(
//            "ThinLensCamera[\n"
//            "  size = %f,\n"
//            "  aspect = %f,\n"
//            "  lensRadius = %f,\n"
//            "  focalDistance = %f\n"
//            "]",
//            m_size, m_aspect, m_lensRadius, m_focalDistance
//        );
//    }
//
//private:
//    float m_size;          ///< Near plane size
//    float m_aspect;        ///< Aspect ratio
//    float m_lensRadius;    ///< Lens aperture radius
//    float m_focalDistance; ///< Distance to focal plane
//    Point3f m_cameraPosition = Point3f(0.0f, 0.0f, 0.0f); ///< Camera origin
//};
//
//NORI_REGISTER_CLASS(ThinLensCamera, "thinLens");
//
//NORI_NAMESPACE_END 