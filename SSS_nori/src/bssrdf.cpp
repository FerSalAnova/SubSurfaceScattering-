#include <nori/bsdf.h>
#include <nori/frame.h>
#include <nori/warp.h>
#include <nori/reflectance.h>
#include <cmath>
#include <vector>
#include <random>

NORI_NAMESPACE_BEGIN

class BSSRDF : public BSDF {
public:
    BSSRDF(const PropertyList &props) {
        // Default values for subsurface scattering coefficients
        m_sigma_a = props.getColor("sigma_a", Color3f(0.1f));  // Absorption coefficient
        m_sigma_s = props.getColor("sigma_s", Color3f(1.0f));  // Scattering coefficient
        m_eta = props.getFloat("eta", 1.5f);  // Refractive index
        m_alpha = props.getFloat("alpha", 0.1f);  // Roughness parameter
        m_albedo = props.getColor("albedo", Color3f(0.8f));  // Albedo for surface color
    }

    /// Evaluate the BSSRDF for a given pair of directions
    Color3f eval(const BSDFQueryRecord &bRec) const override {
        if (bRec.measure != ESolidAngle)
            return Color3f(0.0f);

        // Compute the scattering contribution using single scattering and diffusion
        Color3f singleScattering = computeSingleScattering(bRec);
        Color3f diffusion = computeDiffusion(bRec);

        // Combine the results (scattering + diffusion), applying albedo as a factor
        return m_albedo * (singleScattering + diffusion);

    }
    

    /// PDF for BSSRDF evaluation
    float pdf(const BSDFQueryRecord &bRec) const override {
        if (Frame::cosTheta(bRec.wi) <= 0 || Frame::cosTheta(bRec.wo) <= 0)
            return 0.0f;

        // Lambertian PDF: cos(theta) / pi
        return Frame::cosTheta(bRec.wo) / M_PI;
    }

    /// Sample the BSSRDF
    Color3f sample(BSDFQueryRecord &bRec, const Point2f &sample) const override {
        if (Frame::cosTheta(bRec.wi) <= 0)
            return Color3f(0.0f);

        // Cosine-weighted sampling over the hemisphere
        Vector3f wo = Warp::squareToCosineHemisphere(sample);
        bRec.wo = wo;

        //// Ensure that the outgoing direction is valid
        if (Frame::cosTheta(bRec.wo) <= 0)
            return Color3f(0.0f);

        // Handle extreme grazing angles by clamping cosTheta or modifying Fresnel terms
        if (Frame::cosTheta(bRec.wi) < 0.01f || Frame::cosTheta(bRec.wo) < 0.01f) {
            return Color3f(0.0f);  // Ignore extreme grazing angles where scattering is invalid
        }

        bRec.measure = ESolidAngle;

        return eval(bRec) * Frame::cosTheta(bRec.wo) / pdf(bRec);
    }

    const char* type() const override { return "BSSRDF"; }

    /// Whether the material is diffuse
    bool isDiffuse() const override {
        return true;
    }

    /// Debugging utility to print the properties of the BSSRDF
    std::string toString() const override {
        return tfm::format(
            "BSSRDF[\n"
            "  sigma_a = %s,\n"
            "  sigma_s = %s,\n"
            "  eta = %f,\n"
            "  alpha = %f,\n"
            "  albedo = %s\n"
            "]",
            m_sigma_a.toString(), m_sigma_s.toString(), m_eta, m_alpha, m_albedo.toString()
        );
    }

protected:
    Color3f computeSingleScattering(const BSDFQueryRecord &bRec) const {
        Color3f sigma_t = m_sigma_a + m_sigma_s;
        float attenuation = exp(-sigma_t.average());  // Apply average attenuation for simplicity
        float cosThetaI = Frame::cosTheta(bRec.wi); 
        float Ft = 1.0f - Reflectance::fresnel(cosThetaI, 1.0f, m_eta);  // Fresnel term for single scattering
        float phase = Reflectance::henyeyGreenstein(bRec.wi.dot(bRec.wo), -0.1);  // Phase function

        // Implementing per-channel absorption and scattering for more realistic effects
        Color3f result = m_sigma_s * Ft * phase * attenuation;

        // Apply albedo to the result
        return m_albedo * result;
        //return result;
    }

    // Enhanced diffusion model
    Color3f computeDiffusion(const BSDFQueryRecord &bRec) const {
        Color3f sigma_t = m_sigma_a + m_sigma_s;
        float dr = (bRec.wi - bRec.wo).norm();
        if (dr < 1e-4f) {
        return Color3f(0.0f);  // Avoid computation for extremely small distances
    }
        Color3f sigma_tr = sqrt(3.0f * m_sigma_a * (m_sigma_a + m_sigma_s));
        
        // Diffusion term based on the Henyey-Greenstein phase function
        Color3f term1 = exp(-sigma_t * dr) / (4.0f * M_PI * dr);
        Color3f term2 = exp(-sigma_t * sqrt(dr * dr + sigma_tr * sigma_tr)) / 
                        (4.0f * M_PI * sqrt(dr * dr + sigma_tr * sigma_tr));

        // Apply albedo to the diffusion result
        return m_albedo * m_sigma_s * (term1 - term2);
        //return  m_sigma_s * (term1 - term2);
    }

    /// Sample the radial distance for multiple scattering
    float sampleDistance(float sigma_tr, Sampler *sampler) const {
        float xi = sampler->next1D();
        return -std::log(1.0f - xi) / sigma_tr;
    }

private:
    Color3f m_sigma_a;  // Absorption coefficient
    Color3f m_sigma_s;  // Scattering coefficient
    float m_eta;        // Refractive index
    float m_alpha;      // Roughness parameter
    Color3f m_albedo;   // Albedo for surface color
};

NORI_REGISTER_CLASS(BSSRDF, "bssrdf");
NORI_NAMESPACE_END