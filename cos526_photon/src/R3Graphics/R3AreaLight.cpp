/* Source file for the R3 light class */



/* Include files */

#include "R3Graphics.h"
#include<random>
#include<cmath>
#include<chrono>
#include <iostream>
using namespace std;


/* Public variables */

R3AreaLight R3null_area_light;



/* Class type definitions */

RN_CLASS_TYPE_DEFINITIONS(R3AreaLight);



/* Public functions */

int 
R3InitAreaLight()
{
    /* Return success */
    return TRUE;
}



void 
R3StopAreaLight()
{
}



R3AreaLight::
R3AreaLight(void)
{
}



R3AreaLight::
R3AreaLight(const R3AreaLight& light)
    : R3Light(light),
      circle(light.circle),
      constant_attenuation(light.constant_attenuation),
      linear_attenuation(light.linear_attenuation),
      quadratic_attenuation(light.quadratic_attenuation)
{
}



R3AreaLight::
R3AreaLight(const R3Point& position, RNLength radius, const R3Vector& direction, const RNRgb& color,
	        RNScalar intensity, RNBoolean active,
            RNScalar ca, RNScalar la, RNScalar qa)
    : R3Light(color, intensity, active),
      circle(position, radius, direction),
      constant_attenuation(ca),
      linear_attenuation(la),
      quadratic_attenuation(qa)
{
}




void R3AreaLight::
SetPosition(const R3Point& position)
{
    // Set position
    this->circle.Reset(position, Radius(), Direction());
}



void R3AreaLight::
SetDirection(const R3Vector& direction)
{
    // Set direction
    this->circle.Reset(Position(), Radius(), direction);
}



void R3AreaLight::
SetRadius(RNLength radius)
{
    // Set radius
    this->circle.Reset(Position(), radius, Direction());
}



void R3AreaLight::
SetConstantAttenuation(RNScalar ca)
{
  // Set constant coefficient of attenuation
  this->constant_attenuation = ca;
}



void R3AreaLight::
SetLinearAttenuation(RNScalar la)
{
  // Set linear coefficient of attenuation
  this->linear_attenuation = la;
}



void R3AreaLight::
SetQuadraticAttenuation(RNScalar qa)
{
  // Set quadratic coefficient of attenuation
  this->quadratic_attenuation = qa;
}



RNRgb R3AreaLight::
DiffuseReflection(const R3Brdf& brdf, 
    const R3Point& point, const R3Vector& normal) const
{
    // Parameter ????
    const int max_samples = 4;

    // Check if light is active
    if (!IsActive()) return RNblack_rgb;

    // Get material properties
    const RNRgb& Dc = brdf.Diffuse();

    // Get light properties
    const RNRgb& Ic = Color();

    // Get circle axes
    R3Vector direction = circle.Normal();
    RNDimension dim = direction.MinDimension();
    R3Vector axis1 = direction % R3xyz_triad[dim];
    axis1.Normalize();
    R3Vector axis2 = direction % axis1;
    axis2.Normalize();

    // Sample points on light source
    int sample_count = 0;
    RNRgb sample_sum = RNblack_rgb;
    while (sample_count < max_samples) {
        // Sample point in circle
        RNScalar r1 = RNRandomScalar();
        RNScalar r2 = RNRandomScalar();
        if (r1*r1 + r2*r2 > 1) continue;
        R3Point sample_point = Position();
        sample_point += r1 * Radius() * axis1;
        sample_point += r2 * Radius() * axis2;
        sample_count++;
        
        // Compute intensity at point
        RNScalar I = Intensity();
        RNLength d = R3Distance(point, sample_point);
        RNScalar denom = constant_attenuation;
        denom += d * linear_attenuation;
        denom += d * d * quadratic_attenuation;
        if (RNIsPositive(denom)) I /= denom;

        // Compute direction at point
        R3Vector L = sample_point - point;
        L.Normalize();

        // Compute diffuse reflection from sample point
        RNScalar NL = normal.Dot(L);
        if (RNIsNegativeOrZero(NL)) continue;
        RNRgb diffuse = (I * NL) * Dc * Ic;

        // Add to result
        sample_sum += diffuse;
    }

    // Return diffuse reflection from area
    RNArea area = circle.Area();
    RNRgb sample_mean = sample_sum;
    if (sample_count > 0) sample_mean /= sample_count;
    return area * sample_mean;
}



RNRgb R3AreaLight::
SpecularReflection(const R3Brdf& brdf, const R3Point& eye, 
    const R3Point& point, const R3Vector& normal) const
{
    // Parameter ????
    const int max_samples = 8;

    // Check if light is active
    if (!IsActive()) return RNblack_rgb;

    // Get material properties
    const RNRgb& Sc = brdf.Specular();
    RNScalar s = brdf.Shininess();

    // Get light properties
    const RNRgb& Ic = Color();

    // Get circle axes
    R3Vector direction = circle.Normal();
    RNDimension dim = direction.MinDimension();
    R3Vector axis1 = direction % R3xyz_triad[dim];
    axis1.Normalize();
    R3Vector axis2 = direction % axis1;
    axis2.Normalize();

    // Sample points on light source
    int sample_count = 0;
    RNRgb sample_sum = RNblack_rgb;
    for (int i = 0; i < max_samples; i++) {
        // Sample point in circle
        RNScalar r1 = RNRandomScalar();
        RNScalar r2 = RNRandomScalar();
        if (r1*r1 + r2*r2 > 1) continue;
        R3Point sample_point = Position();
        sample_point += r1 * Radius() * axis1;
        sample_point += r2 * Radius() * axis2;
        sample_count++;

        // Compute intensity at point
        RNScalar I = Intensity();
        RNLength d = R3Distance(point, sample_point);
        RNScalar denom = constant_attenuation;
        denom += d * linear_attenuation;
        denom += d * d * quadratic_attenuation;
        if (RNIsPositive(denom)) I /= denom;

        // Compute direction at point
        R3Vector L = sample_point - point;
        L.Normalize();

        // Compute specular reflection from sample_point
        RNScalar NL = normal.Dot(L);
        if (RNIsNegativeOrZero(NL)) return RNblack_rgb;
        R3Vector R = (2.0 * NL) * normal - L;
        R3Vector V = eye - point;
        V.Normalize();
        RNScalar VR = V.Dot(R);
        if (RNIsNegativeOrZero(VR)) continue;

        // Return specular component of reflection
        sample_sum += (I * pow(VR,s)) * Sc * Ic;
    }

    // Return specular reflection from area
    RNArea area = circle.Area();
    RNRgb sample_mean = sample_sum;
    if (sample_count > 0) sample_mean /= sample_count;
    return area * sample_mean;
}



RNRgb R3AreaLight::
Reflection(const R3Brdf& brdf, const R3Point& eye, 
    const R3Point& point, const R3Vector& normal) const
{
    // Return total reflection
    RNRgb diffuse = DiffuseReflection(brdf, point, normal);
    RNRgb specular = SpecularReflection(brdf, eye, point, normal);
    return diffuse + specular;
}


// Get ray from light source to given point
const R3Ray R3AreaLight::
LightToPointRay(R3Point point) const
{
    return R3Ray(Position(), point);
}




// Give a randomly sampled array from a point light
const R3Ray R3AreaLight::
RandomlySampledRay(void) const
{
    // Get circle axes
    R3Vector direction = circle.Normal();
    RNDimension dim = direction.MinDimension();
    R3Vector axis1 = direction % R3xyz_triad[dim];
    axis1.Normalize();
    R3Vector axis2 = direction % axis1;
    axis2.Normalize();

    // Sample point uniformly from circular light source
    RNScalar r1 = RNRandomScalar() * 2 - 1;
    RNScalar r2 = RNRandomScalar() * 2 - 1; 
    
    while (r1*r1 + r2*r2 > 1) {
        r1 = RNRandomScalar() * 2 - 1;
        r2 = RNRandomScalar() * 2 - 1;
    }

    R3Point sample_point = Position();
    sample_point += r1 * Radius() * axis1;
    sample_point += r2 * Radius() * axis2;

    R3Ray ray = R3Ray(sample_point, Direction());

    return ray;
}


const RNRgb R3AreaLight::
PowerGivenDistance(R3Point reference_point, R3Vector normal) const
{
    RNRgb resulting_power; 

    double x = reference_point.X() - Position().X();
    double y = reference_point.Y() - Position().Y();
    double z = reference_point.Z() - Position().Z();

    // Distance between the two points
    double d = sqrt( x*x + y*y + z*z );

    // Difference in angle between the light direction and the direction to the surface
    R3Vector A = normal; 
    R3Vector B = Direction();

    double Ax = A.X(); 
    double Ay = A.Y();
    double Az = A.Z();
    double Bx = B.X(); 
    double By = B.Y();
    double Bz = B.Z();

    double A_dot_B = A.Dot(B);
    double magnitude_A = sqrt( Ax*Ax + Ay*Ay + Az*Az );
    double magnitude_B = sqrt( Bx*Bx + By*By + Bz*Bz );

    // angle = inverse cosine (arccos) of dot product divided by multiplied magnitudes
    double _theta_ = acos( A_dot_B / ( magnitude_A * magnitude_B ) ); // angle between surface normal and direction of light
    
    // Plus PI/2, then negate that result 
    double theta  = -( _theta_ + ( M_PI / 2 ) ); // angle between surface and direction of light

    double ca = ConstantAttenuation();
    double la = LinearAttenuation();
    double qa = QuadraticAttenuation();

    // Calculate: (r,g,b) * 1.0/ (ca + la*d + qa*d*d)
    resulting_power = Color() * ( cos(theta) / ( ca + (la*d) + (qa*d*d) ) );

    return resulting_power;
}



void R3AreaLight::
Draw(int i) const
{
    // Draw light
    GLenum index = (GLenum) (GL_LIGHT2 + i);
    if (index > GL_LIGHT7) return;
    GLfloat buffer[4];
    buffer[0] = Intensity() * Color().R();
    buffer[1] = Intensity() * Color().G();
    buffer[2] = Intensity() * Color().B();
    buffer[3] = 1.0;
    glLightfv(index, GL_DIFFUSE, buffer);
    glLightfv(index, GL_SPECULAR, buffer);
    buffer[0] = Position().X();
    buffer[1] = Position().Y();
    buffer[2] = Position().Z();
    buffer[3] = 1.0;
    glLightfv(index, GL_POSITION, buffer);
    buffer[0] = 90.0;
    glLightf(index, GL_SPOT_CUTOFF, buffer[0]);
    glEnable(index);
    buffer[0] = ConstantAttenuation();
    buffer[1] = LinearAttenuation();
    buffer[2] = QuadraticAttenuation();
    glLightf(index, GL_CONSTANT_ATTENUATION, buffer[0]);
    glLightf(index, GL_LINEAR_ATTENUATION, buffer[1]);
    glLightf(index, GL_QUADRATIC_ATTENUATION, buffer[2]);
}



