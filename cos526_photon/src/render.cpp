// Source file for photonmap renderer
// This implementation is a simple raycaster.
// Replace it with your own code.
  


////////////////////////////////////////////////////////////////////////
// Include files
////////////////////////////////////////////////////////////////////////

// #include "R3Graphics/R3Graphics.h"
// #include "fglut/fglut.h"
// #include <iostream>
// using namespace std;

#include "R3Graphics/R3Graphics.h"
#include "fglut/fglut.h"
#include "render.h"
#include<random>
#include<cmath>
#include<chrono>
#include <iostream>
using namespace std;

 
////////////////////////////////////////////////////////////////////////
// Function to render image with photon mapping
////////////////////////////////////////////////////////////////////////

// Handy global variables
enum SurfaceInteraction { DIFFUSE_REFLECTED, SPECULAR_REFLECTED, TRANSMITTED, ABSORBED };

int CLOSEST_PHOTONS_COUNT = 200;

//int MAX_BOUNCE_COUNT = 50;
int MAX_BOUNCE_COUNT = 500;
float P_RUSSIAN_ROULETTE = 0.5;

float MIN_DIST = 0;
float GLOBAL_MAX_DIST = 0.5;
float CAUSTIC_MAX_DIST = 0.02;
float MIN_BLEEDING_DIST = 0.1;

int RAYS_PER_PIXEL = 16; 
 
/* 
 
PHOTON TRACING:

[*****] Photon emission: Implement code to emit photons in random directions from every 
light source in a scene. The total number of photons emitted for each light source 
should be proportional to the power of the light source (so that each photon carries 
approximately equal power), and the distribution of photons should be proportional 
to the power in each direction -- e.g., for spot lights (section 2.1.1 in Jensen01).

[*****] Photon scattering: Trace photons via reflections and transmissions through the scene. 
At each ray-surface intersection, randomly generate a secondary ray along a direction 
of diffuse reflection, specular reflection, transmission, or absorption with probability 
proportional to kd, ks, kt, and (1 - kd+ks+kt), respectively (section 2.1.2 in Jensen01).

[*****] Russian Roulette: At each surface intersection, terminate rays with probability p 
(e.g., p=0.5) and multiply the power of surviving rays by 1.0/p (section 2.1.2 in Jensen01).
See section 8.5 of these siggraph course notes for details.

[*****] Photon storage: Store photon-surface intersections in a kd-tree, retaining the position, 
incident direction, and power of each photon. (section 2.1.3 in Jensen01). You can use your code 
from assignment 2, or the R3Kdtree class in R3Shapes to implement this feature.

[TODO] BRDF importance sampling: Select the directions of reflected and transmitted rays with probabilities 
proportional to the Phong BRDF of the scattering surface. See Jason Lawrence's notes for details.

[*****] Multiple photon maps: Implement separate photon maps for global (L{S|D}*D) and 
caustic (LS+D) ray paths (section 2.1.5 in Jensen01).


RENDERING:

[Need to do with importance sampling] Camera ray tracing: Generate a ray(s) from the camera eye point through each pixel. 
Trace them through the scene with reflections and transmissions at surface intersections -- 
i.e., at each ray-surface intersection, randomly generate a secondary ray along a direction 
of diffuse reflection, specular reflection, transmission, or absorption using importance sampling. 
(section 2.4 in Jensen01).

[*****] Radiance estimation: Use the kd-tree to find the N closest photons for each ray-surface intersection. 
Estimate the radiance traveling along the ray towards the camera from the power of those photons. 
(section 2.3.1 in Jensen01).

[TODO] Pixel integration: Trace multiple rays per pixel and average the radiance computed for all rays to 
estimate the radiance to store in the output image for each pixel. Compare the results with different 
numbers of rays per pixel (N).


VISUALIZATION: (helpful for debugging!)

Photon map visualization: Visualize photons stored in your photon map(s) -- e.g., show positions, 
normals, and powers of photons.

[TODO] Ray tracing visualization: Visualize ray paths traced from the camera -- e.g., show line segments 
between the camera and successive surface intersections for a random sampling of rays.

*/



/* PHOTON DATA STRUCTURE */

// Default Constructor
Photon::Photon(void)
{
  position  = R3Point(0, 0, 0);
  direction = R3Ray( R3Point(0, 0, 0), R3Point(1, 1, 1) );
  power     = RNRgb(0, 0, 0);
}

// Parametrized Constructor
Photon::Photon(R3Point _position, R3Ray _direction, RNRgb _power)
{
  position  = _position;
  direction = _direction;
  power     = _power;
}

// Member Functions() 
R3Point Photon::GetPosition()  { return position; } 
R3Ray   Photon::GetDirection() { return direction; } 
RNRgb   Photon::GetPower()     { return power; } 

void Photon::SetPosition(R3Point newPosition) { position = newPosition; } 
void Photon::SetDirection(R3Ray newDirection) { direction = newDirection; } 
void Photon::SetPower(RNRgb newPower)         { power = newPower; } 





/* Private helper function to generate a random number from [ min_value, max_value ) */ 
float 
GenerateRandomValue(float min_value, float max_value)
{
    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    std::mt19937 generator (seed);
    std::uniform_real_distribution<double> uniform01(min_value, max_value);

    return uniform01(generator);
}



double 
DistanceBetweenPoints(R3Point pt1, R3Point pt2) 
{
  double x = pt1.X() - pt2.X(); 
  double y = pt1.Y() - pt2.Y(); 
  double z = pt1.Z() - pt2.Z(); 

  double dist = sqrt(x*x + y*y + z*z);

  return dist; 
}




SurfaceInteraction
GetSurfaceInteraction( R3Scene *scene, 
                       R3Ray incident_ray)
{
  R3SceneNode *node;
  R3SceneElement *element;
  R3Shape *shape;
  R3Point point;
  R3Vector normal;
  RNScalar t;

  SurfaceInteraction surface_interaction = ABSORBED;

  scene->Intersects(incident_ray, &node, &element, &shape, &point, &normal, &t);

  const R3Material *material = (element) ? element->Material() : &R3default_material;
  const R3Brdf *brdf = (material) ? material->Brdf() : &R3default_brdf;

  // TODO: What if specular & transparent? Is this valid
  if      ( brdf->IsSpecular() )    { surface_interaction = SPECULAR_REFLECTED; }
  else if ( brdf->IsTransparent() ) { surface_interaction = TRANSMITTED; }
  else if ( brdf->IsDiffuse() )     { surface_interaction = DIFFUSE_REFLECTED; }

  return surface_interaction;
}








SurfaceInteraction
GetSurfaceInteraction( R3Scene *scene, 
                       Photon* Old_Photon,
                       RNRgb &new_photon_power)
{
  R3SceneNode *node;
  R3SceneElement *element;
  R3Shape *shape;
  R3Point point;
  R3Vector normal;
  RNScalar t;

  SurfaceInteraction surface_interaction;

  const R3Ray old_photon_ray = Old_Photon->GetDirection();

  double P_r = Old_Photon->GetPower().R();
  double P_g = Old_Photon->GetPower().G();
  double P_b = Old_Photon->GetPower().B(); 

  scene->Intersects(old_photon_ray, &node, &element, &shape, &point, &normal, &t);

  const R3Material *material = (element) ? element->Material() : &R3default_material;
  const R3Brdf *brdf = (material) ? material->Brdf() : &R3default_brdf;

  RNRgb Diffuse_Reflection_Coefficient  = brdf->Diffuse();
  RNRgb Specular_Reflection_Coefficient = brdf->Specular();
  RNRgb Transmission_Coefficient        = brdf->Transmission();

  double d_r = Diffuse_Reflection_Coefficient.R();
  double d_g = Diffuse_Reflection_Coefficient.G();
  double d_b = Diffuse_Reflection_Coefficient.B(); 

  double s_r = Specular_Reflection_Coefficient.R();
  double s_g = Specular_Reflection_Coefficient.G();
  double s_b = Specular_Reflection_Coefficient.B(); 

  double t_r = Transmission_Coefficient.R();
  double t_g = Transmission_Coefficient.G();
  double t_b = Transmission_Coefficient.B(); 

  double Max_Power_Incident_Ray = max( P_r, max( P_g, P_b ) );

  // TODO: Should all of these probabilities add up to 1?
  double P_Diffuse      = max( d_r * P_r, max( d_g * P_g, d_b * P_b ) ) / Max_Power_Incident_Ray;
  double P_Specular     = max( s_r * P_r, max( s_g * P_g, s_b * P_b ) ) / Max_Power_Incident_Ray;
  double P_Transmission = max( t_r * P_r, max( t_g * P_g, t_b * P_b ) ) / Max_Power_Incident_Ray;
  double P_Absorption   = 1 - ( P_Diffuse + P_Specular + P_Transmission );

  float p = GenerateRandomValue( 0, max( 1.0, P_Diffuse + P_Specular + P_Transmission ) );

  // Calculate the power of the incident photon, assume absorption by default
  double P_reflected_r = 0;
  double P_reflected_g = 0;
  double P_reflected_b = 0;

  if ( p < P_Diffuse ) // diffuse reflection
  {
    P_reflected_r = ( d_r * P_r ) / P_Diffuse;
    P_reflected_g = ( d_g * P_g ) / P_Diffuse;
    P_reflected_b = ( d_b * P_b ) / P_Diffuse;

    surface_interaction = DIFFUSE_REFLECTED;
  } 
  else if ( p < P_Diffuse + P_Specular ) // specular reflection
  {
    P_reflected_r = ( s_r * P_r ) / P_Specular;
    P_reflected_g = ( s_g * P_g ) / P_Specular;
    P_reflected_b = ( s_b * P_b ) / P_Specular;

    surface_interaction = SPECULAR_REFLECTED;
  } 
  else if ( p < P_Diffuse + P_Specular + P_Transmission ) // transmission through transparency
  {
    P_reflected_r = ( t_r * P_r ) / P_Transmission;
    P_reflected_g = ( t_g * P_g ) / P_Transmission;
    P_reflected_b = ( t_b * P_b ) / P_Transmission;

    surface_interaction = TRANSMITTED;
  } 
  else 
  {
    surface_interaction = ABSORBED;
  }

  new_photon_power = RNRgb( P_reflected_r, 
                            P_reflected_g, 
                            P_reflected_b);

  return surface_interaction;
}







R3Vector
SpecularReflectionBRDFImportanceSampling(const R3Brdf brdf, 
                                         R3Point intersection_point, 
                                         R3Vector normal, 
                                         R3Vector reflected_ray)
{
  // Random variables
  // float u1 = GenerateRandomValue(-100, 100);
  // float u2 = GenerateRandomValue(-100, 100);
  // float u3 = GenerateRandomValue(-100, 100);

  float u1 = GenerateRandomValue(0, 1);
  float u2 = GenerateRandomValue(0, 1);

  // Specular exponent (higher values give a sharper specular reflection)
  RNScalar n = brdf.Shininess();

  // Transformation in spherical coordinates with respect to direction of perfect specular reflection
  float alpha = acos( powf( u1, 1 / (n+1) ) );
  float phi = 2 * M_PI * u2;

  // float x = cos(phi) * cos(theta); 
  // float y = sin(phi) * cos(theta);
  // float z = sin(theta);

  float x = cos(phi) * cos(alpha); 
  float y = sin(phi) * cos(alpha);
  float z = sin(alpha);

  return reflected_ray + R3Vector(x, y, z);
  //return reflected_ray;
  //return R3Vector(x, y, z); 
}








R3Vector
DiffuseReflectionBRDFImportanceSampling(const R3Brdf brdf, 
                                         R3Point intersection_point, 
                                         R3Vector normal, 
                                         R3Vector incident_ray)
{
  // Random variables
  float u1 = GenerateRandomValue(0, 1);
  float u2 = GenerateRandomValue(0, 1);

  // Transformation in spherical coordinates with respect to direction of perfect specular reflection
  float theta = acos( sqrt(u1) );
  float phi = 2 * M_PI * u2;

  float x = cos(phi) * cos(theta); 
  float y = sin(phi) * cos(theta);
  float z = sin(theta);

  return normal + R3Vector(x, y, z);

  //return incident_ray;
}











// // const R3Brdf& brdf, const R3Point& eye, 
// //     const R3Point& point, const R3Vector& normal

/* BRDF IMPORTANCE SAMPLING */

RNRgb 
GetBRDF(R3Brdf brdf, R3Point origin, R3Point intersection_point, R3Vector normal) 
{
  RNRgb rgb;

  float u1 = GenerateRandomValue(0, 1);
  float u2 = GenerateRandomValue(0, 1);

  // Get material properties
  const RNRgb Dc = brdf.Diffuse();
  const RNRgb Sc = brdf.Specular();

  // alpha = the angle between the perfect specular reflective direction and the outgoing ray direction; 
  // values larger than PI/2 are clamped to PI/2 in order to prevent any negative values of the cosine factor

  // k_d = the diffuse reflectivity (the fraction of the incoming energy that is reflected diffusely)

  // k_s = the specular reflectivity (the fraction of the perpendicularly incoming energy that is reflected specularly)

  // n = the specular exponent (higher values give a sharper specular reflection)
  RNScalar n = brdf.Shininess();

  // x = point = the point of intersection with the current ray and the scene

  // N = normal = the vector normal to the differential surface patch at the point of intersection
  R3Vector N = R3Vector(normal); 
  N.Normalize();

  // w_o = a unit vector in the direction of the current outgoing ray of light 
  // (recall that we are tracing paths backwards from the direction of light propagation: from the eye to the light sources)

  // f_r(x, theta_i, theta_o) = the Phong reflectance BRDF with diffuse / specular parameters given by the material of the intersected surface

  // L = direction from light to point
  R3Vector incident_vector = origin - intersection_point;
  incident_vector.Normalize();

  // theta_i = angle between normal and incident ray (ray from light to point)
  RNScalar NL = N.Dot(incident_vector);
  float theta_i = acos(NL);

  // Calculate direction of reflected vector
  // r = perfect specular reflective direction
  R3Vector r = Reflect(incident_vector, normal);
  R3Ray new_photon_ray = R3Ray(intersection_point - incident_vector * 0.00001, r);

  return rgb;
}


    // // Get material properties
    // const RNRgb& Dc = brdf.Diffuse();
    // const RNRgb& Sc = brdf.Specular();
    // RNScalar s = brdf.Shininess();

    // // Get light properties
    // RNScalar I = IntensityAtPoint(point);
    // R3Vector L = DirectionFromPoint(point);
    // const RNRgb& Ic = Color();

    // // Compute geometric stuff
    // RNScalar NL = normal.Dot(L);
    // if (RNIsNegativeOrZero(NL)) return RNblack_rgb;
    // R3Vector R = (2.0 * NL) * normal - L;
    // R3Vector V = eye - point;
    // V.Normalize();
    // RNScalar VR = V.Dot(R);

    // // Compute diffuse reflection
    // RNRgb rgb = (I * NL) * Dc * Ic;

    // // Compute specular reflection
    // if (RNIsPositive(VR)) rgb += (I * pow(VR,s)) * Sc * Ic;

    // // Return total reflection
    // return rgb;



  // if (brdf) {
  //   //Loop over all lights
  //   for (int k = 0; k < scene->NLights(); k++) {
  //     R3Light *light = scene->Light(k);
  //     color += light->Reflection(*brdf, eye, reference_point, normal);
  //   }
  // }

















/* PHOTON EMISSION */

Photon* 
ScatterPhoton(R3Scene *scene, 
              Photon* Old_Photon,
              RNRgb new_photon_power,
              SurfaceInteraction incident_surface_interaction,
              RNRgb &next_photon_power,
              SurfaceInteraction &new_surface_interaction)
{
  R3SceneNode *node;
  R3SceneElement *element;
  R3Shape *shape;
  R3Point point;
  R3Vector normal;
  RNScalar t;

  const R3Point old_photon_position = Old_Photon->GetPosition();
  const R3Ray   old_photon_ray      = Old_Photon->GetDirection(); 

  scene->Intersects(old_photon_ray, &node, &element, &shape, &point, &normal, &t);

  const R3Material *material = (element) ? element->Material() : &R3default_material;
  const R3Brdf *brdf = (material) ? material->Brdf() : &R3default_brdf;

  // Calculate direction of reflected vector
  R3Vector r = Reflect(old_photon_ray.Vector(), normal);
  R3Vector r_sampled = DiffuseReflectionBRDFImportanceSampling(*brdf, point, normal, r);
  R3Ray new_photon_ray = R3Ray(old_photon_position - old_photon_ray.Vector() * 0.001, r_sampled);
  //R3Ray new_photon_ray = R3Ray(point, r_sampled);

  /* RUSSIAN ROULETTE: 
  At each surface intersection, terminate rays with probability p 
  (e.g., p=0.5) and multiply the power of surviving rays by 1.0/p (section 2.1.2 in Jensen01).
  See section 8.5 of these siggraph course notes for details. */

  float p_Termination = GenerateRandomValue(0, 1.0);

  if (p_Termination < P_RUSSIAN_ROULETTE)
  {
    // Terminate photon in Russian Roulette by setting power to 0
    new_photon_power = RNRgb(0, 0, 0);
  }
  else 
  {
    // Calculate secondary ray along direction of diffuse reflection, specular reflection, transmission, or absorption 

    /* BRDF IMPORTANCE SAMPLING:   [TODO]
    Select the directions of reflected and transmitted rays with probabilities proportional 
    to the Phong BRDF of the scattering surface. See Jason Lawrence's notes for details. */

    if ( incident_surface_interaction == DIFFUSE_REFLECTED )
    {
      // Diffuse reflection
      //new_photon_ray = R3Ray(old_photon_position - old_photon_ray.Vector() * 0.00001, r);

      R3Vector reflected_ray = Reflect(old_photon_ray.Vector(), normal);
      R3Vector r_sampled = DiffuseReflectionBRDFImportanceSampling(*brdf, point, normal, reflected_ray);
      new_photon_ray = R3Ray(old_photon_position - old_photon_ray.Vector() * 0.001, r_sampled);
    } 
    else if ( incident_surface_interaction == SPECULAR_REFLECTED )
    {
      // Specular reflection
      //new_photon_ray = R3Ray(old_photon_position - old_photon_ray.Vector() * 0.00001, r);

      R3Vector reflected_ray = Reflect(old_photon_ray.Vector(), normal);
      R3Vector r_sampled = SpecularReflectionBRDFImportanceSampling(*brdf, point, normal, reflected_ray);
      new_photon_ray = R3Ray(old_photon_position - old_photon_ray.Vector() * 0.001, r_sampled);
    } 
    else if ( incident_surface_interaction == TRANSMITTED )
    {
      // Transmission through transparency
      float index_of_refraction = brdf->IndexOfRefraction();
      R3Vector new_photon_direction = Refract(old_photon_ray.Vector(), normal, index_of_refraction);
      new_photon_ray = R3Ray(old_photon_position, new_photon_direction);
    } 

    float Russian_Roulette_factor = 1.0 / max( P_RUSSIAN_ROULETTE, 0.1f ) ;
    new_photon_power *= Russian_Roulette_factor;
  }

  R3Point new_photon_position = R3Point(0, 0, 0);
  if (!scene->Intersects(new_photon_ray, &node, &element, &shape, &new_photon_position, &normal, &t)) { return NULL; }

  Photon* New_Photon = new Photon(new_photon_position, 
                                  new_photon_ray,
                                  new_photon_power); 
  // Set the surface interaction
  new_surface_interaction = GetSurfaceInteraction(scene, 
                                                  New_Photon,
                                                  next_photon_power);
  return New_Photon;
}









/* Main method for creating the PHOTON MAP */

R3Kdtree<Photon *>
CreatePhotonMap(R3Scene *scene, 
                int photon_count,
                bool isCausticMap,
                RNArray<Photon *>* All_Photons)
{
  // Convenient variables
  //const R3Point& eye = scene->Camera().Origin();
  R3SceneNode *node;
  R3SceneElement *element;
  R3Shape *shape;
  R3Point point;
  R3Vector normal;
  RNScalar t;

  /* PHOTON EMISSION: 
  Implement code to emit photons in random directions from every light source in a scene. 
  The total number of photons emitted for each light source should be proportional to 
  the power of the light source (so that each photon carries approximately equal power), 
  and the distribution of photons should be proportional to the power in each direction 
  -- e.g., for spot lights (section 2.1.1 in Jensen01). */

  int total_intensity = 0;                      // Sum of max intensities of all lights in the scene                                                                                                                                                                         
  int light_max_intensities [scene->NLights()]; // Max intensity for respective light

  // Iterate through every light in the scene to find total light intensity
  for (int i = 0; i < scene->NLights(); i++) 
  {
    R3Light *light = scene->Light(i);
    RNRgb color = light->Color();

    int rIntensity = color.R();
    int gIntensity = color.G();
    int bIntensity = color.B();

    light_max_intensities[i] = max( rIntensity, max( gIntensity, bIntensity ) );
    total_intensity += light_max_intensities[i];
  }

  int photon_index = 0;

  // Iterate through every light in the scene again and this time emit photons
  for (int i = 0; i < scene->NLights(); i++) 
  {
    R3Light *light = scene->Light(i);
    RNRgb color = light->Color();

    int start_photon_index = photon_index;
    int photon_allotment = round( ( (float) light_max_intensities[i] / total_intensity ) * photon_count );
    
    // Iterate through all the photons in photon_count to emit photons from light sources
    for (int j = start_photon_index; j < start_photon_index + photon_allotment; j++)
    {
      // Initializing photon position and power variables with dummies values
      R3Point photon_position = R3Point(11, 11, 11);
      RNRgb   photon_power    = RNRgb(10, 10, 10);

      R3Ray photon_ray;

      R3Material *material;
      const R3Brdf *brdf;

      while (photon_position == R3Point(11, 11, 11)) 
      {
        // Generate a random ray from the light to determine direction vector of photon from light
        // Uniform in all directions as a default for the first part of this assignment
        // All lights currently treated like point lights
        photon_ray = light->RandomlySampledRay();

        if ( scene->Intersects(photon_ray, &node, &element, &shape, &point, &normal, &t) ) 
        {
          material = (element) ? element->Material() : &R3default_material;
          brdf     = (material) ? material->Brdf() : &R3default_brdf;

          // For caustic maps, check that initial contact is a spot for specular reflection or transmission 
          if ( isCausticMap ) 
          {
            if ( brdf->IsSpecular() || 
                 brdf->IsTransparent() ) 
            {
              photon_position = R3Point( point.X(), 
                                         point.Y(), 
                                         point.Z() );
            } // otherwise, disregard this photon
          }
          else { 
            photon_position = R3Point( point.X(), 
                                       point.Y(), 
                                       point.Z() );
          }
        }
      }

      photon_power = light->PowerGivenDistance(photon_position, normal);

      float clamping_factor = 1.0f / max(photon_power.R(), max(photon_power.G(), photon_power.B()));
      if (clamping_factor < 1.0) { photon_power *= clamping_factor; }

      //photon_power = photon_power / 10;

      // Make a photon with previously calculated position, direction, and power values
      Photon* Current_Photon = new Photon(photon_position, 
                                          photon_ray,
                                          photon_power); 

      RNRgb new_photon_power;
      SurfaceInteraction incident_surface_interaction = GetSurfaceInteraction( scene, 
                                                                               Current_Photon,
                                                                               new_photon_power );

      //Store photon if at a diffuse surface
      // if ( incident_surface_interaction == DIFFUSE_REFLECTED && !isCausticMap ) 
      // {
      //   All_Photons->InsertKth(Current_Photon, photon_index);
      //   photon_index++;
      // }

      // All_Photons->InsertKth(Current_Photon, photon_index);
      // photon_index++;

      int bounce_count = 0;



      while ( Current_Photon->GetPower() != RNRgb(0,0,0) && 
              bounce_count < MAX_BOUNCE_COUNT )
      {
        RNRgb next_photon_power;
        SurfaceInteraction new_surface_interaction;

        Photon* Scattered_Photon = ScatterPhoton( scene, 
                                                  Current_Photon, 
                                                  new_photon_power, 
                                                  incident_surface_interaction, 
                                                  next_photon_power,
                                                  new_surface_interaction );

        if ( Scattered_Photon != NULL && 
             Scattered_Photon->GetPower() != RNRgb(0,0,0) ) 
        {

          All_Photons->InsertKth(Scattered_Photon, photon_index);
          Current_Photon = Scattered_Photon;
          photon_index++;
          bounce_count++;

          if ( isCausticMap && 
             ( new_surface_interaction == SPECULAR_REFLECTED || new_surface_interaction == TRANSMITTED ) ) 
          {
            // store scattered photon & continue propogating the light chain
            Current_Photon = Scattered_Photon;
            bounce_count++;
            new_photon_power = next_photon_power;
            incident_surface_interaction = new_surface_interaction;
          } 

          else if ( isCausticMap && 
                    new_surface_interaction == DIFFUSE_REFLECTED ) 
          {
            if ( bounce_count != 0 ) // store scattered photon and terminate the light chain
            {
              All_Photons->InsertKth(Scattered_Photon, photon_index);
              bounce_count = MAX_BOUNCE_COUNT;
              photon_index++; 
            } 
            else { bounce_count = MAX_BOUNCE_COUNT; } // disregard scattered photon and terminate light chain
          }

          else if (!isCausticMap) 
          {
            // Store photon if at a diffuse surface and not too close to the corner to prevent color bleeding
            if ( new_surface_interaction == DIFFUSE_REFLECTED && 
                 DistanceBetweenPoints(Current_Photon->GetPosition(), Scattered_Photon->GetPosition()) > MIN_BLEEDING_DIST ) 
            {
              All_Photons->InsertKth(Scattered_Photon, photon_index);
              photon_index++;
            }
            Current_Photon = Scattered_Photon;
            bounce_count++;
            new_photon_power = next_photon_power;
            incident_surface_interaction = new_surface_interaction;
          }
          else { bounce_count = MAX_BOUNCE_COUNT; }
        }

        else { bounce_count = MAX_BOUNCE_COUNT; }  // disregard scattered photon and terminate the light chain

      }
    }
  }

  /* PHOTON STORAGE: 
  Store photon-surface intersections in a kd-tree, retaining the position, 
  incident direction, and power of each photon. (section 2.1.3 in Jensen01). 
  Let's use R3Kdtree class in R3Shapes to implement this feature. */

  // PHOTON MAP, STORING PHOTONS IN KD-TREE
  R3Kdtree<Photon*> Photon_Map(*All_Photons);

  return Photon_Map;
}




/* RADIANCE ESTIMATION: 
Use the kd-tree to find the N closest photons for each ray-surface intersection. 
Estimate the radiance traveling along the ray towards the camera from the power 
of those photons. (section 2.3.1 in Jensen01). */

RNRgb
EstimateRadiance( R3Scene* scene, 
                  R3Point reference_point,
                  const R3Brdf *brdf,
                  R3Vector normal,
                  int closest_points_count, 
                  bool isCausticMap, 
                  R3Kdtree<Photon*> Photon_Map )
{
  const R3Point& eye = scene->Camera().Origin();
  double MAX_DIST = (isCausticMap) ? CAUSTIC_MAX_DIST : GLOBAL_MAX_DIST;
  RNArray<Photon *> Closest_Photons;

  int total_photon_count = Photon_Map.FindClosest( reference_point, 
                                                   MIN_DIST, 
                                                   MAX_DIST, 
                                                   closest_points_count, 
                                                   Closest_Photons );
  RNRgb total_power = RNRgb(0, 0, 0);

  if ((total_photon_count < 0.05 * closest_points_count) && isCausticMap) { return RNRgb(0, 0, 0); }

  float k = 1000; // cone filter constant characterizing the filter, k >= 1

  for (int i = 0; i < total_photon_count; i++) 
  {
    Photon *Current_Photon = Closest_Photons.Kth(i);

    float d_p = powf( reference_point.X() - Current_Photon->GetPosition().X(), 2 ) + 
                powf( reference_point.Y() - Current_Photon->GetPosition().Y(), 2 ) + 
                powf( reference_point.Z() - Current_Photon->GetPosition().Z(), 2 );

    // CONE FILTERING 
    // Weight of cone filter
    float w_cone = 1 - ( d_p / (k * MAX_DIST) );

    // GAUSSIAN FILTERING
    // Weight of Gaussian filter
    float alpha = 0.918; 
    float beta  = 1.953; 
    float w_gaussian = alpha * ( 1 - ( (1-exp(-beta*((d_p*d_p)/(2*MAX_DIST*MAX_DIST)))) / (1-exp(-beta)) ));

    total_power += Current_Photon->GetPower() * w_cone;
    //total_power += Current_Photon->GetPower() * w_gaussian;
  }

  total_power = total_power / total_photon_count;

  // Area from which the photons are gathered
  double surface_area = M_PI * MAX_DIST * MAX_DIST;

  // Nora Willett write-up: 
  // To calculate radiance, sum over the power times the BRDF divided by the area of that the photons were gathered from

  RNRgb color = RNRgb(0, 0, 0);

  if (brdf) {
    //Loop over all lights
    for (int k = 0; k < scene->NLights(); k++) {
      R3Light *light = scene->Light(k);
      color += light->Reflection(*brdf, eye, reference_point, normal);
    }
  }

  total_power = total_power / 4;

  // Check that color channels are not maxed out, clamp if so

  if (!isCausticMap) {
    float clamping_factor = 0.25f / max(color.R(), max(color.G(), color.B()));
    if (clamping_factor < 1.0) { color *= clamping_factor; }
  }

  RNRgb radiance = (total_power * color) / surface_area;
  //RNRgb radiance = (total_power) / surface_area;

  // Cone filter normalization
  float cone_filter_normalization_factor = 1 - ( 2 / (3*k) );
  radiance = radiance / cone_filter_normalization_factor;

  // Gaussian filter needs no normalization

  return radiance;
}













/* CAMERA RAY TRACING: 
Generate a ray(s) from the camera eye point through each pixel. Trace them through 
the scene with reflections and transmissions at surface intersections -- i.e., at each 
ray-surface intersection, randomly generate a secondary ray along a direction of diffuse 
reflection, specular reflection, transmission, or absorption using importance sampling. 
(section 2.4 in Jensen01). */

R2Image *
RenderImagePhotonMapping(R3Scene *scene, 
                         int width, 
                         int height, 
                         int photon_count,
                         R3Kdtree<Photon *> Global_Photon_Map, 
                         R3Kdtree<Photon *> Caustic_Photon_Map)
{
  // Allocate image
  R2Image *image = new R2Image(width, height);
  if (!image) {
    fprintf(stderr, "Unable to allocate image\n");
    return NULL;
  }

  // Convenient variables
  const R3Point& eye = scene->Camera().Origin();
  R3SceneNode *node;
  R3SceneElement *element;
  R3Shape *shape;
  R3Point point;
  R3Vector normal;
  RNScalar t;

  // R3Kdtree<Photon *> Photon_Map = Global_Photon_Map;

  // RNArray<Photon *> Caustic_Photons;
  // int caustic_photon_count = Photon_Map.FindAll( R3Point(0, 0, 0), 0, 10000, Caustic_Photons );

  // for (int i = 0; i < Caustic_Photon_Map.NPoints(); i++)
  // {
  //   Photon_Map.InsertPoint(Caustic_Photons.Kth(i));
  // }

  /* CAMERA RAY TRACING: 
  Generate a ray(s) from the camera eye point through each pixel. Trace them through 
  the scene with reflections and transmissions at surface intersections -- i.e., at each 
  ray-surface intersection, randomly generate a secondary ray along a direction of diffuse 
  reflection, specular reflection, transmission, or absorption using importance sampling. 
  (section 2.4 in Jensen01). */

  for (int i = 0; i < width; i++) {
    for (int j = 0; j < height; j++) {

      R3Ray ray = scene->Viewer().WorldRay(i, j);

      /* At each ray-surface intersection, randomly generate a secondary ray along a direction 
      of diffuse reflection, specular reflection, transmission, or absorption using importance sampling. */

      if (scene->Intersects(ray, &node, &element, &shape, &point, &normal, &t)) {

        // Get intersection information
        const R3Material *material = (element) ? element->Material() : &R3default_material;
        const R3Brdf *brdf = (material) ? material->Brdf() : &R3default_brdf;
        bool inShadow = FALSE;

        // Compute color
        RNRgb color = scene->Ambient();
        if (brdf) 
        {
          color += brdf->Emission();

          //Loop over all lights
          for (int k = 0; k < scene->NLights(); k++) {
            R3Light *light = scene->Light(k);

            R3Ray Light_To_Point_Ray = light->LightToPointRay(point);

            R3SceneNode    *_node;
            R3SceneElement *_element;
            R3Shape        *_shape;
            R3Point         _point;
            R3Vector        _normal;
            RNScalar        _t;

            //color += light->Reflection(*brdf, eye, point, normal);

            // Adding shadow rays
            if ( scene->Intersects(Light_To_Point_Ray, &_node, &_element, &_shape, &_point, &_normal, &_t) )
            {
              double diff_X = point.X() - _point.X(); 
              double diff_Y = point.Y() - _point.Y(); 
              double diff_Z = point.Z() - _point.Z(); 

              double dist = sqrt( diff_X * diff_X + 
                                  diff_Y * diff_Y + 
                                  diff_Z * diff_Z );

              if (dist < 0.005) { color += light->Reflection(*brdf, eye, point, normal); }
              else { 
                SurfaceInteraction surface_interaction = GetSurfaceInteraction(scene, Light_To_Point_Ray);
                if (!brdf->IsTransparent()) { inShadow = TRUE; }
              }
            } 
            else { color += light->Reflection(*brdf, eye, point, normal); }
          }
        }

        color = color / 2;



        /* PIXEL INTEGRATION: Trace multiple rays per pixel and average the radiance computed for 
        all rays to estimate the radiance to store in the output image for each pixel. 
        Compare the results with different numbers of rays per pixel (N). */

        // Determine the surface interaction to determine properties of secondary ray
        SurfaceInteraction surface_interaction = GetSurfaceInteraction( scene, ray );
        RNRgb radiance = RNRgb(0, 0, 0);
        bool assignedRadiance = FALSE;

          //ray = scene->Viewer().WorldRay(i, j);
          //surface_interaction = GetSurfaceInteraction( scene, ray );

        while ( surface_interaction != DIFFUSE_REFLECTED )
        {
          R3Material *material = (element) ? element->Material() : &R3default_material;
          const R3Brdf *brdf   = (material) ? material->Brdf() : &R3default_brdf;

          if ( surface_interaction == TRANSMITTED)
          {
            float index_of_refraction = brdf->IndexOfRefraction();
            R3Vector new_photon_direction = Refract(ray.Vector(), normal, index_of_refraction);
            ray = R3Ray(point + new_photon_direction * 0.00001, new_photon_direction);
          }
          else if ( surface_interaction == SPECULAR_REFLECTED )
          {
            // Calculate direction of reflected vector
            R3Vector reflected_ray = Reflect(ray.Vector(), normal);
            //R3Vector r_sampled = SpecularReflectionBRDFImportanceSampling(*brdf, point, normal, reflected_ray);
            //ray = R3Ray(point - ray.Vector() * 0.001, r_sampled);
            ray = R3Ray(point - ray.Vector() * 0.001, reflected_ray);
          }
          else if ( surface_interaction == ABSORBED )  { break; }

          if (!scene->Intersects(ray, &node, &element, &shape, &point, &normal, &t)) { break; } //TODO: make sure color rendered at pixel is black
          surface_interaction = GetSurfaceInteraction( scene, ray );

          if ( surface_interaction == DIFFUSE_REFLECTED ) 
          {
            R3Material *material = (element) ? element->Material() : &R3default_material;
            const R3Brdf *brdf   = (material) ? material->Brdf() : &R3default_brdf;
            radiance += EstimateRadiance(scene, point, brdf, normal, CLOSEST_PHOTONS_COUNT, FALSE, Global_Photon_Map);
            assignedRadiance = TRUE;
          }
        }

        if ( !assignedRadiance ) 
        {
          for (int k = 0; k < RAYS_PER_PIXEL; k++) 
          {
            if ( surface_interaction == DIFFUSE_REFLECTED ) 
            {
              R3Material *material = (element) ? element->Material() : &R3default_material;
              const R3Brdf *brdf   = (material) ? material->Brdf() : &R3default_brdf;

              R3Vector reflected_ray = Reflect(ray.Vector(), normal);
              R3Vector r_sampled = DiffuseReflectionBRDFImportanceSampling(*brdf, point, normal, reflected_ray);
              ray = R3Ray(point - ray.Vector() * 0.001, r_sampled);

              if (scene->Intersects(ray, &node, &element, &shape, &point, &normal, &t)) 
              {
                radiance += EstimateRadiance(scene, point, brdf, normal, CLOSEST_PHOTONS_COUNT, FALSE, Global_Photon_Map);
              }
            }
          }
          radiance = radiance / RAYS_PER_PIXEL;
          //radiance = radiance / 4;
        }

        //radiance = radiance / 4;

        //if (inShadow) { radiance = radiance / 5; }
        color += radiance;




        RNRgb caustic_radiance = EstimateRadiance(scene, point, brdf, normal, CLOSEST_PHOTONS_COUNT, TRUE, Caustic_Photon_Map);
        color += caustic_radiance;


        // Check that color channels are not maxed out, clamp if so
        float clamping_factor = 1.0f / max(color.R(), max(color.G(), color.B()));
        if (clamping_factor < 1.0) { color *= clamping_factor; }

        // Set pixel color
        image->SetPixelRGB(i, j, color);
      }
    }
  }

  // Return image
  return image;
}









R2Image *
RenderImageRaytracing(R3Scene *scene, 
  int width, int height, 
  int print_verbose)
{
  // Start statistics
  RNTime start_time;
  start_time.Read();
  int ray_count = 0;

  // Allocate image
  R2Image *image = new R2Image(width, height);
  if (!image) {
    fprintf(stderr, "Unable to allocate image\n");
    return NULL;
  }

  // Convenient variables
  const R3Point& eye = scene->Camera().Origin();
  R3SceneNode *node;
  R3SceneElement *element;
  R3Shape *shape;
  R3Point point;
  R3Vector normal;
  RNScalar t;

  // Draw intersection point and normal for some rays
  for (int i = 0; i < width; i++) {
    for (int j = 0; j < height; j++) {
      R3Ray ray = scene->Viewer().WorldRay(i, j);

      // intersection point
      if (scene->Intersects(ray, &node, &element, &shape, &point, &normal, &t)) {

        // Get intersection information
        const R3Material *material = (element) ? element->Material() : &R3default_material;
        const R3Brdf *brdf = (material) ? material->Brdf() : &R3default_brdf;

        // Compute color
        RNRgb color = scene->Ambient();
        if (brdf) 
        {
          color += brdf->Emission();

          //Loop over all lights
          for (int k = 0; k < scene->NLights(); k++) {
            R3Light *light = scene->Light(k);

            R3Ray Light_To_Point_Ray = light->LightToPointRay(point);

            R3SceneNode    *_node;
            R3SceneElement *_element;
            R3Shape        *_shape;
            R3Point         _point;
            R3Vector        _normal;
            RNScalar        _t;

            //color += light->Reflection(*brdf, eye, point, normal);

            // Adding shadow rays
            if ( scene->Intersects(Light_To_Point_Ray, &_node, &_element, &_shape, &_point, &_normal, &_t) )
            {
              double diff_X = point.X() - _point.X(); 
              double diff_Y = point.Y() - _point.Y(); 
              double diff_Z = point.Z() - _point.Z(); 

              double dist = sqrt( diff_X * diff_X + 
                                  diff_Y * diff_Y + 
                                  diff_Z * diff_Z );

              if (dist < 0.01) { color += light->Reflection(*brdf, eye, point, normal); }
            } 
            else { color += light->Reflection(*brdf, eye, point, normal); }
          }
        }

        // Need secondary rays, shadowing
        // similar for second phase of photon mapping
        // Dataa structure 
        // Loop over all lights in a scene and emit photons

        color = color / 2;

        // Check that color channels are not maxed out, clamp if so
        float clamping_factor = 1.0f / max(color.R(), max(color.G(), color.B()));
        if (clamping_factor < 1.0) { color *= clamping_factor; }

        // Set pixel color
        image->SetPixelRGB(i, j, color);

        // Update ray count
        ray_count++;
      }
    }
  }

  // Print statistics
  if (print_verbose) {
    printf("Rendered image ...\n");
    printf("  Time = %.2f seconds\n", start_time.Elapsed());
    printf("  # Rays = %d\n", ray_count);
    fflush(stdout);
  }

  // Return image
  return image;
}



