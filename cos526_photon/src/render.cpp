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

/* 
 
PHOTON TRACING:

[*****] Photon emission: Implement code to emit photons in random directions from every 
light source in a scene. The total number of photons emitted for each light source 
should be proportional to the power of the light source (so that each photon carries 
approximately equal power), and the distribution of photons should be proportional 
to the power in each direction -- e.g., for spot lights (section 2.1.1 in Jensen01).

Photon scattering: Trace photons via reflections and transmissions through the scene. 
At each ray-surface intersection, randomly generate a secondary ray along a direction 
of diffuse reflection, specular reflection, transmission, or absorption with probability 
proportional to kd, ks, kt, and (1 - kd+ks+kt), respectively (section 2.1.2 in Jensen01).

Russian Roulette: At each surface intersection, terminate rays with probability p 
(e.g., p=0.5) and multiply the power of surviving rays by 1.0/p (section 2.1.2 in Jensen01).
See section 8.5 of these siggraph course notes for details.

[*****] Photon storage: Store photon-surface intersections in a kd-tree, retaining the position, 
incident direction, and power of each photon. (section 2.1.3 in Jensen01). You can use your code 
from assignment 2, or the R3Kdtree class in R3Shapes to implement this feature.

BRDF importance sampling: Select the directions of reflected and transmitted rays with probabilities 
proportional to the Phong BRDF of the scattering surface. See Jason Lawrence's notes for details.

Multiple photon maps: Implement separate photon maps for global (L{S|D}*D) and 
caustic (LS+D) ray paths (section 2.1.5 in Jensen01).


RENDERING:

Camera ray tracing: Generate a ray(s) from the camera eye point through each pixel. 
Trace them through the scene with reflections and transmissions at surface intersections -- 
i.e., at each ray-surface intersection, randomly generate a secondary ray along a direction 
of diffuse reflection, specular reflection, transmission, or absorption using importance sampling. 
(section 2.4 in Jensen01).

[*****] Radiance estimation: Use the kd-tree to find the N closest photons for each ray-surface intersection. 
Estimate the radiance traveling along the ray towards the camera from the power of those photons. 
(section 2.3.1 in Jensen01).

Pixel integration: Trace multiple rays per pixel and average the radiance computed for all rays to 
estimate the radiance to store in the output image for each pixel. Compare the results with different 
numbers of rays per pixel (N).


VISUALIZATION: (helpful for debugging!)

Photon map visualization: Visualize photons stored in your photon map(s) -- e.g., show positions, 
normals, and powers of photons.

Ray tracing visualization: Visualize ray paths traced from the camera -- e.g., show line segments 
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






/* PHOTON EMISSION */

Photon* 
ScatterPhoton(R3Scene *scene, 
              Photon* Old_Photon)
{
  R3SceneNode *node;
  R3SceneElement *element;
  R3Shape *shape;
  R3Point point;
  R3Vector normal;
  RNScalar t;

  enum Path_Type { DIFFUSE_REFLECTED, SPECULAR_REFLECTED, TRANSMITTED, ABSORBED };

  const R3Point old_photon_position = Old_Photon->GetPosition();
  const R3Ray   old_photon_ray      = Old_Photon->GetDirection(); 

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

  // cout << "Probability, Diffuse: " << P_Diffuse << endl;
  // cout << "Probability, Specular: " << P_Specular << endl;
  // cout << "Probability, Transmission: " << P_Transmission << endl;
  // cout << "Probability, Absorption: " << P_Absorption << endl;

  float p = GenerateRandomValue( 0, max( 1.0, P_Diffuse + P_Specular + P_Transmission ) );

  // Calculate the power of the incident photon
  // Assume absorption by default
  double P_reflected_r = 0;
  double P_reflected_g = 0;
  double P_reflected_b = 0;

  Path_Type light_path_type = ABSORBED;

  if ( p < P_Diffuse ) // diffuse reflection
  {
    P_reflected_r = ( d_r * P_r ) / P_Diffuse;
    P_reflected_g = ( d_g * P_g ) / P_Diffuse;
    P_reflected_b = ( d_b * P_b ) / P_Diffuse;

    light_path_type = DIFFUSE_REFLECTED;
  } 
  else if ( p < P_Diffuse + P_Specular ) // specular reflection
  {
    P_reflected_r = ( s_r * P_r ) / P_Specular;
    P_reflected_g = ( s_g * P_g ) / P_Specular;
    P_reflected_b = ( s_b * P_b ) / P_Specular;

    light_path_type = SPECULAR_REFLECTED;
  } 
  else if ( p < P_Diffuse + P_Specular + P_Transmission ) // transmission through transparency
  {
    P_reflected_r = ( t_r * P_r ) / P_Transmission;
    P_reflected_g = ( t_g * P_g ) / P_Transmission;
    P_reflected_b = ( t_b * P_b ) / P_Transmission;

    light_path_type = TRANSMITTED;
  } // otherwise, absorption

  RNRgb new_photon_power = RNRgb( P_reflected_r, 
                                  P_reflected_g, 
                                  P_reflected_b );

  /* BRDF IMPORTANCE SAMPLING: 
  Select the directions of reflected and transmitted rays with probabilities proportional 
  to the Phong BRDF of the scattering surface. See Jason Lawrence's notes for details. */

  //R3Point old_photon_position = R3Point(point); 

  // r = d - 2 (d * n) n, where (d * n) is the dot product and n is normalized
  // https://math.stackexchange.com/questions/13261/how-to-get-a-reflection-vector

  R3Vector d = R3Vector(old_photon_ray.Vector());
  R3Vector n = R3Vector(normal);

  n.Normalize();
  float d_n = d.Dot(n);

  R3Vector r = (d -  2 * d_n * n);

  //R3Ray new_photon_ray = R3Ray(old_photon_position, old_photon_ray.Vector());
  R3Ray new_photon_ray = R3Ray(old_photon_position - old_photon_ray.Vector() * 0.00001, r);

  //R3Plane surface_plane       = R3Plane(old_photon_position, normal);
  
  //surface_plane.Align(normal);

  // if (light_path_type == DIFFUSE_REFLECTED){
  //   new_photon_ray = R3Ray(old_photon_position, r);
  // }
  // else if (light_path_type == SPECULAR_REFLECTED) {
  //   new_photon_ray = R3Ray(old_photon_position, r);
  // }
  // else if (light_path_type == TRANSMITTED) {
  //   //http://asawicki.info/news_1301_reflect_and_refract_functions.html
  //   //float index_of_refraction = brdf->IndexOfRefraction();
  //   new_photon_ray = R3Ray(old_photon_position, old_photon_ray.Vector());
  // }

  // Calculate secondary ray along direction of diffuse reflection, specular reflection, transmission, or absorption 
  //R3Ray new_photon_ray = R3Ray( old_photon_position, new_photon_position );
  R3Point new_photon_position = R3Point(0, 0, 0);

  R3SceneNode *node_;
  R3SceneElement *element_;
  R3Shape *shape_;
  R3Vector normal_;
  RNScalar t_;

  scene->Intersects(new_photon_ray, &node_, &element_, &shape_, &new_photon_position, &normal_, &t_);

  // if ( new_photon_position.X() != old_photon_position.X() &&
  //      new_photon_position.Y() != old_photon_position.Y() &&
  //      new_photon_position.Z() != old_photon_position.Z() ) 
  // {
  //   cout << "DIDN'T INTERSECT" << endl;
  // }
  //scene->Intersects(new_photon_ray, &node_, &element_, &shape_, &new_photon_position, &normal_, &t_);



  // cout << "New Position x: " << new_photon_position.X() << endl;
  // cout << "New Position y: " << new_photon_position.Y() << endl;
  // cout << "New Position z: " << new_photon_position.Z() << endl;

  //R3Point new_photon_position = R3Point(point);

  Photon* New_Photon = new Photon(new_photon_position, 
                                  new_photon_ray,
                                  new_photon_power); 

  return New_Photon;

}









/* PHOTON EMISSION */

R3Kdtree<Photon *>
EmitPhotons(R3Scene *scene, 
            int photon_count,
            RNArray<Photon *>* All_Photons)
{
  // Convenient variables
  // const R3Point& eye = scene->Camera().Origin();
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

  // Iterate through every light in the scene again and this time emeit photons
  for (int i = 0; i < scene->NLights(); i++) 
  {
    R3Light *light = scene->Light(i);
    RNRgb    color = light->Color();

    int start_photon_index = photon_index;
    int photon_allotment = round( ( (float) light_max_intensities[i] / total_intensity ) * photon_count );
    
    // Iterate through all the photons in photon_count
    for (int j = start_photon_index; j < start_photon_index + photon_allotment; j++)
    {
      // Initializing photon position and power variables with dummies values
      R3Point photon_position = R3Point(11, 11, 11);
      RNRgb   photon_power    = RNRgb(10, 10, 10);

      R3Ray photon_ray;

      photon_power = RNRgb( color.R(), 
                            color.G(), 
                            color.B() );

      while (photon_position == R3Point(11, 11, 11)) 
      {
        // Generate a random ray from the light to determine direction vector of photon from light
        // Uniform in all directions as a default for the first part of this assignment
        // All lights currently treated like point lights
        photon_ray = light->RandomlySampledRay();

        if ( scene->Intersects(photon_ray, &node, &element, &shape, &point, &normal, &t) ) 
        {
          photon_position = R3Point( point.X(), 
                                     point.Y(), 
                                     point.Z() );
        }
      }

      // Make a photon with previously calculated position, direction, and power values
      Photon* photon = new Photon(photon_position, 
                                  photon_ray,
                                  photon_power); 

      //Photon* New_Photon = EmitPhotonsFromLight(photon_index)
      All_Photons->InsertKth(photon, j);
      photon_index++;
    }

  }

  // For all unassign photons assign a dummy photon
  for (int i = photon_index; i < photon_count; i++)
  {
    Photon* photon = new Photon();
    All_Photons->InsertKth(photon, i);
    photon_index++; 
  }


  /* PHOTON SCATTERING: 
  Trace photons via reflections and transmissions through the scene. At each 
  ray-surface intersection, randomly generate a secondary ray along a direction 
  of diffuse reflection, specular reflection, transmission, or absorption with 
  probability proportional to kd, ks, kt, and (1 - kd+ks+kt), respectively 
  (section 2.1.2 in Jensen01). */

  for (int i = 0; i < photon_count; i++)
  {
    Photon* Current_Photon = All_Photons->Kth(i);

    if ( Current_Photon->GetPower().R() != 0 && 
         Current_Photon->GetPower().G() != 0 && 
         Current_Photon->GetPower().B() != 0 )
    {
      Photon* Scattered_Photon = ScatterPhoton(scene, Current_Photon);
      All_Photons->InsertKth(Scattered_Photon, photon_count + i);
      photon_index++; 
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
EstimateRadiance(R3Point reference_point,
                 int closest_points_count, 
                 RNArray<Photon *> Closest_Photons, 
                 R3Kdtree<Photon*> Photon_Map)
{
    float MIN_DIST = 0;
    float MAX_DIST = 0.005;

    // int FindClosest(const R3Point& query_position, 
    // RNLength min_distance, RNLength max_distance, int max_points, 
    // RNArray<PtrType>& points, RNLength *distances = NULL) const;

    int total_photon_count = Photon_Map.FindClosest(reference_point, 
                                                    MIN_DIST, 
                                                    MAX_DIST, 
                                                    closest_points_count, 
                                                    Closest_Photons);

    // cout << "GREETINGS FROM THE OTHER SIDE " << endl;
    // cout << "Closest Photon: " << Closest_Photons.Kth(0)->GetPosition().X() << endl;
    // cout << "Closest Photon: " << Closest_Photons.Kth(1)->GetPosition().X() << endl;
    // cout << "Closest Photon: " << Closest_Photons.Kth(2)->GetPosition().X() << endl;
    // cout << "Final: " << total_photon_count << endl;

    double R_total = 0;
    double G_total = 0;
    double B_total = 0;

    for (int i = 0; i < total_photon_count; i++) 
    {
      R_total += Closest_Photons.Kth(i)->GetPower().R();
      G_total += Closest_Photons.Kth(i)->GetPower().G();
      B_total += Closest_Photons.Kth(i)->GetPower().B();
    }

    // Return the average color 
    return RNRgb( R_total / total_photon_count, 
                  R_total / total_photon_count, 
                  B_total / total_photon_count );
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
                         R3Kdtree<Photon *> Photon_Map)
{
  // Allocate image
  R2Image *image = new R2Image(width, height);
  if (!image) {
    fprintf(stderr, "Unable to allocate image\n");
    return NULL;
  }

  // Convenient variables
  // const R3Point& eye = scene->Camera().Origin();
  R3SceneNode *node;
  R3SceneElement *element;
  R3Shape *shape;
  R3Point point;
  R3Vector normal;
  RNScalar t;

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

        RNArray<Photon *> Closest_Photons;

        int CLOSEST_PHOTONS_COUNT = 5;

        RNRgb color = EstimateRadiance(point, CLOSEST_PHOTONS_COUNT, Closest_Photons, Photon_Map);

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
        if (brdf) {
          color += brdf->Emission();

          //Loop over all lights
          for (int k = 0; k < scene->NLights(); k++) {
            R3Light *light = scene->Light(k);
            color += light->Reflection(*brdf, eye, point, normal);
          }
        }

        // Need secondary rays, shadowing
        // similar for second phase of photon mapping
        // Dataa structure 
        // Loop over all lights in a scene and emit photons


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



