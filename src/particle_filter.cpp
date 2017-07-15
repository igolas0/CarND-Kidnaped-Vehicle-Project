/*
 * particle_filter.cpp
 *
 *  Created on: Dec 12, 2016
 *      Author: Tiffany Huang
 */

#include <random>
#include <algorithm>
#include <iostream>
#include <numeric>
#include <math.h> 
#include <iostream>
#include <sstream>
#include <string>
#include <iterator>

#include "particle_filter.h"

using namespace std;

/*
 * Constructor.
 */
ParticleFilter::ParticleFilter() {
  
  //num_particles = 100;
  particles(num_particles);
  particles2(num_particles);
  weights(num_particles);
}

void ParticleFilter::init(double x, double y, double theta, double std[]) {
	// TODO: Set the number of particles. Initialize all particles to first position (based on estimates of 
	//   x, y, theta and their uncertainties from GPS) and all weights to 1. 
	// Add random Gaussian noise to each particle.
	// NOTE: Consult particle_filter.h for more information about this method (and others in this file).

	default_random_engine gen;
	 

	// This line creates a normal (Gaussian) distribution for x
	normal_distribution<double> dist_x(x, std[0]);
	
	// Create normal distributions for y and psi
	normal_distribution<double> dist_y(y, std[1]);
	normal_distribution<double> dist_theta(theta, std[2]);

	for (int i = 0; i < num_particles; ++i) {
		// sample particles from normal distributions
		// where "gen" is the random engine initialized earlier.
                particles[i].id = i;
		particles[i].x = dist_x(gen);
		particles[i].y = dist_y(gen);
		particles[i].theta = dist_theta(gen);	 
                particles[i].weight = 1.0;

        is_initialized = true;

}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/
        
	default_random_engine gen;

	// This line creates a normal (Gaussian) distribution for the veloxity
	normal_distribution<double> dist_x(0, std_pos[0]);
	
	// Create normal distributions for y and psi
	normal_distribution<double> dist_y(0, std[1]);
	normal_distribution<double> dist_theta(0, std[2]);

        if yaw_rate==0 {
            
	   for (int i = 0; i < num_particles; ++i) {
		// sample particles from normal distributions
		// where "gen" is the random engine initialized earlier.
		particles[i].x += velocity * delta_t * cos(particles[i].theta) + dist_x(gen);
		particles[i].y += velocity * delta_t * sin(particles[i].theta) + dist_y(gen);
                particles[i].theta += dist_theta(gen);
        }   


        else {
            
	   for (int i = 0; i < num_particles; ++i) {
		// sample particles from normal distributions
		// where "gen" is the random engine initialized earlier.
		particles[i].x += velocity/yaw_rate * (sin(particles[i].theta + yaw_rate * delta_t) - sin(particles[i].theta)) + dist_x(gen);
		particles[i].y += velocity/yaw_rate * (cos(particles[i].theta) - cos(particles[i].theta + yaw_rate * delta_t)) + dist_y(gen);
                particles[i].theta += dist_theta(gen);
        }   

}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the 
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to 
	//   implement this method and use it as a helper during the updateWeights phase.

}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
		std::vector<LandmarkObs> observations, Map map_landmarks) {
	// TODO: Update the weights of each particle using a mult-variate Gaussian distribution. You can read
	//   more about this distribution here: https://en.wikipedia.org/wiki/Multivariate_normal_distribution
	// NOTE: The observations are given in the VEHICLE'S coordinate system. Your particles are located
	//   according to the MAP'S coordinate system. You will need to transform between the two systems.
	//   Keep in mind that this transformation requires both rotation AND translation (but no scaling).
	//   The following is a good resource for the theory:
	//   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
	//   and the following is a good resource for the actual equation to implement (look at equation 
	//   3.33
	//   http://planning.cs.uiuc.edu/node99.html

	for (int i = 0; i < num_particles; i++) {
               

            particles[i].weight=1.0;

            for (int j = 0; j<observations.size(); ++j) { 
               

               //set x and y transformation for each observation to get global coordinates
               double x_trans = particles[i].x + observations[j].x * cos(particles[i].theta) - observations[j].y * sin(particles[i].theta);
               double y_trans = particles[i].y + observations[j].x * sin(particles[i].theta) + observations[j].y * cos(particles[i].theta);

               //using particle.id variable to associate to (nearest) map landmark
               particles[i].id = 0;
               double best_dist = dist(x_trans, y_trans, map_landmarks[0].x_f, map_landmarks[0].y_f]

               //search over all map landmarks to find closest/nearest to observation
               for (int x = 1; x<map_landmarks.size(); ++x) {

                   double new_dist = dist(x_trans, y_trans, map_landmarks[x].x_f, map_landmarks[x].y_f]

                   if (new_dist < best_dist) {

                   best_dist = new_dist
                   //assign landmark ID corresponding to nearest landmark
                   particles[i].id = x;

                   }

               }
               //compute particle weight based on multivariate Gaussian probability
               particles[i].weight *= 1/(2*M_PI*std_landmark[0]*std_landmark[1]) * exp(-(pow(x_trans - map_landmarks[particles[i].id].x_f, 2) 
                                                                                      + pow(y_trans - map_landmarks[particles[i].id].y_f, 2))
                                                                                        /(2 * M_PI *std_landmark[0]*std_landmark[1]));

            }
        }
}

void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight. 
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution

	for (int i = 0; i < num_particles; i++) {

          weights[i] = particles[i].weight;
      
        }
        //calculate max value in weights vector
        double mw = *max_element(begin(weights), std::end(weights));           

        default_random_engine gen;

        uniform_int_distribution <int> dist_n(0,num_particles);
        uniform_real_distribution <double> dist_mw(0.0, 2*mw);

        int index = dist_n(gen);
        double beta = 0.0;

        for (int i=0; i < num_particles; ++i) {

           beta += dist_mw(gen);
           
           while (beta > weights[index]) {

             beta -= weights[index];
             index = (index + 1)%num_particles;
           }
        particles2[i] = particles[index];
        }
        particles = particles2;
}

Particle ParticleFilter::SetAssociations(Particle particle, std::vector<int> associations, std::vector<double> sense_x, std::vector<double> sense_y)
{
	//particle: the particle to assign each listed association, and association's (x,y) world coordinates mapping to
	// associations: The landmark id that goes along with each listed association
	// sense_x: the associations x mapping already converted to world coordinates
	// sense_y: the associations y mapping already converted to world coordinates

	//Clear the previous associations
	particle.associations.clear();
	particle.sense_x.clear();
	particle.sense_y.clear();

	particle.associations= associations;
 	particle.sense_x = sense_x;
 	particle.sense_y = sense_y;

 	return particle;
}

string ParticleFilter::getAssociations(Particle best)
{
	vector<int> v = best.associations;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<int>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseX(Particle best)
{
	vector<double> v = best.sense_x;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseY(Particle best)
{
	vector<double> v = best.sense_y;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
