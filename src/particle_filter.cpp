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



void ParticleFilter::init(double x, double y, double theta, double std[]) {

	// Initialize all particles to first position (based on estimates of 
	//   x, y, theta and their uncertainties from GPS) and all weights to 1. 
	// Add random Gaussian noise to each particle.
        //Initialize randome engine to generate random numbers
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
        }
        is_initialized = true;

}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {

	// Add measurements to each particle and add random Gaussian noise.
        
	default_random_engine gen;

	// This line creates a normal (Gaussian) distribution for the velocity
	normal_distribution<double> dist_x(0, std_pos[0]);
	
	// Create normal distributions for y and psi
	normal_distribution<double> dist_y(0, std_pos[1]);
	normal_distribution<double> dist_theta(0, std_pos[2]);

        //prediction at t+1 for particle position and orientation in the case of yaw_rate=0
        if (yaw_rate==0) {
            
	   for (int i = 0; i < num_particles; ++i) {
		// sample particles from normal distributions
		// where "gen" is the random engine initialized earlier.
		particles[i].x += velocity * delta_t * cos(particles[i].theta) + dist_x(gen);
		particles[i].y += velocity * delta_t * sin(particles[i].theta) + dist_y(gen);
                particles[i].theta += dist_theta(gen);
           }
        }   


        //prediction at t+1 for particle position and orientation in the case of yaw_rate!=0
        else {
            
	   for (int i = 0; i < num_particles; ++i) {
		// sample particles from normal distributions
		// where "gen" is the random engine initialized earlier.
		particles[i].x += velocity/yaw_rate * (sin(particles[i].theta + yaw_rate * delta_t) - sin(particles[i].theta)) + dist_x(gen);
		particles[i].y += velocity/yaw_rate * (cos(particles[i].theta) - cos(particles[i].theta + yaw_rate * delta_t)) + dist_y(gen);
                particles[i].theta += yaw_rate * delta_t + dist_theta(gen);
           }
        }   

}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
     //method not used
}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
		std::vector<LandmarkObs> observations, Map map_landmarks) {

	//@Updates the weights of each particle using a multi-variate Gaussian distribution. 
	// The observations coming in are given in the VEHICLE'S coordinate system and are hence 
	// transformed into the MAP'S (global) coordinate system. 

        //Loop over every particle setting weight to 1 and checking which landmark fits best with observations
	for (int i = 0; i < num_particles; i++) {

            particles[i].weight=1.0;

            //For every particle loop over observations and check which landmark fits best to each observation
            // (Afterwards the multivariate Gaussian will compute which particles have the highest probability
            // of matching the "real" location and orientation of the vehicle based on the observations).
            for (int j = 0; j<observations.size(); ++j) { 

               //set x and y transformation for each observation to get global coordinates
               double x_trans = particles[i].x + observations[j].x * cos(particles[i].theta) - observations[j].y * sin(particles[i].theta);
               double y_trans = particles[i].y + observations[j].x * sin(particles[i].theta) + observations[j].y * cos(particles[i].theta);

               //using particle.id variable to associate to (nearest) map landmark
               particles[i].id = 0;
               double best_dist = dist(x_trans, y_trans, map_landmarks.landmark_list[0].x_f, map_landmarks.landmark_list[0].y_f);

               //search over all map landmarks to find closest/nearest to observation
               for (int k = 1; k < map_landmarks.landmark_list.size(); ++k) {

                   //only consider landmarks within sensor range
                   if (((map_landmarks.landmark_list[k].x_f - particles[i].x) < sensor_range) && ((map_landmarks.landmark_list[k].y_f - particles[i].y) < sensor_range))
                   {
                      double new_dist = dist(x_trans, y_trans, map_landmarks.landmark_list[k].x_f, map_landmarks.landmark_list[k].y_f);

                      if (new_dist < best_dist) 
                      {
                      best_dist = new_dist;
                      //assign landmark ID corresponding to nearest landmark
                      particles[i].id = k;
                      }
                   }
               }

               //compute particle weight based on multivariate Gaussian probability
               particles[i].weight *= 1/(2*M_PI*std_landmark[0]*std_landmark[1]) * exp(-(pow(x_trans - map_landmarks.landmark_list[particles[i].id].x_f, 2) 
                                                                                      + pow(y_trans - map_landmarks.landmark_list[particles[i].id].y_f, 2))
                                                                                        /(2 * M_PI *std_landmark[0]*std_landmark[1]));
            }

            weights[i] = particles[i].weight;
        }
}

void ParticleFilter::resample() {

	// @Resamples particles with replacement with probability proportional to their weight (using discrete integer distribution)

        default_random_engine gen;

        discrete_distribution<int> dist_d(weights.begin(), weights.end());


        for (int i=0; i < num_particles; ++i) {

           particles2[i] = particles[dist_d(gen)];

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
