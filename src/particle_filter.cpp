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
	// TODO: Set the number of particles. Initialize all particles to first position (based on estimates of 
	//   x, y, theta and their uncertainties from GPS) and all weights to 1. 
	// Add random Gaussian noise to each particle.
	// NOTE: Consult particle_filter.h for more information about this method (and others in this file).
	// std::cout << "INIT" << std::endl;

	// Create default random generator
	default_random_engine gen;
	
	// Set number of particle
	// There doesn't seem to be any benefit above 100
	num_particles = 100;

	// create normal distributions to draw from
	normal_distribution<double> dist_x(x, std[0]);
	normal_distribution<double> dist_y(y, std[1]);
	normal_distribution<double> dist_theta(theta, std[2]);

	// set initial particles based on initial input
	for (int i = 0; i < num_particles; i++) {
		Particle p;
		p.id = i;
		p.x = dist_x(gen);
		p.y = dist_y(gen);
		p.theta = dist_theta(gen);
		p.weight = 1.0;

		particles.push_back(p);
		weights.push_back(1.0);

		// std::cout << "[" << p.id << "] x:" << p.x << " y:" << p.y << " theta:" << p.theta << std::endl;	
	}

	is_initialized = true;

}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/

	// Generator for Guassian sensor noise
	default_random_engine gen;

	// Iterate through all particles
	for (std::vector<Particle>::iterator i = particles.begin(); i != particles.end(); i++) {

		// Handle yaw_rate == 0 
		if (yaw_rate == 0.0) {
			(*i).x = (*i).x + velocity*delta_t*cos((*i).theta);
			(*i).y = (*i).y + velocity*delta_t*sin((*i).theta);
		} else {
			(*i).x = (*i).x + (velocity/yaw_rate)*(sin((*i).theta + (yaw_rate*delta_t)) - sin((*i).theta));
			(*i).y = (*i).y + (velocity/yaw_rate)*(cos((*i).theta) - cos((*i).theta + yaw_rate*delta_t));
			(*i).theta = (*i).theta + yaw_rate*delta_t;
		}


		// add Gaussian sensor noise
		normal_distribution<double> dist_x((*i).x, std_pos[0]);
		normal_distribution<double> dist_y((*i).y, std_pos[1]);
		normal_distribution<double> dist_theta((*i).theta, std_pos[2]);

		(*i).x = dist_x(gen);
		(*i).y = dist_y(gen);
		(*i).theta = dist_theta(gen);

	}

}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the 
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to 
	//   implement this method and use it as a helper during the updateWeights phase.

	// Perform nearest neighbor association for observations vs. pruned map landmark list
	for (std::vector<LandmarkObs>::iterator obs = observations.begin(); obs != observations.end(); obs++) {
		// set an initial very large distance
		double closest_dist = 100000000.0;
		double curr_dist = 0.0;
		int closest_id = -1;

		// search through map landmarks for closest match
		for (std::vector<LandmarkObs>::iterator pred = predicted.begin(); pred != predicted.end(); pred++) {
			//curr_dist = sqrt((pow(((*pred).x - (*obs).x), 2.0) + pow(((*pred).y - (*obs).y), 2.0)));
			curr_dist = dist((*pred).x, (*pred).y, (*obs).x, (*obs).y);
			if (curr_dist < closest_dist) {
				closest_dist = curr_dist;
				closest_id = (*pred).id;
			}
		}

		(*obs).id = closest_id;
	}
}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
		const std::vector<LandmarkObs> &observations, const Map &map_landmarks) {
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


	// std::cout << "UPDATE WEIGHTS for " << observations.size() << "observations" << std::endl;

	// foreach particle
	//for (std::vector<Particle>::iterator p = particles.begin(); p != particles.end(); p++) {
	for (int p_i = 0; p_i < num_particles; p_i++) {

		double theta = particles[p_i].theta;
		// not sure if this is (p.x,p.y) or (p.sense_x, p.sense_y)
		double xp = particles[p_i].x;
		double yp = particles[p_i].y;

		// Convert observation vehicle coordinates to map coordinates using Homogenous Transform
		std::vector<LandmarkObs> t_obs;
		for (std::vector<LandmarkObs>::const_iterator obs = observations.begin(); obs != observations.end(); obs++) {
			
			LandmarkObs o;
			double xc = (*obs).x;
			double yc = (*obs).y;

			o.id = (*obs).id;
			o.x = xc*cos(theta) - yc*sin(theta) + xp;
			o.y = xc*sin(theta) + yc*cos(theta) + yp;

			t_obs.push_back(o);
		}	

		// create vector of predicted landmarks (filters by sensor distance)
		std::vector<LandmarkObs> predicted_landmarks;
		for (std::vector<Map::single_landmark_s>::const_iterator lm = map_landmarks.landmark_list.begin(); lm != map_landmarks.landmark_list.end(); lm++) {
			double xl = (double)(*lm).x_f;
			double yl = (double)(*lm).y_f;
			if (dist(xp, yp, xl, yl) <= sensor_range) {
				LandmarkObs l;
				l.id = (*lm).id_i;
				l.x = xl;
				l.y = yl;

				predicted_landmarks.push_back(l);
			}
		}

		// Create Associations
		dataAssociation(predicted_landmarks, t_obs);

		// Set associations for debugging (will show blue lines in simulator)
		std::vector<int> associations;
		std::vector<double> a_sense_x;
		std::vector<double> a_sense_y;
		for (unsigned int i = 0; i < t_obs.size(); i++) {

			associations.push_back(t_obs[i].id);
			a_sense_x.push_back(t_obs[i].x);
			a_sense_y.push_back(t_obs[i].y);
		}
		SetAssociations(particles[p_i], associations, a_sense_x, a_sense_y);

		// calculate new weight
		std::vector<double> obs_weights;
		for (std::vector<LandmarkObs>::iterator obs = t_obs.begin(); obs != t_obs.end(); obs++) {
			// cache some values for convenience
			double x_obs = (*obs).x;
			double y_obs = (*obs).y;
			double mu_x = (double)map_landmarks.landmark_list[(*obs).id - 1].x_f;
			double mu_y = (double)map_landmarks.landmark_list[(*obs).id - 1].y_f;
			double sig_x = std_landmark[0];
			double sig_y = std_landmark[1];
			
			// set intermediate terms
			double gauss_norm = 1/(2*M_PI*sig_x*sig_y);
			double exponent = pow(x_obs - mu_x, 2.0)/(2*pow(sig_x, 2.0)) + pow(y_obs - mu_y, 2.0)/(2*pow(sig_y, 2.0));

			// calculate weight
			double new_weight = gauss_norm*exp(-1.0*exponent);

			obs_weights.push_back(new_weight);
		}
			

		// multiply all weights together
		double weight = std::accumulate(obs_weights.begin(), obs_weights.end(), 1.0, std::multiplies<double>());
		// Sum all weights for normalization
		double weight_norm = std::accumulate(obs_weights.begin(), obs_weights.end(), 0.0);

		// set particle weight and weight list
		particles[p_i].weight = weight/weight_norm;
		weights[p_i] = weight/weight_norm;
	}

}

void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight. 
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution

	// use discrete distribution to perform re-sampling
	std::vector<Particle> new_particles;
	default_random_engine gen;
	std::discrete_distribution<> d(weights.begin(), weights.end());

	for (int n=0; n < num_particles; n++) {
		new_particles.push_back(particles[d(gen)]);
	}

	particles.clear();
	particles = new_particles;
}

Particle ParticleFilter::SetAssociations(Particle& particle, const std::vector<int>& associations, 
                                     const std::vector<double>& sense_x, const std::vector<double>& sense_y)
{
    //particle: the particle to assign each listed association, and association's (x,y) world coordinates mapping to
    // associations: The landmark id that goes along with each listed association
    // sense_x: the associations x mapping already converted to world coordinates
    // sense_y: the associations y mapping already converted to world coordinates

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
