/*
 * particle_filter.cpp
 *
 *  Created on: Dec 12, 2016
 *      Author: Tiffany Huang
 */

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


	num_particles = 50;

	double std_x, std_y, std_theta; // Standard deviations for x, y, and theta

	// Set standard deviations for x, y, and theta.
	 std_x = std[0];
	 std_y = std[1];
	 std_theta = std[2];


	// This line creates a normal (Gaussian) distribution for x.
	double gps_x = x;
	normal_distribution<double> dist_x(gps_x, std_x);

	// Create normal distributions for y and theta.
	double gps_y = y;
	normal_distribution<double> dist_y(gps_y, std_x);
	normal_distribution<double> dist_theta(theta, std_theta);


	for (int i = 0; i < num_particles; ++i) {
		double sample_x, sample_y, sample_theta;

		// Sample from normal distrubtions
		// where "gen" is the random engine initialized in the header.
        sample_x = dist_x(gen_);
        sample_y = dist_y(gen_);
        sample_theta = dist_theta(gen_);

				Particle p;
				p.id = i;
				p.x = sample_x;
				p.y = sample_y;
				p.theta = sample_theta;
				p.weight = 1.0;

				//weights.push_back(1);
				particles.push_back(p);

		// Print your samples to the terminal.
		//cout << "Sample " << i + 1 << " " << sample_x << " " << sample_y << " " << sample_theta << endl;
	}

	is_initialized = true;

}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/

	// define normal distributions for sensor noise
	normal_distribution<double> N_x(0, std_pos[0]);
	normal_distribution<double> N_y(0, std_pos[1]);
	normal_distribution<double> N_theta(0, std_pos[2]);

	Particle p;
	double x, y, theta;
	for (int i = 0; i < num_particles; ++i) {
		p = particles[i];
		x = p.x;
		y = p.y;
		theta = p.theta;

		//avoid division by zero
    if (fabs(yaw_rate) > 0.0001) {
        particles[i].x = x + velocity/yaw_rate * (sin (theta + yaw_rate*delta_t) - sin(theta));
        particles[i].y = y + velocity/yaw_rate * (cos(theta) - cos(theta+yaw_rate*delta_t));
				particles[i].theta = theta + yaw_rate * delta_t;
    }
    else {
        particles[i].x = x + velocity * delta_t*cos(theta);
        particles[i].y = y + velocity * delta_t*sin(theta);
    }

		// add noise with zero mean
		particles[i].x += N_x(gen_);
		particles[i].y += N_y(gen_);
		particles[i].theta += N_theta(gen_);
	}
}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to
	//   implement this method and use it as a helper during the updateWeights phase.


    double min_dist = numeric_limits<double>::max();

	for (int o = 0; o < observations.size(); ++o) {
		// grab current observation
		x_obs = observations[o].x;
		y_obs = observations[o].y;

		// init minimum distance to maximum possible
		double min_dist = std::numeric_limits<double>::max();
		double cur_dist = 0.0;

		for (int p = 0; p < predicted.size(); ++p) {
			x_pred = predicted[p].x;
			y_pred = predicted[p].y;

			// Calculate Euclidean distance between the predicted and observed measurement
			cur_dist = dist(x_obs,y_obs,x_pred,y_pred);
			if (min_dist > cur_dist) {
				min_dist = cur_dist;
				// associate landmark from map with the observation
				observations[o].id = predicted[p].id;
			}
		}
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

	// foreach particle ...
	for (int i = 0; i < num_particles; ++i) {
		// get its x and y map coordinates
		p_x = particles[i].x;
		p_y = particles[i].y;
		p_theta = particles[i].theta;

		// create a vector to hold the map landmark locations predicted to be within sensor range of the particle
    std::vector<LandmarkObs> predictions;

    // for each map landmark...
    for (int l = 0; j < map_landmarks.landmark_list.size(); l++) {
			// get id and x,y coordinates
			lm_x = map_landmarks.landmark_list[l].x;
			lm_y = map_landmarks.landmark_list[l].y;
			lm_id = map_landmarks.landmark_list[l].id;

			 // only consider landmarks within sensor range of the particle
			if (dist(p_x, p_y, lm_x, lm_y) < sensor_range) {
				// add prediction to vector
        predictions.push_back(LandmarkObs{lm_id, lm_x, lm_y});
			}
		}
		// create and populate a copy of the list of observations transformed from vehicle/particle coordinates to map coordinates
	    vector<LandmarkObs> transformed_os;
	    for (int j = 0; j < observations.size(); j++) {
	      double t_x = cos(p_theta)*observations[j].x - sin(p_theta)*observations[j].y + p_x;
	      double t_y = sin(p_theta)*observations[j].x + cos(p_theta)*observations[j].y + p_y;
	      transformed_os.push_back(LandmarkObs{ observations[j].id, t_x, t_y });
	    }

	    // perform dataAssociation for the predictions (landmarks within particle range) and transformed observations on current particle
			// predictions are landmarks within sensor range (of the current particle) and therefore given in map coordinates
			// transformed observations are the sensor measurements of the vehicle but transformed from particle to map coordinates
	    dataAssociation(predictions, transformed_os);

			// reinit weight
	    particles[i].weight = 1.0;

			// foreach transformed observation ...
	    for (int j = 0; j < transformed_os.size(); j++) {

	      // placeholders for observation and associated prediction coordinates
	      double o_x, o_y, pr_x, pr_y;
	      o_x = transformed_os[j].x;
	      o_y = transformed_os[j].y;

	      int associated_prediction = transformed_os[j].id;

	      // get the x,y coordinates of the prediction associated with the current observation
	      for (int k = 0; k < predictions.size(); k++) {
	        if (predictions[k].id == associated_prediction) {
	          pr_x = predictions[k].x;
	          pr_y = predictions[k].y;
	        }
	      }

	      // calculate weight for this observation with multivariate Gaussian
	      double s_x = std_landmark[0];
	      double s_y = std_landmark[1];
	      double obs_w = (1/(2*M_PI*s_x*s_y)) * exp( -( pow(pr_x-o_x,2)/(2*pow(s_x, 2)) + (pow(pr_y-o_y,2)/(2*pow(s_y, 2))) ) );

	      // product of this obersvation weight with total observations weight
	      particles[i].weight *= obs_w;
	    }
	  }
}

void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight.
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution

	vector<Particle> resampled_particles;

  // get all of the current weights
  vector<double> weights;
  for (int i = 0; i < num_particles; ++i) {
    weights.push_back(particles[i].weight);
  }

	std::random_device rd;
  std::mt19937 gen(rd());
  std::discrete_distribution<> d(weights.begin(), weights.end());

  for(int i = 0; i < num_particles; ++i) {
      resampled_particles.push_back(particles[d(gen)]);
  }

	particles = resampled_particles;

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
