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
#include <sstream>
#include <string>
#include <iterator>

#include "particle_filter.h"

using namespace std;




normal_distribution<double> distribution(0,1);

void ParticleFilter::init(double x, double y, double theta, double std[]) {
	// TODO: Set the number of particles. Initialize all particles to first position (based on estimates of 
	//   x, y, theta and their uncertainties from GPS) and all weights to 1. 
	// Add random Gaussian noise to each particle.
	// NOTE: Consult particle_filter.h for more information about this method (and others in this file).

	 num_particles = 200;
	 weights.resize(num_particles);

    random_device rd_engine;
    mt19937 generator(rd_engine());

    //initializing randomly all the particles

	 for (int i = 0; i < num_particles; i++) {
		 		Particle p_i;
				p_i.id = i;
				p_i.x = x + std[0]*distribution(generator);
				p_i.y = y + std[1]*distribution(generator);
				p_i.theta = theta + std[2]*distribution(generator);
                p_i.weight = 1.0f;
                weights[i] = 1.0f;
				particles.push_back(p_i);

			}
      is_initialized = true;


}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/
        random_device rd_engine;
        mt19937 generator(rd_engine());

     for (int i = 0; i < num_particles; i++) {

             double th = particles[i].theta;
             double v = velocity;
             double dth = yaw_rate;

            //checking to make sure yaw not zero

            cout<<"dth is"<<dth<<endl;
            if (fabs(dth)>0.001){
                    particles[i].x +=   v/dth*(sin(th+dth*delta_t) - sin(th))+ std_pos[0]*distribution(generator);
                    particles[i].y +=  -v/dth*(cos(th+dth*delta_t) - cos(th))+ std_pos[1]*distribution(generator);
                    particles[i].theta += dth*delta_t+ std_pos[2]*distribution(generator);
                } else {
                    particles[i].x += v*cos(th)*delta_t+ std_pos[0]*distribution(generator);
                    particles[i].y += v*sin(th)*delta_t+ std_pos[1]*distribution(generator);
                    particles[i].theta += dth*delta_t+ std_pos[2]*distribution(generator);
            }
            cout<<"particle  is"<<particles[i].x<<" " <<particles[i].y  <<particles[i].theta<<endl;


	}
}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the 
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to 
	//   implement this method and use it as a helper during the updateWeights phase.

    double dist_1 = 0;
	for (int i=0;i<observations.size();i++){ // loop over observations
		double min_dist = 100000;
		int closest_id = -1;
		for (int j=0;j<predicted.size();j++){

			dist_1 = dist(observations[i].x, observations[i].y,
										predicted[j].x, predicted[j].y);
			if (dist_1<min_dist){
				min_dist = dist_1;
				closest_id = predicted[j].id;
			}
		}
		//putting nearest observation
		observations[i].id  = closest_id;
	}
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
    weights.clear();


	for (int i_p=0;i_p<particles.size();i_p++){

        vector<LandmarkObs> observations_1;
        for (int i=0;i<observations.size();i++){
          LandmarkObs obsv_i;
          obsv_i.x = 0;
                obsv_i.x += observations[i].x * cos(particles[i_p].theta);
                obsv_i.x += -observations[i].y * sin(particles[i_p].theta);
                obsv_i.x += particles[i_p].x;


          obsv_i.y = 0;
          obsv_i.y += observations[i].x * sin(particles[i_p].theta);
                obsv_i.y += observations[i].y * cos(particles[i_p].theta);
                obsv_i.y += particles[i_p].y;

                obsv_i.id = -1; // Temporary ID.
          observations_1.push_back(obsv_i);
                }



        std::vector<LandmarkObs> predicted_meas;
            // Compute predicted measurements
            for (int i=0;i<map_landmarks.landmark_list.size();i++){
                double distance_particle_obs;
                distance_particle_obs = dist(particles[i_p].x,particles[i_p].y,
                                                                        map_landmarks.landmark_list[i].x_f,
                                                                        map_landmarks.landmark_list[i].y_f);
                if (distance_particle_obs <= sensor_range){
                    LandmarkObs predicted_i;
                    predicted_i.id = map_landmarks.landmark_list[i].id_i;
                    predicted_i.x = map_landmarks.landmark_list[i].x_f;
                    predicted_i.y = map_landmarks.landmark_list[i].y_f;

                    predicted_meas.push_back(predicted_i);
                }
             }





             dataAssociation(predicted_meas, observations_1);


             double prob = 1.0;
         double prob_i;



             for (int i =0; i<predicted_meas.size(); i++){
                 int ind_min = -1;
                 double dist_min = 1000000;


                 for (int j =0; j<observations_1.size(); j++){

                     if (predicted_meas[i].id == observations_1[j].id ){
                         double check_dist = dist(predicted_meas[i].x,
                                      predicted_meas[i].y,
                                                                        observations_1[j].x,
                                      observations_1[j].y);
                        if (check_dist<dist_min){
                            ind_min = j;
                            dist_min = check_dist;
                        }
                     }
                 }
                 if (ind_min!=-1){

                    prob_i =    exp(-((predicted_meas[i].x-observations_1[ind_min].x)*(predicted_meas[i].x-observations_1[ind_min].x)/(2*std_landmark[0]*std_landmark[0])
                            + (predicted_meas[i].y-observations_1[ind_min].y)*(predicted_meas[i].y-observations_1[ind_min].y)/(2*std_landmark[1]*std_landmark[1]))) / (2.0*3.14159*std_landmark[0]*std_landmark[1]);

                    prob = prob*prob_i;
                 }


             }
             weights.push_back(prob);
             particles[i_p].weight = prob;
	}
}

void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight. 
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution

        random_device rd_wts;
	    mt19937 generator_wts(rd_wts());

    // Creates a discrete distribution for weight.
    discrete_distribution<int> distribution_wts(weights.begin(), weights.end());
    vector<Particle> reSampled;

    for(int i=0;i<num_particles;i++){
				Particle particles_i = particles[distribution_wts(generator_wts)];
        reSampled.push_back(particles_i);
    }
    particles = reSampled;
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
