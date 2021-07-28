#include <iostream>
#include <cmath>

const double pi = 4.0*atan(1.);
const double mu = 39.4234021;       // gravitational parameter G*M_sun. AU^3/Years^2
const double a = 0.38709893;        // semimajor axis
const double alpha = 1.09778201e-8; // units of AU^2, 3L^2/(mc)^2, where m is the reduced mass and L the angular momentum
const double e = 0.20563069;        // eccentricity

const double tmax = 10.0;
const double dt_max = 1e-6;
const double dt_min = 1e-12;

double rad_to_arcsec(const double& radians){ return radians*180/pi * 60 * 60; }

struct Point {
	double x,y;
	Point(){}
	Point(double x0, double y0){ x = x0; y = y0; }
	double norm() const { return sqrt(x*x+y*y); }
	double angle() const {
		double temp = atan2(y,x);
		if (temp < 0){
			return temp+2*pi;
		} else {
			return temp;
		}
	}
};

const Point operator + (const Point& p1, const Point& p2){ return Point(p1.x+p2.x, p1.y+p2.y); }
const Point operator * (const Point& point, const double& scalar){ return Point(point.x*scalar, point.y*scalar); }
const Point operator * (const double& scalar, const Point& point){return point * scalar;}
const Point operator / (const Point& point, const double& scalar){return Point(point.x/scalar, point.y/scalar);}

double dot_product(const Point& a, const Point& b){ return a.x*b.x + a.y*b.y; }

Point midpoint(const Point p1, const Point p2){ return (p1+p2)*0.5; }

Point acceleration(const Point& r_vec){
	const double r = r_vec.norm();
	const double beta = -mu/(r*r*r)*(1+alpha/(r*r));
	return r_vec * beta;
}

Point r_now,r_past, v_now,v_past, a_now,a_past;
double event_now, event_past;
double dt,t;

int main(){
	r_past = Point(a*(1+e),0.);
	v_past = Point(0.,sqrt((mu/a)*(1-e)/(1+e)) );
	event_past = dot_product(r_past,v_past);
	
	dt = dt_max;
	t = 0.0;
	while (t<tmax) {
		r_now = r_past + v_past*dt + a_past*dt*dt/2.0;
		a_now = acceleration(r_now);
		v_now = v_past + dt*(a_now+a_past)/2.0;
		event_now = dot_product(r_now,v_now);
		
		if (event_past<0 and event_now>0) {
			if (dt < dt_min) {
				std::cout << t << " " << rad_to_arcsec( midpoint(r_now,r_past).angle()-pi ) << std::endl;
				dt = dt_max;
			} else {
				// rewind time
				t -= dt;
				r_now = r_past;
				v_now = v_past;
				a_now = a_past;
				event_now = event_past;
				// adjust dt
				dt /= 2.0;
			}
		}
		
		r_past = r_now;
		v_past = v_now;
		a_past = a_now;
		event_past = event_now;
		t += dt;
	}
	
	return 0;
}
