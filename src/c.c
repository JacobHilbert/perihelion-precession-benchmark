#include <stdio.h>
#include <math.h>

static const double pi = 3.141592653589793;
static const double mu = 39.4234021;       // gravitational parameter G*M_sun. AU^3/Years^2
static const double a = 0.38709893;        // semimajor axis
static const double alpha = 1.09778201e-8; // units of AU^2, is 3L^2/(mc)^2, where m is the reduced mass and L the angular momentum
static const double e = 0.20563069;        // eccentricity

static const double dt_max = 1e-6;
static const double dt_min = 1e-12;
static const double tmax = 10.0;

double rad_to_arcsec(const double radians) {
	return radians*180/pi*60*60;
}

struct Point { double x,y; };
typedef struct Point point;
double point_norm(const point p);
double point_angle(const point p);
point point_add(const point p1, const point p2);
point point_scalar_mul(const point p, const double scalar);
double point_dot_product(const point p1, const point p2);
point point_midpoint(const point p1, const point p2);

point acceleration(const point r_vec);


int main(){
	point r_now;
	point r_past;
	point v_now;
	point v_past;
	point a_now;
	point a_past;
	double event_now;
	double event_past;
	double dt = dt_max;
	double t;
	r_past = (point){a*(1+e),0.};
	v_past = (point){0.,sqrt((mu/a)*(1-e)/(1+e))};
	a_past = acceleration(r_past);
	event_past = point_dot_product(r_past,v_past);
	
	double angle;
	t = 0.0;
	while (t<tmax){
		r_now = point_add(point_add(r_past,point_scalar_mul(v_past,dt)),point_scalar_mul(a_past,0.5*dt*dt));
		a_now = acceleration(r_now);
		v_now = point_add(v_past,point_scalar_mul(point_midpoint(a_now,a_past),dt));
		
		event_now = point_dot_product(r_now,v_now);
		
		if (event_past<0 && event_now>0) {
			if (dt < dt_min){
				angle = point_angle(point_midpoint(r_now,r_past));
				printf("%.8f\t%.8f\n",t,rad_to_arcsec(angle-pi));
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

/*******************************/

double point_norm(const point p){
	return hypot(p.x,p.y);
}

// returns the angle in radians, adjusted to be in [0,2pi)
double point_angle(const point p){
	double result = atan2(p.y,p.x);
	if (result < 0)
		result += 2*pi;
	return result;
}

point point_add(const point p1, const point p2){
	point result;
	result.x = p1.x + p2.x;
	result.y = p1.y + p2.y;
	return result;
}
point point_scalar_mul(const point p, const double scalar){
	point result;
	result.x = p.x * scalar;
	result.y = p.y * scalar;
	return result;
}
double point_dot_product(const point p1, const point p2){
	return p1.x*p2.x + p1.y*p2.y;
}
point point_midpoint(const point p1, const point p2){
	return point_scalar_mul(point_add(p1,p2),0.5);
}

point acceleration(const point r_vec){
	const double r = point_norm(r_vec);
	const double beta = -mu/(r*r*r)*(1+alpha/(r*r));
	return point_scalar_mul(r_vec,beta);
}