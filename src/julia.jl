using Base: Float64
const π = 3.1415926535
const μ = 39.4234021
const a = 0.38709893
const α = 1.09778201e-8
const e = 0.20563069

rad_to_arcsec(radians::Float64)::Float64 = radians*180/π*60*60

struct Point
	x::Float64
	y::Float64
	Point(x,y) = new(x,y)
	Point() = new(0.0,0.0)
end

norm(p::Point)::Float64 = hypot(p.x, p.y)
angle(p::Point)::Float64 = mod(atan(p.y, p.x),2π)
dot(p1::Point,p2::Point)::Float64 = p1.x*p2.x + p1.y*p2.y

import Base: +,*,/
(p1::Point + p2::Point)::Point = Point(p1.x+p2.x, p1.y+p2.y)
(p::Point * scalar)::Point = Point(p.x*scalar,p.y*scalar)
(scalar * p::Point)::Point = p*scalar
(p::Point / scalar)::Point = Point(p.x/scalar,p.y/scalar)
	
function acceleration(r_vec::Point)::Point
	r::Float64 = norm(r_vec)
	β::Float64 = -μ/r^3 * (1+α/r^2)
	return r_vec * β
end

function verlet(tmax::Float64,dt_max::Float64,dt_min::Float64)::Nothing
	r_past = Point(a*(1+e),0.)
	v_past = Point(0.,sqrt((μ/a)*(1-e)/(1+e)))
	a_past = acceleration(r_past)
	
	r_now = Point()
	v_now = Point()
	a_now = Point()
	
	event_past = dot(r_past,v_past)
	event_now = 0.0
	
	dt = dt_max
	t = 0.0
	while t < tmax
		r_now = r_past + v_past*dt + a_past*dt^2/2
		a_now = acceleration(r_now)
		v_now = v_past + dt*(a_now+a_past)/2
		event_now = dot(r_now,v_now)
		
		if event_past<0 && event_now>0
			if dt<dt_min
				println(t," ",rad2deg(angle((r_now+r_past)/2)-π)*3600)
				dt = dt_max
			else
				# rewind time
				t -= dt
				r_now = r_past
				v_now = v_past
				a_now = a_past
				event_now = event_past
				# adjust dt
				dt /= 2
			end
		end
		
		r_past = r_now
		v_past = v_now
		a_past = a_now
		event_past = event_now
		t += dt
	end
	return nothing
end

# compile (?)
verlet(1e-1,1e-2,1e-3)

println()

verlet(10.0,1e-6,1e-12)