#include <iostream>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <vector>
#include <fstream>
#include <numeric>
#include <algorithm>
#include <ctime>
#include <random>
#include <unordered_map>

using namespace std;

const float dt = 1e-4;
const int monitor = 5000;
const int predation_monitor = 10000;
float t_final = 2000.0;
float t_initial = 0.0;
//float satisfy_timesteps = 0.0;
//float time_BELT_starts = 0.0;

// Update in case of self-propelled particles
float time_SYS_starts = 0.0;
float time_HUNT_starts = 2500.0;
int hunt_timelimit = 0;		// In terms of iteration count
int satisfy_timelimit = 0;	// In terms of iteration count
int refocus_timelimit = 0;	// In terms of iteration count
/**********************************************/

float Cv = 0.0;
const float SPRING = 20000;
const float SPRING_WALL = 5000;
const float eta = 2;
const int n_ex = 800;
float exp_power;
float exf[n_ex];

//float gaussStdDevPercRAD = 0.05; // Uncomment this line to specify the standard deviation of the Gaussian w.r.t. mean radius (RAD)

const float SIGMA = 2.9256;
const float RHO = 1000; //1.22;

const float GRAV = 0; //9.81;

float RAD = 0.0136;
float Pred_RAD = RAD*4.0;
//float Secondary_RAD = 0.00125;
//float SMALL_RAD = 0.001;
//float BIG_RAD = 0.01;

//const float Wall_RAD = 0.0025;

const float X0 = 0.0;
const float Y0 = 0.0;
const float WIDTH = 2.5;
//const float HEIGHT = 0.4;
const float HEIGHT = WIDTH;
//const float xCentre = WIDTH/2.0;
//const float yCentre = HEIGHT/2.0;

// Zones for prey
float h_align = 5.0;	// r_align = h_align*d
float h_separate = 0.5;	// r_separate = h_separate*d
float h_attract = 20.0;	// r_attract = h_attract*d
float h_escape = 20.0;	// r_escape = h_escape*d

// Zones for predator
float h_detect = 5.0; 	// r_detect = h_detect*D
float r_kill = RAD + Pred_RAD;

float r_skin = 0.5*RAD;
float r_cutoff = h_attract*2.0*RAD;

float V_surr = M_PI*pow(r_cutoff,2);

float Xmin = X0;		//-2.0*Wall_RAD;
float Ymin = Y0;		//-2.0*Wall_RAD;
float Xmax = WIDTH;		// + 2.0*Wall_RAD;
float Ymax = HEIGHT;	// + 2.0*Wall_RAD;

//float CONC = 0.00;
const int CORE_PARTICLES = 508;

//int bot_Core_small = 0;
//int bot_Core_big = 0;
//int top_Core_small = 0;
//int top_Core_big = 0;
int m = ceil((Ymax - Ymin)/(r_cutoff+r_skin));
int n = ceil((Xmax - Xmin)/(r_cutoff+r_skin));
int ncells = m*n;

vector<vector<int> > nb_part(CORE_PARTICLES);

//float BELT_ini_vel = 0.0;
//float BELT_vel_x_bot = 0.0;
//float BELT_vel_y_right = 0.0;
//float BELT_vel_x_top = 0.0;
//float BELT_vel_y_left = 0.0;

vector<vector<int>> hunt_list;
vector<vector<int>> predating_list;

unordered_map<int,char> escape_list;	// unordered map conisting of prey IDs with escape activated

float beta_hunt = 0.08;
float beta_escape = 15.00;
float beta_repulse = 0.25;

float fric_factor = 0.1;
float noise = 0.0;
vector<float> noise_val(CORE_PARTICLES);
float blind_angle = 0.0;
float F_atr_max = 1e-6;	// Maximum attraction force for danios = 1.2e-3
float wall_avoid_dist = 8*RAD;
float pred_avoid_dist = 0.05*(WIDTH+HEIGHT)/2;
float rot_angle_avoid = M_PI/4;

float vel_limit = 6.0*2.0*RAD;
float vmag_init = 6.0*2.0*RAD;
vector<float> vel_init(CORE_PARTICLES);
//int obs_p_id = 1 + rand() % CORE_PARTICLES;

//ofstream ofp_obs("obs_data.txt", ios::out);

//int satisfy_list[CORE_PARTICLES];

inline float rand01() {return static_cast<double> (rand())/RAND_MAX;}
inline double rand02(double range)
{
	double c = 1 + (rand()/(RAND_MAX+1.0))*(1+range);
	return c;
}
void gaussian_white(float mag) // returns Gaussian white random number between min and max
{
    int seed = clock();
    mt19937 gen(seed); //mersenne twister
    for (unsigned short i = 0; i < noise_val.size(); i++)
    {
        normal_distribution <float> dist(0, sqrt(mag*mag/12.0));
        noise_val[i] = dist(gen);
    }
}
void uniform_random(float lower_lim, float upper_lim) // returns uniform random number between the limits
{
    int seed = clock();
    mt19937 gen(seed); //mersenne twister
    for (unsigned short i = 0; i < vel_init.size(); i++)
    {
        uniform_real_distribution <float> dist(lower_lim, upper_lim);
        vel_init[i] = dist(gen);
    }
}
int signum(auto x)
{
    return (x > 0) ? 1 : ((x < 0) ? -1 : 0); 
}
/*inline float randGauss(float sigma)
{
    float phi = rand01()*2*M_PI;
    float Upsilon = rand01();
    float Psi = -log(Upsilon);
    float r = sigma*sqrt(2.0*Psi);
    float x = r*cos(phi);
    float y = r*sin(phi);
    if (rand01() > 0.5)
        return x;
    else
        return y;
}

// for a more accurate calculation of neighbouring particles' volume within the cutoff radius
float circle_int_area(float r1, float r2, float d)
{
    float r = r2;
    float R = r1;
    if(R < r)
    {
        // swap
        r = r1;
        R = r2;
    }
    if (d < R - r)
        return M_PI*r*r;
    else
    {
        float a1 = r*r*acos((d*d + r*r - R*R)/(2*d*r));
        float a2 = R*R*acos((d*d + R*R - r*r)/(2*d*R));
        float a3 = 0.5*sqrt((-d+r+R)*(d+r-R)*(d-r+R)*(d+r+R));
        return a1 + a2 - a3;
    }
}*/

//ofstream ofp_energy("energy.txt", ios::out);
//ofstream ofp_monitor("monitor.txt", ios::out);
//float wall_energy_input = 0.0;
//float sp_pot_energy = 0.0;

class vect
{
public:
    float x, y;
    vect()      //default constructor
        {
        x = y = 0.0f;
        }
    vect(float tx, float ty) //passing x and y position values at time t
        {
        x = tx;
        y = ty;
        }
    ~vect() {}      //destructor
    //Verlet-velocity algorithm start
    vect distance_calc(vect & d_other) //address of d_other as argument
    {
        /* Uncomment below code for confined domain
         * vect d;
         * d.x = d_other.x - x;
         * d.y = d_other.y - y;
         * return sqrt(d.x * d.x + d.y * d.y);
         */
        /**** Added distance calculation algorithm for periodic domain. (Updated on 08/08/20) ****/
        vect d1;
        float hx = WIDTH/2.0, hy = HEIGHT/2.0;
        d1.x = d_other.x - x;
        d1.y = d_other.y - y;
        
        if (d1.x > hx)	d1.x -= WIDTH;
        else if (d1.x < -hx) d1.x += WIDTH;
        
        if (d1.y > hy)	d1.y -= HEIGHT;
        else if (d1.y < -hy)	d1.y += HEIGHT;
        
        //if (*dz>hL)       *dz-=L;
        //else if (*dz<-hL) *dz+=L;
        return d1;
        //	Returning the vector difference between the particle positions
    }
    float dotProduct(vect & other)	// Returns the dot product between two vectors
    {
    	return (x*other.x + y*other.y);
    }
    void add_pos_vel(const vect & v, const vect & f, const float m, float dt = 0.0f)
        {
        if (dt == 0.0f)
            {
            x += v.getX();
            y += v.getY();
            }
        else
            {
            x += v.getX() * dt + 0.5*(f.getX()/m)*dt*dt;
            y += v.getY() * dt + 0.5*(f.getY()/m)*dt*dt;
            }
        }
    void add_pos_vel(const vect & f, const float m, float dt = 0.0f)
        {
        if (dt == 0.0f)
            {
            x += f.getX();
            y += f.getY();
            }
        else
            {
            x += 0.5*(f.getX()/m) * dt;
            y += 0.5*(f.getY()/m) * dt;
            }
        }
    //Verlet-velocity algorithm end
/*    void addc(float f)
        {
        x += f;
        y += f;
        }
    void subtract(const vect & v)
        {
        x -= v.getX();
        y -= v.getY();
        }
    void scale(float mul)
        {
        x *= mul;
        y *= mul;
        }
    void normalize()
        {
        float ln = magnitude();
        x /= ln;
        y /= ln; // the present particle is in the sphere of influence of the 'other'
        }*/
    float magnitude() const {return sqrt( (x)*(x) + (y)*(y) );}

    float getX() const {return x;}
    float getY() const {return y;}

    void setX(float v) {x = v;}
    void setY(float v) {y = v;}
    
    void rotate_vect(bool dirn)	// true for CCW and false for CW
    {
    	float x_tmp = x;
    	if (dirn)
    	{
    		x = x * cos(rot_angle_avoid) - y * sin(rot_angle_avoid);
    		y = x_tmp * sin(rot_angle_avoid) + y * cos(rot_angle_avoid);
    	}
    	else
    	{
    		x = x * cos(-rot_angle_avoid) - y * sin(-rot_angle_avoid);
    		y = x_tmp * sin(-rot_angle_avoid) + y * cos(-rot_angle_avoid);
		}
	}
};

class particle : public vect
    {
    public:
        vect position, position_prev;
        float displace;
        vect velocity, vf;//, velocity_prev;
        vect force;//, force_prev;
        float radius, mass, beta;
//        int species;
//        unsigned short wall_tag;
//        unsigned short wall_layer;
        float vf_0, v_surr_part, V_surr;
//        vect local_order_param;
//        float mass_surr = 0, mass_bot = 0, mass_top = 0;
//        short unsigned num_surr = 0, num_bot = 0, num_top = 0;
        int id, wall_tag; //, wall_layer, kind;
        
        particle()
            {
            radius = RAD; //creates disks of radius=RAD.
            position.x = rand01()*WIDTH;
            position.y = rand01()*HEIGHT;
            position_prev.x = 0.0;
            position_prev.y = 0.0;
            displace = 0.0;
            mass = SIGMA*M_PI*radius*radius;
            beta = 0.0;
            velocity.x = 0.0; //initial velocity and force is zero.
            velocity.y = 0.0;
//            velocity_prev.x = 0.0; //initial velocity and force is zero.
//            velocity_prev.y = 0.0;
            force.x = 0.0;
            force.y = 0.0;
//            force_prev.x = 0.0;
//            force_prev.y = 0.0;
            v_surr_part = 0;
            vf_0 = 0.0;
            V_surr = M_PI*pow(h_attract*2.0*radius,2);
//            species = 0;
            wall_tag = 0;
/*            wall_layer = 0;
            local_order_param.x = 0.0;
            local_order_param.y = 0.0;
            mass_surr = 0.0; num_surr = 0;
            mass_bot = 0.0; num_bot = 0;
            mass_top = 0.0; num_top = 0;            
			kind = 0;
*/			
            id = 0;
            }

        particle(float u)
            {
            radius = RAD;
            position.x = rand01()*WIDTH;
            position.y = rand01()*u;
            position_prev.x = 0.0;
            position_prev.y = 0.0;
            displace = 0.0;
            mass = SIGMA*M_PI*radius*radius;
            beta = 0.0;
            velocity.x = 0.0;
            velocity.y = 0.0;
//            velocity_prev.x = 0.0; //initial velocity and force is zero.
//            velocity_prev.y = 0.0;
            force.x = 0.0;
            force.y = 0.0;
//            force_prev.x = 0.0;
//            force_prev.y = 0.0;
            v_surr_part = 0;
            vf_0 = 0.0;
            V_surr = M_PI*pow(h_attract*2.0*radius,2);
//            species = 0;
            wall_tag = 0;
/*            wall_layer = 0;
            local_order_param.x = 0.0;
            local_order_param.y = 0.0;
            mass_surr = 0.0; num_surr = 0;
            mass_bot = 0.0; num_bot = 0;
            mass_top = 0.0; num_top = 0;            
			kind = 0;
*/
            id = 0;
        }

        particle(float u,float v)
            {
            radius = RAD;
            position.x = u;
            position.y = v;
            position_prev.x = u;
            position_prev.y = v;
            displace = 0.0;
            mass = SIGMA*M_PI*radius*radius;
            beta = 0.0;
            velocity.x = 0.0;
            velocity.y = 0.0;
//            velocity_prev.x = 0.0; //initial velocity and force is zero.
//            velocity_prev.y = 0.0;
            force.x = 0.0;
            force.y = 0.0;
//            force_prev.x = 0.0;
//            force_prev.y = 0.0;
            v_surr_part = 0;
            vf_0 = 0.0;
            V_surr = M_PI*pow(h_attract*2.0*radius,2);
//            species = 0;
            wall_tag = 0;
/*            wall_layer = 0;
            local_order_param.x = 0.0;
            local_order_param.y = 0.0;
            mass_surr = 0.0; num_surr = 0;
            mass_bot = 0.0; num_bot = 0;
            mass_top = 0.0; num_top = 0;            
			kind = 0;
*/			
            id = 0;
        }

        particle(float u,float v, float r)
            {
            radius = r;
            position.x = u;
            position.y = v;
            position_prev.x = u;
            position_prev.y = v;
            displace = 0.0;
            mass = SIGMA*M_PI*radius*radius;
            beta = 0.0;
            velocity.x = 0.0;
            velocity.y = 0.0;
//            velocity_prev.x = 0.0; //initial velocity and force is zero.
//            velocity_prev.y = 0.0;
            force.x = 0.0;
            force.y = 0.0;
//            force_prev.x = 0.0;
//            force_prev.y = 0.0;
            v_surr_part = 0;
            vf_0 = 0.0;
            V_surr = M_PI*pow(h_attract*2.0*radius,2);
//            species = 0;
            wall_tag = 0;
/*            wall_layer = 0;
            local_order_param.x = 0.0;
            local_order_param.y = 0.0;
            mass_surr = 0.0; num_surr = 0;
            mass_bot = 0.0; num_bot = 0;
            mass_top = 0.0; num_top = 0;            
			kind = 0;
*/
            id = 0;
        }
        
        particle(float u,float v, float r, float vel_x, float vel_y)
            {
            radius = r;
            position.x = u;
            position.y = v;
            position_prev.x = u;
            position_prev.y = v;
            displace = 0.0;
            mass = SIGMA*M_PI*radius*radius;
            beta = 0.0;
            velocity.x = vel_x;
            velocity.y = vel_y;
//            velocity_prev.x = 0.0; //initial velocity and force is zero.
//            velocity_prev.y = 0.0;
            force.x = 0.0;
            force.y = 0.0;
//            force_prev.x = 0.0;
//            force_prev.y = 0.0;
            v_surr_part = 0;
            vf_0 = 0.0;
            V_surr = M_PI*pow(h_attract*2.0*radius,2);
//            species = 0;
            wall_tag = 0;
/*            wall_layer = 0;
            local_order_param.x = 0.0;
            local_order_param.y = 0.0;
            mass_surr = 0.0; num_surr = 0;
            mass_bot = 0.0; num_bot = 0;
            mass_top = 0.0; num_top = 0;            
			kind = 0;
*/
            id = 0;
        }

        particle(float u,float v, float r, float vel_x, float vel_y, int part_type)
        {
            radius = r;
            position.x = u;
            position.y = v;
            position_prev.x = u;
            position_prev.y = v;
            displace = 0.0;
            mass = SIGMA*M_PI*radius*radius;
            beta = 0.0;
            velocity.x = vel_x;
            velocity.y = vel_y;
//            velocity_prev.x = 0.0; //initial velocity and force is zero.
//            velocity_prev.y = 0.0;
            force.x = 0.0;
            force.y = 0.0;
//            force_prev.x = 0.0;
//            force_prev.y = 0.0;
            v_surr_part = 0;
            vf_0 = 0.0;
            V_surr = M_PI*pow(h_attract*2.0*radius,2);
//            species = 0;
            wall_tag = part_type;
/*            wall_layer = 0;
            local_order_param.x = 0.0;
            local_order_param.y = 0.0;
            mass_surr = 0.0; num_surr = 0;
            mass_bot = 0.0; num_bot = 0;
            mass_top = 0.0; num_top = 0;            
			kind = part_type;
*/
            id = 0;
        }

        ~particle() {}

        void setPrevPosx(float x_prev) { position_prev.x = x_prev; }
        void setPrevPosy(float y_prev) { position_prev.y = y_prev; }
        void setDisplace(float d) { displace = d; }
        void setVelx(float velx) { velocity.x = velx; }
        void setVely(float vely) { velocity.y = vely; }
        void addCx(float c) { position.x += c; }
        void addCy(float c) { position.y += c; }
        void setRadius(float r) { radius = r; }
        void setMass(float m) { mass = m; }
        void setBeta(float b) { beta = b; }
        void setVsurr(float vsurr) { V_surr = vsurr; }
//        void setSpecies(int tag) { species = tag; }
        void setWall(unsigned short tag) //, unsigned short layer)
        {
            wall_tag = tag;
//            wall_layer = layer;
        }
        bool vision_check(particle * other)
        {
        	// Calculating the angle between the velocity and the pos diff vectors using dot product
			// acos returns in range of [0,pi]
			
        	vect pos_diff;
        	float vmag, dot_prod, angle_diff;
        	
        	pos_diff = position.distance_calc(other->position);	// \vec_{pos2} - \vec_{pos1}
        	dot_prod = pos_diff.dotProduct(velocity);
        	
        	vmag = velocity.magnitude();
        	if (vmag == 0.0)	vmag = 1.0;
        	
        	angle_diff = acos(dot_prod/(pos_diff.magnitude() * vmag));
        	
        	if (angle_diff < M_PI - blind_angle)	return true;
        	else return false;	
        }
        float angle_diff_calc(particle * other)
        {
        	// Calculating the angle between the velocity and the pos diff vectors using dot product
			// acos returns in range of [0,pi]
			
        	vect pos_diff;
        	float vmag, dot_prod, angle_diff;
        	
        	pos_diff = position.distance_calc(other->position);	// \vec_{pos2} - \vec_{pos1}
        	dot_prod = pos_diff.dotProduct(velocity);
        	
        	vmag = velocity.magnitude();
        	if (vmag == 0.0)	vmag = 1.0;
        	
        	angle_diff = acos(dot_prod/(pos_diff.magnitude() * vmag));
        	
        	return angle_diff;	
        }
		void weight_other(particle * other)		// alignment
		{
			float tmpexf, masstmp, r_align, r_sep;
			int tmp_id;
            float distance = position.distance_calc(other->position).magnitude();
            
	    	// the 'other' particle is in the sphere of influence of the present
	    	r_align = h_align*2.0*radius;
	    	r_sep = h_separate*2.0*(radius + other->radius);
			if (distance > r_sep && distance < r_align && this->vision_check(other))
			{	 
			   	tmp_id = 1+n_ex*(distance/r_align);
			    tmpexf = exf[tmp_id];
				masstmp = other->mass*tmpexf;	
				v_surr_part += M_PI*pow(other->radius,2);
			    vf.x += masstmp * other->velocity.x; //total x-weighted velocity calculated
			    vf.y += masstmp * other->velocity.y; //total y-weighted velocity calculated
			    vf_0 += masstmp; //total weights calculated
			}
			
		    // the present particle is in the sphere of influence of the 'other'
			r_align = h_align*2.0*other->radius;
				
			if (distance > r_sep && distance < r_align && other->vision_check(this)) 
			{
			    tmp_id = 1+n_ex*(distance/r_align);
			    tmpexf = exf[tmp_id];
				masstmp = mass*tmpexf;
				other->v_surr_part += M_PI*pow(radius,2);
				other->vf.x += masstmp * velocity.x;
				other->vf.y += masstmp * velocity.y;
				other->vf_0 += masstmp;   
			}
		}
		
        void separate_other(particle * other, bool VerletUpdate = false)     // repulsion.
        {   
			//other is a pointer to an object
            
			vect f, dist;  //dist and f both have x and y members
			dist = position.distance_calc(other->position);
            float distance = dist.magnitude(); //pythagoras distance
            float r_sep = h_separate * 2.0*(other->radius + radius);            

            // Short range repulsion force function
            if (distance < r_sep)
            {
            	if ((wall_tag > 0 && other->wall_tag <= 0) || (wall_tag <= 0 && other->wall_tag > 0))
				//if ((kind != -1 && other->kind == -1) || (kind == -1 && other->kind != -1))
				{
					f.x = -(r_sep/distance -1) * SPRING_WALL * dist.x;  //Spring force in case either is a wall particle
					f.y = -(r_sep/distance -1) * SPRING_WALL * dist.y;
//					if (VerletUpdate)
//						sp_pot_energy += 0.5*SPRING_WALL*pow(minDist - distance,2);		
				}
				else if (wall_tag <= 0 && other->wall_tag <= 0)
				{             
	        	    f.x = -(r_sep/distance -1) * SPRING * dist.x;
	                f.y = -(r_sep/distance -1) * SPRING * dist.y;
	//              if (VerletUpdate)
	//              	sp_pot_energy += 0.5*SPRING*pow(minDist - distance,2); //0.5kx^2
				}
        	}
            //Force transfer
            force.x += f.x;
            force.y += f.y;

            other->force.x -= f.x; //Twice acting
            other->force.y -= f.y;
        }
/*        void orderCalcOther(particle * other) //Calculating the Order parameter, not required
        {
            vect a, dist;
            float angle;
            dist.x = other->position.x - position.x;
            dist.y = other->position.y - position.y;
            float distance = position.distance_calc(other->position).magnitude();
            float minDist = (other->radius + radius);

            if (distance < minDist*1.2)
            {
                angle = atan2(dist.y, dist.x);

                local_order_param.x += cos(6*angle);
                local_order_param.y += sin(6*angle);

                angle = atan2(-1.0*dist.y, -1.0*dist.x);
                other->local_order_param.x += cos(6*angle);
                other->local_order_param.y += sin(6*angle);

                num_surr +=1;
                other->num_surr += 1;
                mass_surr += other->mass;
                other->mass_surr += mass;

                if (other->species == 1)
                {
                    num_bot += 1;
                    mass_bot += other->mass;
                }
                else
                {
                    num_top += 1;
                    mass_top += other->mass;
                }

                if (species == 1)
                {
                    other->num_bot += 1;
                    other->mass_bot += mass;
                }
                else
                {
                    other->num_top += 1;
                    other->mass_top += mass;
                }
            }
        }
*/
        void update_pos(const float dt)
        {	
            position.add_pos_vel(velocity,force,mass,dt);
            float part_x = position.getX();
            float part_y = position.getY();
		
		//	Uncomment below lines for a confined domain
        /*            if ( part_x < X0)  {position.setX(RAD*2.0);} //if particle moves to left of origin
            if ( part_x > X0+WIDTH)  {position.setX(WIDTH - 2.0*RAD);} //if particle moves to right of width
 
            if ( part_y < Y0)  {position.setY(RAD*2.0);}
            if ( part_y > Y0+HEIGHT)  {position.setY(HEIGHT - 2.0*RAD);}
         */
        /**** Below code is for periodic domain. Added on 08/08/20. ****/
        if ( part_x < X0)  {position.setX(WIDTH + part_x);} //if particle moves to left of origin
        if ( part_x > X0+WIDTH)  {position.setX(part_x - WIDTH);} //if particle moves to right of width
        
        if ( part_y < Y0)  {position.setY(HEIGHT + part_y);}
        if ( part_y > Y0+HEIGHT)  {position.setY(part_y - HEIGHT);}
        /**** Update for periodic domain ends here. ****/
        
	// Uncomment below for circular domain
/*            float theta = atan2(part_y - yCentre, part_x - xCentre);
            float r = sqrt(pow((yCentre - part_y),2) + pow((xCentre - part_x),2));
            
            if ( r > WIDTH/2.0) //if particle moves outside domain
            {
            	position.setX((WIDTH/2 - 2.0*RAD)*cos(theta) + xCentre);
            	position.setY((WIDTH/2 - 2.0*RAD)*sin(theta) + yCentre);
            }
*/
        }

/*        void update_wall_pos(const float dt)
        {
            position.add_pos_vel(velocity,1.0, dt);
            float part_x = position.getX();
            float part_y = position.getY();

            if (wall_tag == 1 || wall_tag == 3) // Top or bottom wall
            {
                if (wall_layer == 1 || wall_layer == 3) // 1st or 3rd layer
                {
                    if ( part_x > WIDTH + 2*Wall_RAD)  {position.setX(part_x - (WIDTH + 4*Wall_RAD));}
                    if ( part_x < -2*Wall_RAD)  {position.setX(part_x + (WIDTH + 4*Wall_RAD));}
                }
                else                                    // 2nd layer
                {
                    if ( part_x > WIDTH + Wall_RAD)  {position.setX(part_x - (WIDTH + 2*Wall_RAD));}
                    if ( part_x < -Wall_RAD)  {position.setX(part_x + (WIDTH + 2*Wall_RAD));}
                }
            }

            if (wall_tag == 2 || wall_tag == 4) // Side walls (left or right)
            {
                if (wall_layer == 1 || wall_layer == 3) // 1st or 3rd layer
                {
                    if ( part_y > HEIGHT + 2*Wall_RAD)  {position.setY(part_y - (HEIGHT + 4*Wall_RAD));}
                    if ( part_y < -2*Wall_RAD)  {position.setY(part_y + (HEIGHT + 2*Wall_RAD));}
                }
                else                                    // 2nd layer
                {
                    if ( part_y > HEIGHT + Wall_RAD)  {position.setY(part_y - (HEIGHT + 2*Wall_RAD));}
                    if ( part_y < -Wall_RAD)  {position.setY(part_y + (HEIGHT + 2*Wall_RAD));}
                }

            }

        }
*/        
        void update_vel(const float dt) //Updates velocity by 1/2 timestep each time it is called
        {
        	velocity.add_pos_vel(force, mass, dt);
        	
        	// Comment below if there is no limit on velocity
        	float theta = atan2(velocity.y, velocity.x);
        	if (velocity.magnitude() > vel_limit)
        	{
        		velocity.x = vel_limit * cos(theta);
        		velocity.y = vel_limit * sin(theta);
        	}
        }

    };

vector<particle> disks; //empty vector array of particle type

class simulate : public particle
    {
    public:
        int stop = clock();
        vect gravity;
        vector<int> wall_limits;
        float time = t_initial;
//        vect displace_max;
//        float power_sum = 0;

        simulate()
            {
            gravity.x = 0.0;
            gravity.y = -GRAV;
            }
        ~simulate() {}

 /*       void boundaryParticlesGenerate()
        {
            ofstream ofp_wall("wall_number.txt", ios::out);
            ofstream ofp_wall_part_stat_1("wall_data_stat_1.txt", ios::out);
            ofstream ofp_wall_part_stat_2("wall_data_stat_2.txt", ios::out);
            ofstream ofp_wall_part_stat_3("wall_data_stat_3.txt", ios::out);
            ofstream ofp_wall_part_stat_4("wall_data_stat_4.txt", ios::out);

            wall_limits.push_back(disks.size());
            //            / Bottom wall (belt)  /
            for( float i = -Wall_RAD; i <= WIDTH + Wall_RAD ; i+= RAD * 2 )    // Layer 1
                {
                particle disk = particle(i, -Wall_RAD, Wall_RAD, BELT_ini_vel, 0.0); // Moving wall
                disk.setWall(1,1);    // bottom wall tag = 1
                disks.push_back(disk);
                if (fabs(BELT_vel_x_bot) == 0.0)
                ofp_wall_part_stat_1 << disk.position.x << "\t" << disk.position.y << "\t"<< disk.radius << "\t"<< disk.velocity.x << "\t" << disk.velocity.y << endl;
                }
            for( float i = 0; i <= WIDTH ; i+= RAD * 2 )    // Layer 2
                {
                particle disk = particle(i, -Wall_RAD*(1+2*sin(M_PI/3.0)), Wall_RAD, BELT_ini_vel, 0.0); // Moving wall
                disk.setWall(1,2);    // bottom wall tag = 1
                disks.push_back(disk);
                if (fabs(BELT_vel_x_bot) == 0.0)
                ofp_wall_part_stat_1 << disk.position.x << "\t" << disk.position.y << "\t"<< disk.radius << "\t"<< disk.velocity.x << "\t" << disk.velocity.y << endl;
                }
            for( float i = -Wall_RAD; i <= WIDTH + Wall_RAD ; i+= RAD * 2 )     // Layer 3
                {
                particle disk = particle(i, -Wall_RAD*(1+4*sin(M_PI/3.0)), Wall_RAD, BELT_ini_vel, 0.0); // Moving wall
                disk.setWall(1,3);    // bottom wall tag = 1
                disks.push_back(disk);
                if (fabs(BELT_vel_x_bot) == 0.0)
                ofp_wall_part_stat_1 << disk.position.x << "\t" << disk.position.y << "\t"<< disk.radius << "\t"<< disk.velocity.x << "\t" << disk.velocity.y << endl;
                }
            wall_limits.push_back(disks.size());
            unsigned short bot_wall = wall_limits.back() - wall_limits.front();
            ofp_wall << bot_wall << endl;
            ofp_wall_part_stat_1.close();

            //            /***************** Right wall /
            for( float i = HEIGHT + Wall_RAD; i >= -Wall_RAD ; i-= Wall_RAD * 2)     // Layer 1
                {
                particle disk = particle(WIDTH + Wall_RAD, i, Wall_RAD);
                disk.setWall(2,1);    // right wall tag = 2
                disks.push_back(disk);
                if (fabs(BELT_vel_y_right) == 0.0)
                ofp_wall_part_stat_2 << disk.position.x << "\t" << disk.position.y << "\t"<< disk.radius << "\t"<< disk.velocity.x << "\t" << disk.velocity.y << endl;
                }
            for( float i = HEIGHT; i >= 0.0 ; i-= Wall_RAD * 2)     // Layer 2
                {
                particle disk = particle(WIDTH + Wall_RAD*(1+2*sin(M_PI/3.0)), i, Wall_RAD);
                disk.setWall(2,2);    // right wall tag = 2
                disks.push_back(disk);
                if (fabs(BELT_vel_y_right) == 0.0)
                ofp_wall_part_stat_2 << disk.position.x << "\t" << disk.position.y << "\t"<< disk.radius << "\t"<< disk.velocity.x << "\t" << disk.velocity.y << endl;
                }
            for( float i = HEIGHT + Wall_RAD; i >= -Wall_RAD ; i-= Wall_RAD * 2)     // Layer 3
                {
                particle disk = particle(WIDTH + Wall_RAD*(1+4*sin(M_PI/3.0)), i, Wall_RAD);
                disk.setWall(2,3);    // right wall tag = 2
                disks.push_back(disk);
                if (fabs(BELT_vel_y_right) == 0.0)
                ofp_wall_part_stat_2 << disk.position.x << "\t" << disk.position.y << "\t"<< disk.radius << "\t"<< disk.velocity.x << "\t" << disk.velocity.y << endl;
                }
            wall_limits.push_back(disks.size());
            unsigned short right_wall = wall_limits.back() - bot_wall;
            ofp_wall << right_wall << endl;
            ofp_wall_part_stat_2.close();


            //            /***************** Top wall /
            for( float i = WIDTH + Wall_RAD; i >= -Wall_RAD ; i-= Wall_RAD * 2)     // Layer 1
                {
                particle disk = particle(i, HEIGHT + Wall_RAD, Wall_RAD);
                disk.setWall(3,1);    // top wall tag = 3
                disks.push_back(disk);
                if (fabs(BELT_vel_x_top) == 0.0)
                ofp_wall_part_stat_3 << disk.position.x << "\t" << disk.position.y << "\t"<< disk.radius << "\t"<< disk.velocity.x << "\t" << disk.velocity.y << endl;
                }
            for( float i = WIDTH; i >= 0.0 ; i-= Wall_RAD * 2)      // Layer 2
                {
                particle disk = particle(i, HEIGHT + Wall_RAD*(1+2*sin(M_PI/3.0)), Wall_RAD);
                disk.setWall(3,2);    // top wall tag = 3
                disks.push_back(disk);
                if (fabs(BELT_vel_x_top) == 0.0)
                ofp_wall_part_stat_3 << disk.position.x << "\t" << disk.position.y << "\t"<< disk.radius << "\t"<< disk.velocity.x << "\t" << disk.velocity.y << endl;
                }
            for( float i = WIDTH + Wall_RAD; i >= -Wall_RAD ; i-= Wall_RAD * 2)     // Layer 3
                {
                particle disk = particle(i, HEIGHT + Wall_RAD*(1+4*sin(M_PI/3.0)), Wall_RAD);
                disk.setWall(3,3);    // top wall tag = 3
                disks.push_back(disk);
                if (fabs(BELT_vel_x_top) == 0.0)
                ofp_wall_part_stat_3 << disk.position.x << "\t" << disk.position.y << "\t"<< disk.radius << "\t"<< disk.velocity.x << "\t" << disk.velocity.y << endl;
                }
            wall_limits.push_back(disks.size());
            unsigned short top_wall = wall_limits.back() - bot_wall - right_wall;
            ofp_wall << top_wall << endl;
            ofp_wall_part_stat_3.close();

            //            /***************** Left wall /
            for( float i = -Wall_RAD; i <= HEIGHT + Wall_RAD ; i+= Wall_RAD*2 )     // Layer 1
                {
                particle disk = particle(-Wall_RAD, i, Wall_RAD);
                disk.setWall(4,1);    // left wall tag = 4
                disks.push_back(disk);
                if (fabs(BELT_vel_y_left) == 0.0)
                ofp_wall_part_stat_4 << disk.position.x << "\t" << disk.position.y << "\t"<< disk.radius << "\t"<< disk.velocity.x << "\t" << disk.velocity.y << endl;
                }
            for( float i = 0.0; i <= HEIGHT ; i+= Wall_RAD*2 )     // Layer 2
                {
                particle disk = particle(-Wall_RAD*(1+2*sin(M_PI/3.0)), i, Wall_RAD);
                disk.setWall(4,2);    // left wall tag = 4
                disks.push_back(disk);
                if (fabs(BELT_vel_y_left) == 0.0)
                ofp_wall_part_stat_4 << disk.position.x << "\t" << disk.position.y << "\t"<< disk.radius << "\t"<< disk.velocity.x << "\t" << disk.velocity.y << endl;
                }
            for( float i = -Wall_RAD; i <= HEIGHT + Wall_RAD ; i+= Wall_RAD*2 )     // Layer 3
                {
                particle disk = particle(-Wall_RAD*(1+4*sin(M_PI/3.0)), i, Wall_RAD);
                disk.setWall(4,3);    // left wall tag = 4
                disks.push_back(disk);
                if (fabs(BELT_vel_y_left) == 0.0)
                ofp_wall_part_stat_4 << disk.position.x << "\t" << disk.position.y << "\t"<< disk.radius << "\t"<< disk.velocity.x << "\t" << disk.velocity.y << endl;
                }
            wall_limits.push_back(disks.size());
            unsigned short left_wall = wall_limits.back() - bot_wall - right_wall - top_wall;
            ofp_wall << left_wall << endl;
            ofp_wall.close();
            ofp_wall_part_stat_4.close();
        }
 */
        void boundaryParticlesRead()
        {
            int i = 0, wall_tot;
            float x,y,r,velx,vely;
//	    	int type;
            wall_limits.push_back(disks.size());
            ifstream ifp_wall("wall_ini.txt", ios::in);
            while(!ifp_wall.eof())  //until end-of-file, transfer the position and velocity data stored in the file to particle function
            {
                ifp_wall>>x>>y>>r>>velx>>vely; // >>type;
                particle disk = particle(x,y,r,velx,vely); //,type);
                disks.push_back(disk);
                i++;
            }
            disks.pop_back();
            wall_tot = disks.size();
            // Assigning wall tags in a single statement
            //for (int j=0; j <disks.size(); j++)
              //  disks[j].setWall((j+1)%(disks.size()/4),1);    // The read particles are assigned tags 1, 2, 3, 4

            // The read particles are assigned tags 1, 2, 3, 4
            // Assigning wall tags in 4 statements
            for (int j=0; j <disks.size()/4; j++)
                disks[j].setWall(1); //,1);
            for (int j=disks.size()/4; j <disks.size()/2; j++)
                disks[j].setWall(2); //,1);
            for (int j=disks.size()/2; j <disks.size()*3/4; j++)
                disks[j].setWall(3); //,1);
            for (int j=disks.size()*3/4; j <disks.size(); j++)
                disks[j].setWall(4); //,1);

            wall_limits.push_back(disks.size()/4);
            wall_limits.push_back(disks.size()/2);
            wall_limits.push_back(disks.size()*3/4);
            wall_limits.push_back(disks.size());
	    	cout << "\n ***** Wall particles' initial data imported.\n" << endl;
        }
/*
        void interiorParticlesGenerate(int i)
        {
            // N_s = N* (x / (1-x*(r/R)^2))
            int small_part = floor(CORE_PARTICLES*CONC/(1-CONC*(1-pow(Secondary_RAD/RAD,2))));

            // N_b = N* ( (1-x) / ( 1 - x * (1 - (r/R)^2) ) )
            int big_part = floor(CORE_PARTICLES*(1 - CONC)/(1-CONC*(1-pow(Secondary_RAD/RAD,2))));
            ofstream ofp_int("bot_top_small_big_number.txt", ios::out);

            bot_Core_small = floor(small_part/2);
            bot_Core_big = floor(big_part/2);
            top_Core_small = small_part - bot_Core_small;
            top_Core_big = big_part - bot_Core_big;

            for (int i = 0; i < bot_Core_small; i++)
                {
                particle disk = particle(HEIGHT/2.0);
                disk.setSpecies(1);
                disk.setRadius(Secondary_RAD);
                disk.setMass(SIGMA*M_PI*pow(Secondary_RAD,2));
                disk.setVsurr(M_PI*pow(5*2*Secondary_RAD,2));
                disks.push_back(disk);
                }
            ofp_int << bot_Core_small << endl;
            for (int i = bot_Core_small; i < bot_Core_small+bot_Core_big; i++)
                {
                particle disk = particle(HEIGHT/2.0);
                disk.setSpecies(1);
                disks.push_back(disk);
                }
            ofp_int << bot_Core_big << endl;
            for (int i = bot_Core_small+bot_Core_big; i < bot_Core_small+bot_Core_big + top_Core_small; i++)
                {
                particle disk = particle(HEIGHT/2.0);
                disk.addCy(HEIGHT/2.0);
                disk.setSpecies(2);
                disk.setRadius(Secondary_RAD);
                disk.setMass(SIGMA*M_PI*pow(Secondary_RAD,2));
                disk.setVsurr(M_PI*pow(5*2*Secondary_RAD,2));
                disks.push_back(disk);
                }
            ofp_int << top_Core_small << endl;
            for (int i = bot_Core_small+bot_Core_big + top_Core_small; i < bot_Core_small+bot_Core_big + top_Core_small + top_Core_big; i++)
                {
                particle disk = particle(HEIGHT/2.0);
                disk.addCy(HEIGHT/2.0);
                disk.setSpecies(2);
                disks.push_back(disk);
                }
            ofp_int << top_Core_big << endl;
            ofp_int.close();
            if (i>1)
            {
                float norm_area = CORE_PARTICLES * M_PI*RAD*RAD, tot_area;
                vector<float> r;
                if (i == 2)           // Uniform random
                {
                for (int i = 0; i< CORE_PARTICLES; i++)
                    r.push_back(rand01()*(BIG_RAD - SMALL_RAD) + SMALL_RAD);
                }

                if (i == 3)// Gaussian random
                {

                    float norm_area = CORE_PARTICLES * M_PI*RAD*RAD, tot_area;
                            // Gaussian random
                    float min = 0.0, max = 0.0;
                    for (int i = 0; i< CORE_PARTICLES; i++)
                    {
                        r.push_back(randGauss(1.0));
                        if (r.back() < min) min = r.back();
                        if (r.back() > max) max = r.back();
                    }
                    float mean = accumulate(r.begin(), r.end(), 0.0)/r.size();
                    float std_dev = 0.0;
                    for_each (r.begin(), r.end(), [&](const float rad) { std_dev += (rad - mean) * (rad - mean); });
                    std_dev = sqrt(std_dev/(r.size() - 1));

                    for (int i = 0; i< CORE_PARTICLES; i++)
                        r[i] = (r[i] - mean) / (max - min); // Mean normalization and scaling done

                    mean = accumulate(r.begin(), r.end(), 0.0)/r.size();
                    std_dev = 0.0;
                    for_each (r.begin(), r.end(), [&](const float rad) { std_dev += (rad - mean) * (rad - mean); });
                    std_dev = sqrt(std_dev/(r.size() - 1));

                    for (int i = 0; i< CORE_PARTICLES; i++)
                        r[i] = r[i]*gaussStdDevPercRAD * RAD/std_dev + RAD;
                }
                tot_area = accumulate(r.begin(), r.end(), 0.0, [](float total, float rad){return total + M_PI*rad*rad;});

                float max_RAD = RAD;
                for (int i = 0; i< CORE_PARTICLES; i++)
                {
                    r[i] *= sqrt(norm_area/tot_area);
                    disks[i+wall_limits.back()].setRadius(r[i]);
                    disks[i+wall_limits.back()].setMass(SIGMA*M_PI*r[i]*r[i]);
                    disks[i+wall_limits.back()].setVsurr(M_PI*pow(5*2*r[i],2));
                    if (r[i] > max_RAD)
                        max_RAD = r[i];
                }
                r_cutoff = 5.0*(2.0*max_RAD);
                r_skin = 0.5*RAD;
            }
        }
*/
        void interiorParticlesRead()
        {
            int i = 0;
            float x,y,r,velx,vely;
            int type;
            ifstream ifp_int("interior_ini.txt", ios::in);
            while(!ifp_int.eof())
            {
                ifp_int>>x>>y>>r>>velx>>vely>>type;
                particle disk = particle(x,y,r,velx,vely,type);
/*                if (i < bot_Core_small + bot_Core_big)
                    disk.setSpecies(1);
                else
                    disk.setSpecies(2);
*/                
				disks.push_back(disk);
                i++;
            }
            disks.pop_back();
	    	cout << "\n ***** Interior particles' initial data imported.\n" << endl;
        }
		
		void neighbour_list_simple()
		{
			nb_part.resize(disks.size());
			float dist, r_cutoff_loc;
			for (unsigned short i = 0; i < disks.size(); i++)
			{
				nb_part[i].clear();
			}
			for (unsigned short i = 0; i < disks.size()-1; i++)
			{
				if (disks[i].wall_tag == 0)	r_cutoff_loc = h_attract*2.0*disks[i].radius;
				else if (disks[i].wall_tag == -1)	r_cutoff_loc = h_detect*2.0*disks[i].radius;
				
				for (unsigned short j = i+1; j < disks.size(); j++)
				{
					if (disks[i].wall_tag <= 0 && disks[j].wall_tag <= 0)
					{
						dist = disks[i].position.distance_calc(disks[j].position).magnitude();
						if (dist <= r_cutoff_loc) nb_part[i].push_back(j);
					}
				}
			}
		}		

        void neighbour_list()
        {	
            float r_cutoff_loc;
          
	 		// Un-comment the below 4 lines in case the domain is expanding
            n = ceil((Xmax - Xmin)/(r_cutoff+r_skin)); // n has been updated
            m = ceil((Ymax - Ymin)/(r_cutoff+r_skin)); // m has been updated
            ncells = n*m; // ncells has been updated
            vector<int> head(ncells);
            vector<int> tail(disks.size());
                        
            nb_part.resize(disks.size());
            vector<vector<int>> cell_part(ncells);

            int cell_ix, cell_iy, cell_n;
            int l;

            for (unsigned short j = 0; j < ncells; j++)  head[j] = -1;

            // Head-tail construction
            for (int i = 0 ; i < disks.size(); i++)
                {
                nb_part[i].clear();
                cell_ix = floor((disks[i].position.x - Xmin)/(r_cutoff+r_skin));
                cell_iy = floor((disks[i].position.y - Ymin)/(r_cutoff+r_skin));
                cell_n = cell_ix + cell_iy * n;
                tail[i] = head[cell_n];
                head[cell_n] = i;
                }

            // particles in cell construction
         for (int i=0; i < ncells; i++)
         {
            if (head[i] > -1)
            {
                cell_part[i].clear(); //cell_part.shrink_to_fit();
                
                l=0;
                cell_part[i].push_back(head[i]);
                
                while (tail[cell_part[i][l]] > -1 && tail[cell_part[i][l]] < disks.size())
                {
                    cell_part[i].push_back(tail[cell_part[i][l]]);
                    l++;
                    
                }
                /*
                 * cell_part[i] contains the ids of particles in the i-th cell.
                 * The interior particles have higher ids and are more probable to be heads in a given cell.
                 */
            }
        }


            // Neigbourhood list construction
            int potential_nb_cell;
            vector <int> nb_cells;
            vector<int> nb_add(4);
            // Comment the following line for periodic domain (Updated on 08/08/20)
//            nb_add[0] = n; nb_add[1] = n+1; nb_add[2] = 1; nb_add[3] = -n+1;

            for (int i = 0; i < ncells; i++)
            {
                if (head[i] > -1)
                {

                    for (int j = 0; j< cell_part[i].size(); j++)
                    {
                        for (int k = j+1; k < cell_part[i].size(); k++)
                        {
                            if (disks[cell_part[i][j]].wall_tag <= 0)
							//if (disks[cell_part[i][j]].kind != -1) // In order to avoid checking for interaction between wall particles.
                            {
                            	if (disks[cell_part[i][j]].wall_tag == -1)	r_cutoff_loc = h_detect*2.0*disks[cell_part[i][j]].radius;
                            	else	r_cutoff_loc = h_attract*2.0*disks[cell_part[i][j]].radius;
                                if (disks[cell_part[i][j]].position.distance_calc(disks[cell_part[i][k]].position).magnitude() < r_cutoff_loc + r_skin)
                                    nb_part[cell_part[i][j]].push_back(cell_part[i][k]);
                            }
                        }
                    }
                    /************************************************************************

                        The id of wall_particles in the domain is lower than that of interior ones.
                        Therefore, if cell_part[i][j] is the id of a wall particle,
                        it should not check for interaction with the subsequent wall particles
                    ************************************************************************/

                    nb_cells.clear();
                    
                    // Comment the code below in case of periodic domain
/*                    
                    for (auto potential_nb_cell_add: nb_add)
                    {
                        potential_nb_cell = potential_nb_cell_add + i;
                        if ((potential_nb_cell < ncells) && (potential_nb_cell > -1))
                        {
                            if (i%n == 0)
                            {
                                nb_cells.push_back(potential_nb_cell);
                            }
                            else
                            {
                                if (potential_nb_cell % n != 0)
                                    nb_cells.push_back(potential_nb_cell);
                            }
                        }
                    }
*/
					/****Below code is for periodic domain. Added on 08/08/20. ****/
                
                
                	if (i == ncells-1)	// For the last cell
                	{
                    	nb_add[0] = -n+1; nb_add[1] = -2*n+1; nb_add[2] = -ncells+n; nb_add[3] = -ncells+1;
                	}
	                else if (i == n-1)	// For the last cell of the first row
	                {
	                    nb_add[0] = -n+1; nb_add[1] = n; nb_add[2] = 1; nb_add[3] = ncells-2*n+1;
	                }
	                else if (i != ncells-1 && i != n-1 && (i+1)%n == 0)	// For all the cells of the last column except the cells from first and last rows
	                {
	                    nb_add[0] = -n+1; nb_add[1] = n; nb_add[2] = 1; nb_add[3] = -2*n+1;
	                }
	                else if (i > ncells-n-1 && i != ncells-1)	// For the cells of the last row except the cell of the last column
	                {
	                    nb_add[0] = n+1; nb_add[1] = 1; nb_add[2] = -ncells+n; nb_add[3] = -ncells+n+1;
	                }
	                else if (i < n-1)	// For the cells of the first row except the cell of the last column
	                {
	                    nb_add[0] = n+1; nb_add[1] = n; nb_add[2] = 1; nb_add[3] = ncells-n+1;
	                }
	                else	// For all other cells
	                {
	                    nb_add[0] = n; nb_add[1] = n+1; nb_add[2] = 1; nb_add[3] = -n+1;
	                }
	                for (auto potential_nb_cell_add:nb_add)
	                {
	                    potential_nb_cell = potential_nb_cell_add + i;
	                    if (potential_nb_cell > -1 && potential_nb_cell < ncells)	nb_cells.push_back(potential_nb_cell);
	                }
	                /****Update for periodic domain ends here. ****/
	                
                    for (auto part: cell_part[i])
                    {
                        for (auto nb_cell: nb_cells)
                        {
                            for (auto neighbour_particle: cell_part[nb_cell])
                            {
                                if (!(disks[part].wall_tag > 0 && disks[neighbour_particle].wall_tag > 0))
								//if (!(disks[part].kind == -1 && disks[neighbour_particle].kind == -1))
                                    // In order to avoid checking for interaction between wall particles.
                                    // We require nb_part[part] to be populated basically; but don't this (populate nb_part[part]) only if both
                                    // the particle in the cell (part) and the neighbouring particle belong to the wall.
                                {
                                    if (disks[part].wall_tag == -1)	r_cutoff_loc = h_detect*2.0*disks[part].radius;
                            		else	r_cutoff_loc = h_attract*2.0*disks[part].radius;
//                                    if (r_cutoff_loc < SMALL_RAD + BIG_RAD) r_cut_loc = SMALL_RAD + BIG_RAD; // updating local r_cutoff in order to at least calculate the spring force
                                    if (disks[part].position.distance_calc(disks[neighbour_particle].position).magnitude() < r_cutoff_loc + r_skin)
                                        // The above conditional is for the current implementation of fluid drag force (5 times particle diameter)
                                        // In case, a common r_cutoff is to be used for drag calculation for all particles,
                                        // comment the above conditional and uncomment the below one.
//                                    if (disks[part].position.distance_calc(disks[neighbour_particle].position).magnitude() < r_cutoff + r_skin)
                                        nb_part[part].push_back(neighbour_particle);
                                }
                            }
                        }

                    }
                }
            }
        }
/*        void order_calc_init()
        {
            for(unsigned short i = wall_limits.back() ; i < disks.size(); i++)
            {
                disks[i].num_top = 0;
                disks[i].num_bot = 0;
                disks[i].num_surr = 0;
                disks[i].mass_top = 0.0;
                disks[i].mass_bot = 0.0;
                disks[i].mass_surr = 0.0;
                disks[i].local_order_param.x = 0.0;
                disks[i].local_order_param.y = 0.0;
            }
        }
        void order_calc()
        {
            for(int i = 0; i<disks.size(); i++)
                for (auto nb: nb_part[i])
                    disks[i].orderCalcOther(&disks[nb]);
        }
*/
        void force_initialize()
        {
//            if (VerletUpdate) sp_pot_energy = 0.0;

            for(unsigned short i = 0; i < disks.size(); i++)
            {
                    disks[i].force.x = disks[i].mass*gravity.x;
                    disks[i].force.y = disks[i].mass*gravity.y;
                    disks[i].vf.x = 0.0;
                    disks[i].vf.y = 0.0;
                    disks[i].vf_0 = 0.0;
                    disks[i].v_surr_part = 0.0;
            }
        }

        void force_pp()
        {
            for(int i = 0; i < disks.size(); i++)
            {
                for (auto nb: nb_part[i])
                {
                    disks[i].separate_other(&disks[nb]); //, true);	// For all particles
                    if (disks[i].wall_tag == 0 && disks[nb].wall_tag == 0 && disks[i].id > 0 && disks[nb].id > 0)	disks[i].weight_other(&disks[nb]);	// Only for pair of alive prey
            	}
            }
        }

        void force_align()
        {
            float mfrac = 1.0, m_fluid, m_surr;
            vect fa; double theta_align;
            noise_val.resize(disks.size());
            gaussian_white(noise);
            for(unsigned short i = 0; i < disks.size(); i++)
            {
				if(disks[i].id > 0 && disks[i].wall_tag == 0)
				{
		   			m_fluid = RHO * (disks[i].V_surr - disks[i].v_surr_part);
                    m_surr = SIGMA * disks[i].v_surr_part;
                    if (m_fluid<0)
                    {
                        m_fluid=0.0;
                    }
                    mfrac = m_surr/(m_surr + m_fluid);

                    if (disks[i].vf_0 > 0)
                    {
                        disks[i].vf.x /= disks[i].vf_0;
                        disks[i].vf.y /= disks[i].vf_0;
                    }
                    else
                    {
                        disks[i].vf.x = 0;
                        disks[i].vf.y = 0;
                    }

		    		fa.x = Cv*(2*disks[i].radius)*(mfrac * disks[i].vf.x - disks[i].velocity.x);
            	    fa.y = Cv*(2*disks[i].radius)*(mfrac * disks[i].vf.y - disks[i].velocity.y);
            
            	    theta_align = atan2(fa.y,fa.x) + noise_val[i];
            
            	    disks[i].force.x += fa.magnitude() * cos(theta_align);				
            	    disks[i].force.y += fa.magnitude() * sin(theta_align);
				}
				else if (disks[i].id > 0 && disks[i].wall_tag == -1)
				{
					// Add coordination for predators here or friction in case of non-coordination
				}
            }
        }

	    void force_sp()
        {
            float vmag;
            for(unsigned short i = 0; i < disks.size(); i++)
            {
                if (escape_list.find(i) == escape_list.end() && disks[i].id > 0)
                {
		            vmag = disks[i].velocity.magnitude();
                    if (vmag == 0.0) vmag = 1.0;
                    disks[i].force.x += disks[i].mass*disks[i].beta*disks[i].velocity.x/vmag;
                    disks[i].force.y += disks[i].mass*disks[i].beta*disks[i].velocity.y/vmag;
                }
                force_fric(i);
            }
        }
        void force_sp(auto pred_id)
        {
            float vmag;
		    vmag = disks[pred_id].velocity.magnitude();
            if (vmag == 0.0) vmag = 1.0;
            disks[pred_id].force.x += disks[pred_id].mass*disks[pred_id].beta*disks[pred_id].velocity.x/vmag;
            disks[pred_id].force.y += disks[pred_id].mass*disks[pred_id].beta*disks[pred_id].velocity.y/vmag;
        }
        
		void force_attract()
		{
			vect dist, f;
			float dist_mag, r_attract, r_align;
			for(int i = 0; i < disks.size(); i++)
            {
            	if (disks[i].wall_tag == 0 && disks[i].id > 0)	// If ith particle is a live prey
            	{
		        	r_attract = h_attract*2.0*disks[i].radius;
		        	r_align = h_align*2.0*disks[i].radius;
		            for (auto nb: nb_part[i])
		            {
		            	dist = disks[i].position.distance_calc(disks[nb].position); 
		            	dist_mag = dist.magnitude();
		            	if ((dist_mag > r_align && dist_mag <= r_attract) && (disks[nb].wall_tag == 0 && disks[nb].id > 0))
		            	{
				            f.x = F_atr_max * sqrt(1 - pow(1 - (dist_mag/r_attract),2)) * dist.x / dist_mag;
				            f.y = F_atr_max * sqrt(1 - pow(1 - (dist_mag/r_attract),2)) * dist.y / dist_mag;
				            disks[i].force.x += f.x;
				            disks[i].force.y += f.y;
				            disks[nb].force.x -= f.x;
				            disks[nb].force.y -= f.y;
				        }
		        	}
		        }
            }
		}         
            
        void force_pursuit(auto hunter, auto target)
        {
        	vect dist = disks[hunter].position.distance_calc(disks[target].position);
            disks[hunter].force.x += ( beta_hunt - fric_factor * pow(disks[hunter].velocity.magnitude(),2)) * disks[hunter].mass * dist.x / dist.magnitude();
            disks[hunter].force.y += ( beta_hunt - fric_factor * pow(disks[hunter].velocity.magnitude(),2)) * disks[hunter].mass * dist.y / dist.magnitude();
        }

        void force_escape(auto pursuer, auto own)
        {
        	vect dist = disks[pursuer].position.distance_calc(disks[own].position);
        	vect f;
        	int sign;
            f.x = (beta_escape - fric_factor * pow(disks[own].velocity.magnitude(),2)) * disks[own].mass * dist.x / dist.magnitude();
            f.y = (beta_escape - fric_factor * pow(disks[own].velocity.magnitude(),2)) * disks[own].mass * dist.y / dist.magnitude();
            
            // Algorithm for condition when prey is sandwiched between predator and wall (only for rectangular domain)
/*            sign = signum(dist.x * dist.y);
            if ((disks[own].position.x < (X0 + wall_avoid_dist) || disks[own].position.x > (X0 + WIDTH - wall_avoid_dist)) && fabs(dist.y) < pred_avoid_dist)
            {
            	if (sign >= 0)	f.rotate_vect(true);
            	else	f.rotate_vect(false);
            }
            
            else if ((disks[own].position.y < (Y0 + wall_avoid_dist) || disks[own].position.y > (Y0 + HEIGHT - wall_avoid_dist)) && fabs(dist.x) < pred_avoid_dist)
            {
            	if (sign < 0)	f.rotate_vect(true);
            	else	f.rotate_vect(false);
            }	
            
*/          disks[own].force.x += f.x;
            disks[own].force.y += f.y;
        }

/*
        void force_predator_repulse(int pred1, int pred2)
        {
            vect f, pdist;
            pdist = disks[pred1].position.distance_calc(disks[pred2].position);
            f.x = beta_repulse * (1 - (pdist.magnitude()/r_cutoff)) * pdist.x / pdist.magnitude();
            f.y = beta_repulse * (1 - (pdist.magnitude()/r_cutoff)) * pdist.y / pdist.magnitude();
	
		    disks[pred2].force.x += f.x;
		    disks[pred2].force.y += f.y;
				
		    disks[pred1].force.x -= f.x;
		    disks[pred1].force.y -= f.y;
		}
*/
	void force_fric(auto i)
	{
//		if (disks[i].wall_tag == 0)
//		{
			float circum = 2*M_PI*disks[i].radius;
			disks[i].force.x -= 0.02*RHO*circum*disks[i].velocity.magnitude()*disks[i].velocity.x;
//			disks[i].force.x -= 0.125*Cv*disks[i].radius*disks[i].velocity.x;
			disks[i].force.y -= 0.02*RHO*circum*disks[i].velocity.magnitude()*disks[i].velocity.y;
//			disks[i].force.y -= 0.125*Cv*disks[i].radius*disks[i].velocity.y;
//		}
//		else if (disks[i].wall_tag == -1)
//		{
//			disks[i].force.x -= 0.25*Cv*disks[i].radius*disks[i].velocity.x;
//			disks[i].force.y -= 0.25*Cv*disks[i].radius*disks[i].velocity.y;
//		}
	}
/*
	void calcDisplace()
	{
		displace_max.x = 0.0;		
		for(unsigned short i = wall_limits.back(); i < disks.size(); i++)
		{
			disks[i].setDisplace(disks[i].position.distance_calc(disks[i].position_prev).magnitude());
                	/-------Additional step to calculate the maximum displacement among all the particles------/
                   	if (disks[i].displace > displace_max.x)
                    	{
                        	displace_max.x = disks[i].displace;
                        	displace_max.y = disks[i].id;
                    	}
                	
		}
		//cout<<"Maximum Displacement = "<<displace_max.x<<endl;
	}	

        void force_update(bool VerletUpdate)
        {
                if(VerletUpdate)
                {
                    for(short unsigned i = wall_limits.back(); i < disks.size(); i++)
                    {
                        disks[i].force_prev.x = disks[i].force.x;
                        disks[i].force_prev.y = disks[i].force.y;
                    }
                }
        }
*/
	
	void make_predating_list()
	{
		float dist;
		predating_list.resize(hunt_list.size());
		escape_list.clear();
		for (short unsigned i = 0; i < predating_list.size(); i++)
		{
			predating_list[i].clear();
		}
		for (short unsigned i = 0; i < nb_part.size(); i++)
		{
			for (auto nb:nb_part[i])
			{
				if (disks[nb].wall_tag == -1 && disks[i].wall_tag == 0 && disks[i].id > 0) // i == alive prey and nb == predator
				{
					dist = disks[i].position.distance_calc(disks[nb].position).magnitude();
					if (dist < (h_escape*2.0*disks[i].radius) && disks[i].vision_check(&disks[nb]))
					{
						force_escape(nb, i);	// Only if the predator is inside escape zone of prey
						escape_list.insert(make_pair(i, 'p'));
					}
					if (disks[nb].vision_check(&disks[i]) && dist < (h_detect*2.0*disks[nb].radius))	// Checking if within vision zone of predator
					{
						for (short unsigned k = 0; k < hunt_list.size(); k++)
						{
							if (hunt_list[k][0] == disks[nb].id && hunt_list[k][1] == -1 && hunt_list[k][3] == 0 && hunt_list[k][4] == 0) // matching predator ID, checking for no active target prey and not in satisfaction state or refocus state
							{
								predating_list[k].push_back(i);
							}
						}
					}
				}
				else if (disks[nb].wall_tag == 0 && disks[i].wall_tag == -1 && disks[nb].id > 0) // i == predator and nb == alive prey
				{
					dist = disks[nb].position.distance_calc(disks[i].position).magnitude();
					if (dist < (h_escape*2.0*disks[nb].radius) && disks[nb].vision_check(&disks[i]))
					{
						force_escape(i, nb);	// Only if the predator is inside escape zone of prey
						escape_list.insert(make_pair(nb, 'p'));
					}
					if (disks[i].vision_check(&disks[nb]) && dist < (h_detect*2.0*disks[i].radius))	// Checking if within vision zone of predator
					{
						for (short unsigned k = 0; k < hunt_list.size(); k++)
						{
							if (hunt_list[k][0] == disks[i].id && hunt_list[k][1] == -1 && hunt_list[k][3] == 0 && hunt_list[k][4] == 0) // matching predator ID and checking if any active target prey and not in satisfaction state or refocus state
							{
								predating_list[k].push_back(nb);
							}
						}
					}
				}
				else if (disks[nb].wall_tag == -1 && disks[i].wall_tag == -1) // both are predators
				{
//					force_predator_repulse(i, nb);
//					force_fric(i);
				}
			}
		}
	}
	void target_select()
	{
		int tr_prey = -1;
		for (short unsigned pred = 0; pred < hunt_list.size(); pred++)
		{
			if (hunt_list[pred][1] == -1 && hunt_list[pred][3] == 0 && hunt_list[pred][4] == 0) // No active target prey, not in satisfaction state, not in refocus state
			{
				short unsigned key = strategy_select(hunt_list[pred][0]);
				switch(key)
				{
					case 1:
						tr_prey = central_select(pred);
						if (tr_prey != -1)	
						{
							cout << "Central prey selected.\n" << endl;
							hunt_list[pred][5] = 1;
						}
						break;
					
					case 2:
						tr_prey = isolated_select(pred);
						if (tr_prey != -1)
						{
							cout << "Isolated prey selected.\n" << endl;
							hunt_list[pred][5] = 2;
						}
						else
						{
							tr_prey = nearest_select(pred);
							cout << "Isolated prey not available. Defaulting to nearest prey.\n" << endl;
							hunt_list[pred][5] = 3;
						}
						break;
						
					case 3:
						tr_prey = peripheral_select(pred);
						if (tr_prey != -1)
						{
							cout << "Peripheral prey selected.\n" << endl;
							hunt_list[pred][5] = 3;
						}
						break;
					case 4:
					default:
						tr_prey = nearest_select(pred);
						if (tr_prey != -1)	
						{
							cout << "Nearest prey selected.\n" << endl;
							hunt_list[pred][5] = 4;
						}
				}
				
				if (tr_prey != -1)
				{
					force_pursuit(0 + hunt_list[pred][0] - 1, tr_prey);
					hunt_list[pred][1] = tr_prey;	// Saving active prey info
					hunt_list[pred][2] += 1;	// Setting hunt_timer to zero
				}
			}
			else if (hunt_list[pred][1] != -1 && hunt_list[pred][3] == 0 && hunt_list[pred][4] == 0)  // Actively pursuing a prey
			{
				force_pursuit(0 + hunt_list[pred][0] - 1, hunt_list[pred][1]);
			}
		}
	}
	short unsigned strategy_select(int pred_id)
	{
		if (pred_id == 1 || pred_id == 2)	return 4;	// Both first and second predators choose nearest strategy
		else return 4;
	}
	
	int nearest_select(short unsigned pred)
	{
		// Nearest prey selection routine
		float closest = WIDTH; int tr_prey = -1;
		for (auto prey:predating_list[pred])
		{
			float dist = disks[0 + hunt_list[pred][0] - 1].position.distance_calc(disks[prey].position).magnitude();
			if (dist < closest)
			{
				closest = dist;
				tr_prey = prey;
			}
		}
		return tr_prey;
	}
	int central_select(short unsigned pred)
	{
		// Central prey selection routine
		vect avg_pos; int tr_prey = -1;
		float dist;
		// Finding the average position of the neighbours
		for (auto prey:predating_list[pred])
		{
			avg_pos.x += disks[prey].position.x;
			avg_pos.y += disks[prey].position.y;
		}
	
		avg_pos.x /= predating_list[pred].size();
		avg_pos.y /= predating_list[pred].size();
		
		float closest = WIDTH;
		for (auto prey:predating_list[pred])
		{
			dist = disks[prey].position.distance_calc(avg_pos).magnitude();
			if (dist < closest)
			{
				closest = dist;
				tr_prey = prey;
			}
		}
		return tr_prey;
	}
	int isolated_select(short unsigned pred)
	{
		// Isolated prey selection routine
		int nb_prey, tr_prey = -1;
		for (auto prey:predating_list[pred])
		{
			nb_prey = 0;
			for (short unsigned i = 0; i < disks.size(); i++)
			{
				if (i != prey)
				{
					for (auto nb:nb_part[i])
					{
						if (disks[i].wall_tag == 0 && disks[i].id > 0 && nb == prey)	
						{
							nb_prey++;
							break;
						}
					}
				}
				else
				{
					for (auto nb:nb_part[prey])
					{
						if (disks[nb].wall_tag == 0 && disks[nb].id > 0)
						{
							nb_prey++;
							break;
						}
					}
					break;
				}
			}
			if (nb_prey == 0)
			{
				tr_prey = prey;
				break;
			}
		}
		return tr_prey;
	}
	int peripheral_select(short unsigned pred)
	{
		// Peripheral prey selection routine
		int tr_prey = -1;
		float max_angle = 0.0, angle_diff;
		
		// Finding the max abs angle between the positional difference vector and the predator velocity 
		for (auto prey:predating_list[pred])
		{	
			angle_diff = disks[hunt_list[pred][0] - 1].angle_diff_calc(&disks[prey]);
			
			if (angle_diff > max_angle)
			{
				max_angle = angle_diff;
				tr_prey = prey;
			}
		}
		return tr_prey;
	}
	void hunt_result(int iter)
	{
		for (short unsigned pred = 0; pred < hunt_list.size(); pred++)
		{
			// Actively hunting with non-expired hunting time and not in satisfaction state, not in refocus stat0e
			if (hunt_list[pred][1] != -1 && hunt_list[pred][2] > 0 && hunt_list[pred][2] < hunt_timelimit && hunt_list[pred][3] == 0 && hunt_list[pred][4] == 0) 
			{
				// Target prey within kill radius
				if (disks[0 + hunt_list[pred][0] - 1].position.distance_calc(disks[hunt_list[pred][1]].position).magnitude() <= r_kill)
				{
					
					// Printing to file. Added on 21-01-21 by SM.
								
					ofstream ofp2("predator_state.txt", std::ios_base::app);
					ofp2 << hunt_list[pred][0] << "\t" << (iter*dt) << "\t" << 0 << "\t" << 1 << "\t" << hunt_list[pred][1]+1 << "\t" << hunt_list[pred][5] << endl;  // Pred_ID, time elapsed, prev. state, new state, target ID (if any), strategy select.
					// State value change from 0 to 1 represents entering satisfaction state from hunting state
					ofp2.close();

					disks[hunt_list[pred][1]].id = -iter;	 // Setting the global ID of killed prey to -iter 
					hunt_list[pred][1] = -1;
					hunt_list[pred][2] = 0;		// Hunt timer reset
					hunt_list[pred][3] += 1;	// Satisfaction state started
					hunt_list[pred][5] = 0;		// Strategy select reset to NA
//					hunt_list[pred][6] += 1;	// No. of successful hunts +1
				    	
				}
				// Target prey outside kill radius
				else
				{
					hunt_list[pred][2] += 1;	// Hunt time elapsed += dt
				}
			}
			// Actively hunting with expired hunting time and not in satisfaction state
			else if (hunt_list[pred][1] != -1 && hunt_list[pred][2] >= hunt_timelimit && hunt_list[pred][3] == 0 && hunt_list[pred][4] == 0) 
			{	
				// Printing to file. Added on 21-01-21 by SM.

				ofstream ofp2("predator_state.txt", std::ios_base::app);
				ofp2 << hunt_list[pred][0] << "\t" << (iter*dt) << "\t" << 0 << "\t" << -1 << "\t" << hunt_list[pred][1]+1 << "\t" << hunt_list[pred][5] << endl;  // Pred_ID, time elapsed, prev. state, new state, target ID (if any). Target ID of -1 shows no present target.
				// State value change from 0 to -1 represents entering refocus state from hunting state
				ofp2.close();
				
				hunt_list[pred][1] = -1;
				hunt_list[pred][2] = 0;		// Hunt timer reset
				hunt_list[pred][4] += 1;	// Refocus timer started
				hunt_list[pred][5] = 0;		// Strategy reset to NA
//				hunt_list[pred][7] += 1;	// No. of unsuccessful hunts +1

			}
			// Hunting inactive, and satisfaction state time limit reached, not in refocus state
			else if (hunt_list[pred][1] == -1 && hunt_list[pred][2] == 0 && hunt_list[pred][3] >= satisfy_timelimit && hunt_list[pred][4] == 0) 
			{
				hunt_list[pred][3] = 0;		// Satisfaction timer reset
				hunt_list[pred][2] += 1;	// Hunt phase started

				// Printing to file. Added on 21-01-21 by SM.

				ofstream ofp2("predator_state.txt", std::ios_base::app);
				ofp2 << hunt_list[pred][0] << "\t" << (iter*dt) << "\t" << 1 << "\t" << 0 << "\t" << -1 << "\t" << hunt_list[pred][5] << endl;  // Pred_ID, time elapsed, prev. state, new state, target ID (if any), strategy select (if applicable). Target ID of -1 shows target not yet selected. Strategy select of 0 here shows no strategy selected yet.
				// State value change from 1 to 0 represents entering hunting state from satisfaction state
				ofp2.close();
			}
			// Hunting inactive, and non-expired satisfaction state time limit, not in refocus state
			else if (hunt_list[pred][1] == -1 && hunt_list[pred][2] == 0 && hunt_list[pred][3] > 0 && hunt_list[pred][3] < satisfy_timelimit && hunt_list[pred][4] == 0) 
			{
				hunt_list[pred][3] += 1; 	// Satisfaction time elapsed += dt
				force_sp(0 + hunt_list[pred][0] - 1);
			}
			// Hunting inactive with non-expired refocus time, not in satisfaction state
			else if (hunt_list[pred][1] == -1 && hunt_list[pred][2] == 0 && hunt_list[pred][3] == 0 && hunt_list[pred][4] > 0 && hunt_list[pred][4] < refocus_timelimit) 
			{
				hunt_list[pred][4] += 1;	// Refocus time elapsed += dt
				force_sp(0 + hunt_list[pred][0] - 1);
			}
			// Hunting inactive with expired refocus time, not in satisfaction state
			else if (hunt_list[pred][1] == -1 && hunt_list[pred][2] == 0 && hunt_list[pred][3] == 0 && hunt_list[pred][4] >= refocus_timelimit) 
			{
				hunt_list[pred][4] = 0;		// Refocus timer reset
				hunt_list[pred][2] += 1;	// Hunt phase started

				// Printing to file. Added on 21-01-21 by SM.

				ofstream ofp2("predator_state.txt", std::ios_base::app);
				ofp2 << hunt_list[pred][0] << "\t" << (iter*dt) << "\t" << -1 << "\t" << 0 << "\t" << -1 << "\t" << hunt_list[pred][5] << endl;  // Pred_ID, time elapsed, prev. state, new state, target ID (if any), strategy select (if applicable). Target ID of -1 shows target not yet selected. Strategy select of 0 here shows no strategy selected yet.
				// State value change from -1 to 0 represents entering hunting state from refocus state
				ofp2.close();
			}
		}
	}	
	
    void render(int iter) 
    {
//        int n_active;
	    for(unsigned short i = 0; i < disks.size(); i++)
	    {	    
			disks[i].setPrevPosx(disks[i].position.x);
            disks[i].setPrevPosy(disks[i].position.y);	
			disks[i].update_pos(dt);
	    }
//	    calcDisplace();
	    if (iter >= (ceil(time_HUNT_starts/dt)+1))
	    {
	    	hunt_result(iter);
			for(unsigned short i = 0; i < disks.size(); i++)
			{
				if (disks[i].id != (-1*iter)) 
				{
					disks[i].update_vel(dt);
				}
			}	
	    }
	    else
	    {
			for(unsigned short i = 0; i < disks.size(); i++)
			{
				disks[i].update_vel(dt);
			}
  	    }
	    force_initialize();
//          	if (displace_max.x > 0.99*r_skin) // Update the neighbour list only if the maximum displacement is close to the "skin" thickness for the VERLET list
//	    	if (iter % nb_list_monitor == 0)
//               	{
            	//cout << "Max diplacement = " << displace_max.x <<" at " << displace_max.y << " Wall tag = " << disks[displace_max.y].wall_tag<< endl;
//	    neighbour_list();            	
		neighbour_list_simple();
//		nb_list_calc++;
//          }   
	    force_pp();
	    if (iter >= ceil(time_SYS_starts/dt))
	    {	       
		    force_sp();
//     		force_fric();
        	force_align();
        	force_attract();
	    }
        if (iter >= ceil(time_HUNT_starts/dt))
        {
/*        	check_status();
        	n_active = 0;
        	for(short unsigned i = wall_limits.back(); i < disks.size(); i++)
       		{
       			if (disks[i].id > 0 && disks[i].state == 0) n_active++;
			}
			activate_cann(n_active);
			activate_escape();	
*/
			make_predating_list();
			target_select();
	    }
            
/*            power_sum = 0;
            for(unsigned short i = wall_limits.front(); i < wall_limits.back(); i++)
            {
                if (fabs(disks[i].velocity.x) > 0 || fabs(disks[i].velocity.y) > 0)
                    disks[i].update_wall_pos(dt);
                power_sum += (disks[i].velocity.x * disks[i].force.x + disks[i].velocity.y * disks[i].force.y);
                disks[i].setDisplace(disks[i].position.distance_calc(disks[i].position_prev).magnitude());

            }
*/    
	    for(unsigned short i = 0; i < disks.size(); i++)
        {
            disks[i].update_vel(dt);
        }

    }
    void mainLoop()
    {
//        int n_alive = 0;
    	char time_data[50], pred_data[50];
        int start = 0;
	float tot_simul_time = 0.0;

        double duration;// mix_num = 0.0, mix_mass = 0.0, global_order_param = 0.0, local_param_mag = 0.0;
        int num_iterations = ceil((t_final - t_initial)/dt);
        cout << "Total iterations = " << num_iterations << endl;
	
	// Below lines are for storing the initial velocity of the particles
		uniform_random(0.0, vmag_init);
		int count = 1, pred_count = 0;
		for(short unsigned i = 0; i < disks.size(); i++)
		{
			disks[i].id  = count;	// Particle IDs start from 1 and exist for all interior particles
	//		satisfy_list[k-1]  = rand02(satisfy_timelimit - 1); //All are satisfied initially, but in different stages of satisfaction, which are randomly alloted
			count++;
			disks[i].velocity.x = vel_init[i] * cos(rand02(2*M_PI)); //Giving a velocity of constant magnitude and random direction to the particles initially 
			disks[i].velocity.y = vel_init[i] * sin(rand02(2*M_PI));
			if (disks[i].wall_tag == -1)
			{
				pred_count++;
			}
		}
		vel_init.clear();
		
		hunt_list.resize(pred_count);
		pred_count = 0;
		for(short unsigned i = 0; i < disks.size(); i++)
		{
			if (disks[i].wall_tag == -1)
			{
				hunt_list[pred_count].push_back(disks[i].id);	// Predator ID				
				hunt_list[pred_count].push_back(-1);  // Active target prey global ID. -1 represents no active target prey
				hunt_list[pred_count].push_back(0);	// Hunt time elapsed. 0 when not hunting
				hunt_list[pred_count].push_back(rand02(satisfy_timelimit - 1));	  // Satisfaction time elapsed. 0 when not satisfied
				hunt_list[pred_count].push_back(0);	// Refocus time set to 0
				hunt_list[pred_count].push_back(0);	// Strategy select set to 0 (default meaning NA)
	//			hunt_list[pred_count].push_back(0);	// No. of successful hunts
	//			hunt_list[pred_count].push_back(0);	// No. of unsuccessful hunts
				pred_count++;
			}
		}
		if (pred_count > 0)
		{
			// Printing to file. Added on 21-01-21 by SM.
		
			ofstream ofp2("predator_state.txt", std::ios_base::out);
			ofp2 << "ID\tT\tP.St.\tC.St.\tT_id\tPr.Str.\n" << endl;  // Pred_ID, time elapsed, prev. state, new state, target ID (if any), strategy select (if applicable)
			ofp2.close();
		
			cout << "No. of predators = " << pred_count << endl;
		}
	
/*        neighbour_list();
        force_initialize();
        force_pp();
	if(time_SYS_starts == 0.0) force_sp();
        order_calc_init();
        order_calc();
        nb_list_calc++;
*/

        for (int j = 0; j <= num_iterations; j++)
        {
		
//            wall_energy_input += (-power_sum)*dt;

            if ((Cv>0) && ((time > time_SYS_starts && time < time_SYS_starts + 2*dt) || (t_initial > time_SYS_starts && time < t_initial + 2*dt)))
//	    if ((time > time_SYS_starts && time < time_SYS_starts + 2*dt) || (t_initial > time_SYS_starts && time < t_initial + 2*dt))
            {
                for(unsigned short i = 0; i < disks.size(); i++)
                {
                    if (disks[i].wall_tag == 0)	disks[i].setBeta(beta_escape/2);
                    else if (disks[i].wall_tag == -1)	disks[i].setBeta(beta_hunt/2);
                }
                cout << "Self_propulsion of prey particles starts with BETA = " << beta_escape/2 << endl;
                cout << "Self_propulsion of predator particles starts with BETA = " << beta_hunt/2 << endl;
                cout << "Coordination set to = " << Cv << endl << endl;
            }
	    if (j == ceil(time_HUNT_starts/dt))	cout << "Escape-and-pursuit starts"<<endl<<endl;
/*
            if ((time > time_BELT_starts && time < time_BELT_starts + 2*dt) || (t_initial > time_BELT_starts && time < t_initial + 2*dt))
            {

                for(unsigned short i = wall_limits.front(); i < wall_limits.back(); i++)
                {

                    if (disks[i].wall_tag == 1) // this is redundant if only one wall is moving
                    {

                        disks[i].setVelx(BELT_vel_x_bot);
                    }

                    if (disks[i].wall_tag == 2) // this is redundant if only one wall is moving
                        disks[i].setVely(BELT_vel_y_right);

                    if (disks[i].wall_tag == 3) // this is redundant if only one wall is moving
                        disks[i].setVelx(BELT_vel_x_top);

                    if (disks[i].wall_tag == 4) // this is redundant if only one wall is moving
                        disks[i].setVely(BELT_vel_y_left);

                }
                cout << "Belt starts " << endl;
            }
*/
	   if (j % monitor == 0 )
           {
//                float grav_potential_energy =0.0, tot_kinetic_energy = 0.0;


                snprintf(time_data,sizeof time_data, "./particle_data/interior_%.2f.txt",time);
                ofstream ofp1(time_data, ios::out);  // ofp is the identifier. ios::out means it is for writing
//                order_calc_init();
//                order_calc();
                for(unsigned short i = 0; i < disks.size(); i++)
                {
/*                    grav_potential_energy += disks[i].mass * GRAV * disks[i].position.y;
                    tot_kinetic_energy += 0.5 * disks[i].mass * pow(disks[i].velocity.magnitude(),2);

                    if (disks[i].num_surr > 0)
                    {
                        disks[i].local_order_param.x /= disks[i].num_surr;
                        disks[i].local_order_param.y /= disks[i].num_surr;
                        mix_num += pow( (disks[i].num_top - disks[i].num_bot)/disks[i].num_surr , 2);
                    }
                    if (disks[i].mass_surr > 0.0)
                    {
                        mix_mass += pow( (disks[i].mass_top - disks[i].mass_bot)/disks[i].mass_surr , 2);
                    }
                    local_param_mag = disks[i].local_order_param.magnitude();
*/                    
		    ofp1 << disks[i].id << "\t" << disks[i].position.x << "\t" << disks[i].position.y << "\t"<< disks[i].radius << "\t"<< disks[i].velocity.x << "\t" << disks[i].velocity.y << "\t" << disks[i].force.x << "\t" << disks[i].force.y << endl;
//                    global_order_param += local_param_mag;
                }
                ofp1.close();
/*                mix_num = sqrt(mix_num/(disks.size() - wall_limits.back()));
                mix_mass = sqrt(mix_mass/(disks.size() - wall_limits.back()));
                global_order_param /= (disks.size() - wall_limits.back());

                ofp_monitor<<time<<"\t"<<mix_num<<"\t"<<mix_mass<<"\t"<<global_order_param<<endl;

                snprintf(time_data_1,sizeof time_data_1, "./particle_data/wall_1_%.2f.txt",time);
                snprintf(time_data_2,sizeof time_data_2, "./particle_data/wall_2_%.2f.txt",time);
                snprintf(time_data_3,sizeof time_data_3, "./particle_data/wall_3_%.2f.txt",time);
                snprintf(time_data_4,sizeof time_data_4, "./particle_data/wall_4_%.2f.txt",time);

                ofstream ofp2_1(time_data_1, ios::out);
                ofstream ofp2_2(time_data_2, ios::out);
                ofstream ofp2_3(time_data_3, ios::out);
                ofstream ofp2_4(time_data_4, ios::out);


                for(unsigned short i = wall_limits.front(); i < wall_limits.back(); i++)
                {


                    if (BELT_vel_x_bot != 0)
                    {
                        if (disks[i].wall_tag == 1)
                            ofp2_1 << disks[i].position.x << "\t" << disks[i].position.y << "\t"<< disks[i].radius << "\t"<< disks[i].velocity.x << "\t" << disks[i].velocity.y << endl;
                    }
                    if (BELT_vel_y_right != 0)
                    {
                        if (disks[i].wall_tag == 2)
                            ofp2_2 << disks[i].position.x << "\t" << disks[i].position.y << "\t"<< disks[i].radius << "\t"<< disks[i].velocity.x << "\t" << disks[i].velocity.y << endl;
                    }
                    if (BELT_vel_x_top != 0)
                    {
                        if (disks[i].wall_tag == 3)
                            ofp2_3 << disks[i].position.x << "\t" << disks[i].position.y << "\t"<< disks[i].radius << "\t"<< disks[i].velocity.x << "\t" << disks[i].velocity.y << endl;
                    }
                    if (BELT_vel_y_left != 0)
                    {
                        if (disks[i].wall_tag == 4)
                            ofp2_4 << disks[i].position.x << "\t" << disks[i].position.y << "\t"<< disks[i].radius << "\t"<< disks[i].velocity.x << "\t" << disks[i].velocity.y << endl;
                    }
                }
                ofp2_1.close();
                ofp2_2.close();
                ofp2_3.close();
                ofp2_4.close();

                ofp_energy << time << "\t" << -power_sum << "\t" << wall_energy_input << "\t" << grav_potential_energy << "\t" << tot_kinetic_energy << "\t" << sp_pot_energy << "\n";
*/
		count = 0;
		for (short unsigned  i = 0; i < disks.size(); i++)
		{
			if (disks[i].id > 0 && disks[i].wall_tag == 0)	count++;
		}
                start=stop;
                stop=clock();
                duration = (double)(stop - start) / CLOCKS_PER_SEC;
				tot_simul_time += duration;
                cout << "Time = " << time << "\t" << "Simulation time = " << duration << endl;
                cout << "Total no. of particles = " << disks.size() << endl;
                
                if (pred_count > 0)	cout << "No. of prey agents alive = " << count << endl;
                cout << endl;
//                cout << "Total Verlet list updates = " << nb_list_calc << endl<< endl;
//                nb_list_calc = 0;

		// Printing to file. Added on 21-01-21 by SM.

		if (pred_count > 0 &&  j != 0)
	        {
			for (short unsigned i = 0; i < hunt_list.size(); i++)
			{
				snprintf(pred_data, sizeof pred_data, "predator_%d_stats.txt", (i+1));			
				ofstream ofp2(pred_data, std::ios_base::app);
		    	ofp2 << time << "\t" << hunt_list[i][0] << "\t" << hunt_list[i][1]+1 << endl;
				ofp2.close();
		    	}
		}
		else if (pred_count > 0 &&  j == 0)
		{
			for (short unsigned i = 0; i < hunt_list.size(); i++)
			{
				snprintf(pred_data, sizeof pred_data, "predator_%d_stats.txt", (i+1));			
				ofstream ofp2(pred_data, std::ios_base::out);
		    		ofp2 << "Time\tPr_ID\tTr_ID" << endl << endl;
				ofp2.close();
				cout << "Predator " << (i+1) << " ID is:" << hunt_list[i][0] << endl << endl;
		    	}
		}

            }

            time = (j+1)*dt + t_initial;
            // modify for oscillation

/*		 if ((time > time_BELT_starts )&& (j % monitor == 0 ))
            	{
            		//BELT_vel_x_bot= 0.5 - 1*sin(time/2); // modify the belt velocity
		        for(unsigned short i = wall_limits.front(); i < wall_limits.back(); i++)
		        {

		            if (disks[i].wall_tag == 1) // this is redundant if only one wall is moving
		            {

		                disks[i].setVelx(BELT_vel_x_bot);
		            }

		            if (disks[i].wall_tag == 2) // this is redundant if only one wall is moving
		                disks[i].setVely(BELT_vel_y_right);

		            if (disks[i].wall_tag == 3) // this is redundant if only one wall is moving
		                disks[i].setVelx(BELT_vel_x_top);

		            if (disks[i].wall_tag == 4) // this is redundant if only one wall is moving
		                disks[i].setVely(BELT_vel_y_left);

		        }
//		        cout << "Belt velocity = " << BELT_vel_x_bot << endl;
            	}
*/
            render(j);
         }
         
	 cout << "\nTotal time elapsed = " << (tot_simul_time/3600) << " h" << endl;
    }
};



int main (int argc,char *argv[])
{
    bool interiorParticlesProvided = 0, wallParticlesProvided = 0;
//    short unsigned interiorPartType = 1;
    if( argc > 1 )
    {
      t_initial = atof(argv[1]); //atof converts string to floating point. <stdlib.h>
      cout << "\n ***** Initial time = " << t_initial << " s"<< endl;
      t_final = atof(argv[2]);
      cout << "\n ***** Simulation until time = " << t_final << " s" << endl;
      interiorParticlesProvided = atof(argv[3]);
      wallParticlesProvided = atof(argv[4]);

//      RAD = atof(argv[5]); // Not valid if initial interior particles' file is provided
//      CONC = atof(argv[6]); // Not valid if initial interior particles' file is provided
//      Secondary_RAD = atof(argv[7]); // Not valid if initial interior particles' file is provided

//      interiorPartType = atof(argv[8]); // Not valid if initial interior particles' file is provided. Default value of 1 corresponds to mono-disperse or concentration based distribution
//      gaussStdDevPercRAD = atof(argv[9]); // Not valid if initial interior particles' file is provided

//      SMALL_RAD = atof(argv[10]); // Not valid if initial interior and wall particles' file is provided
//      BIG_RAD = atof(argv[11]); // Not valid if initial interior and wall particles' file is provided

//      BELT_vel_x_bot = atof(argv[12]);
//      BELT_vel_y_left = atof(argv[13]);
//      BELT_vel_x_top = atof(argv[14]);
//      BELT_vel_y_right = atof(argv[15]);
//      time_BELT_starts = atof(argv[16]);
//      cout << "\n ***** Belt starts at = " << time_BELT_starts << " s" << endl;

      Cv = atof(argv[5]);
      time_SYS_starts = atof(argv[6]);
      time_HUNT_starts = atof(argv[7]);
      hunt_timelimit = ceil(atof(argv[8])/dt);
      refocus_timelimit = ceil(atof(argv[9])/dt);
      satisfy_timelimit = ceil(atof(argv[10])/dt);
      blind_angle = M_PI*atof(argv[11])/360;	// Converted to radians and divided by 2 as we need value of half of the blind angle for implementation
      noise = atof(argv[12]);
      
      if (time_SYS_starts < t_final)	cout << "\n ***** SPP starts at t = " << time_SYS_starts << " s" << endl;
      if (time_HUNT_starts < t_final)
      {	
	      cout << "\n ***** Hunt starts at t = " << time_HUNT_starts << " s" << endl;
	      cout << "\n ***** Time limit for hunt = " << hunt_timelimit*dt << " s" << endl;
	      cout << "\n ***** Time limit for refocus = " << refocus_timelimit*dt << " s" << endl;
	      cout << "\n ***** Time limit for satisfaction = " << satisfy_timelimit*dt << " s" << endl;
      }
      cout << "\n ***** Blind angle = " << blind_angle*2.0 << " rad" << endl;
      cout << "\n ***** Variance of zero-mean GWN = " << (noise*noise/12.0) << endl;
    }
//	cout << "\n Observor particle ID = " << obs_p_id << endl;

    srand(static_cast<unsigned int>(time(nullptr))); //Seed selection of random function set to clock time of system

    // Generating look-up table
    for (unsigned int j = 0; j < n_ex; j++)
    {
        exp_power = static_cast<float> (j)/n_ex;
        exf[j] = exp(-eta * exp_power * exp_power);
    }

    simulate start = simulate();

    if (wallParticlesProvided)        start.boundaryParticlesRead();
//    else        start.boundaryParticlesGenerate();


    if (interiorParticlesProvided)    start.interiorParticlesRead();
//    else        start.interiorParticlesGenerate(interiorPartType);

    start.mainLoop();

//    ofp_energy.close();
//    ofp_monitor.close();
//	ofp_obs.close();
    return 0;
}
