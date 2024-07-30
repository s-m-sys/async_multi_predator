#############----------------------#############
#	This is the input file for the code.
#	It does the following:
#		a) Compiles the code
#		b) Creates a directory particle_data. If the directory is already present, it leaves it as it is.
#		c) Executes the program with the inputs specified below
#	WARNING: Please ensure there are no spaces before/after the equal to (=) sign in the following statements
#############----------------------#############
# Compiling and generating the executable (machine level instruction)
g++ SPP_active_vision_periodic.cpp -std=c++14 -O3 -fconcepts

mkdir particle_data

# There are 24 arguments along with the executable: 
t_initial=0.0 #1 Initial time = 0.0 (default)
t_final=1000.0 #2 Final time = 410.0 (default)

# Specify whether initial data file present or not. If then, the files should be named "interior_ini.txt" and "wall_ini.txt", respectively, and kept in the same directory as the executable 'a.out'

IntDataFile_isThere=1 #3 Initial interior particles data file present or not (Yes = 1, No = 0) = 0 (default)
WallDataFile_isThere=0 #4 Initial wall particles data file present or not (Yes = 1, No = 0) = 0 (default) # If provided, the walls must be the 4 sides of a rectangular enclosure. The bottom and top walls can only have x-velocity while the side walls can have only the y-component of velocity

########### Ignore these inputs if the initial interior and wall particles files are provided ################

#R_m=0.0025 #5 Mean radius of the particles (m) = 0.0025 (default)
#N_conc=0.0 #6 Number concentration of secondary particles = 0.0
#r_m=0.0 #7 Mean radius of the small (secondary, impurity) particles (m) = 0.00125 (default) 

#interiorParticleType=1 #8 Type of poly-dispersity for the particles in the domain. 
#1 --> Mono/bi-disperse (CONC. based)
#2 --> Uniform random within r_small and r_big
#3 --> Gaussian random distribution, with standard deviation percentage below
#gauss_stdDev_percRm=0.01 #9 Standard deviation of the Gaussian distribution for radius (in terms of % mean radius, R_m). Relevant only if Gaussian distribution is active

#r_small=0.0025 #10 Minimum radius of the particles in the domain
#r_big=0.0025 #11 Maximum radius of the particles in the domain

#vel_BOT_wall=0.0 #12 Bottom belt's x-velocity (m/s) = 0.0 (default)
#vel_RIGHT_wall=0.0 #13 Right wall's y-velocity (m/s) = 0.0 (default)
#vel_TOP_wall=0.0 #14 Top wall's x-velocity (m/s) = 0.0 (default)
#vel_LEFT_wall=0.0 #15 Left wall's y-velocity (m/s) = 0.0 (default)
#time_BELT_starts=1000.00 #16 Belt starts at time (s) = 100.0 (default)

########### *************************************************************** ################

Cv=45 #17 CV = 0.0 (default)
time_SPP_starts=0.001 #18 Self-propelling particles get active at time (s) = 0.0 (default)
time_HUNT_starts=0.002  #19 Hunting starts at time (s)= 1000.0 (default)
hunt_time=20.0	#20 Maximum time for a single chase
refocus_time=5.0	#21 Time taken by predator to recover after a failed hunt
satisfy_time=10.0	#22 Time taken by predator to start hunting after a successful hunt
blind_angle=60.0	#23 Blind angle of each particle (in degrees)
noise_mag=1.0	#24 Magnitude of std dev of white noise = 0.0 (default)

./a.out $t_initial $t_final $IntDataFile_isThere $WallDataFile_isThere $Cv $time_SPP_starts $time_HUNT_starts $hunt_time $refocus_time $satisfy_time $blind_angle $noise_mag >&1 | tee logfile.txt
###################Following part is for graph and graphical results purpose, not needed at the moment
#gnuplot Post_eg.gp
#gnuplot Post_graph.gp

#for f in *.eps; do j=${f%.*}; convert -density 200 $j.eps $j.png; done

#int_part=`cat ./particle_data/interior_0.05.txt| wc -l`
#tot_files=`ls particle_data/interior* -1|wc -l`
#sed -i 's/#define core.*/#define core '${int_part}'/' ./my_svd.c
#sed -i 's/#define total_files.*/#define total_files '${tot_files}'/' ./my_svd.c

#gcc my_svd.c -lm -llapack -O2 -o svd

#./svd

