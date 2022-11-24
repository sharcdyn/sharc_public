1.- Compile:
	Directory src --> make COMP=gfortran MKL=no STATIC=no
	Directory src --> make COMP=gfortran MKL=no STATIC=no MODEL=IBr
	Put the src directory in the PATH --> export PATH=$PWD:$PATH

2.- Create initial conditions --> init
	File geom at equilibrium --> xyz2.exe
	Obtain Force Constant Matrix --> model_IBr.exe
	Create the normal modes --> sharc_init0.exe < init0.in | tee init0.out
	Sampling normal mode 6 -->  sharc_init1.exe < init1.in | tee init1.out
	Mix normal modes -->        sharc_init2.exe < init2.in | tee init2.out

3.- Create a semiclassical spectra --> spec
	Copy files and geoms.init from init
	Create spectra selecting trajs between 2 and 3 --> sharc_spectra.exe < spec.in | tee spec.out
	Analyze the spectrum.dat --> gnuplot
		plot "spectrum.dat" u ($1*27.211):2 w l

4.- Run trajectories in the queue --> TRAJs
	Copy the file geoms from spec
	Run one trajectory --> sbatch --array=1 run_traj.sh
	Analyze the trajectory --> gnuplot
		plot "energy.out" u 1:($2+$3) w l lw 3,"" u 1:3 w p lt 1, for [i=1:3] "" u 1:(column(i+3)) w l lt 3
		plot [][-.01:1.01] "pop_a.out" u 1:2 w l,"" u 1:3 w l,"" u 1:4 w l
		plot [][-.01:1.01] "pop_raw.out" u 1:2 w l,"" u 1:3 w l,"" u 1:4 w l

5.- Run trajectories in the queue --> TRAJs_model
	Copy the file geoms from spec
	Run all trajectories --> sbatch --array=1-100 run_traj.sh
	Analyze all trajectories --> python3 ../../analysis/do_norm.py 1 100
	Go to the directory steps to check the results --> gnuplot
		p "norm.out" u 1:3 w l,"" u 1:4 w l,"" u 1:5 w l
		p "norm.out" u 1:($3+$4) w l,"" u 1:5 w l
		p "norm_raw.out" u 1:2 w l,"" u 1:3 w l,"" u 1:4 w l
		p "norm_raw.out" u 1:2 w l,"" u 1:3 w l,"" u 1:4 w l,"norm.out" u 1:($3+$4) w l,"" u 1:5 w l
