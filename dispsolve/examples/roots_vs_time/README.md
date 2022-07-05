Here we attempt to solve the dispersion relation in order to obtain roots (values of 
complex omega which solve the dispersion relation) as a function of time in the simulation,
for a given k0 of interest. In the simulation which generated the data found in `some_data`,
nonlinear effects result in a significant modification to the distribution function. Here 
we attempt to observe the effect of this on the damping of plasma waves at a particular k.

To observe this, we first plot the distribution function at three different times throughout
the simulation (f_vs_time.png). Notice how at the phase velocity corresponding to our chosen
k, the distribution function develops a pertubation.

Next, we observe how this pertubation affects the damping of plasma waves at the chosen k
(roots_vs_time.png). We see that initially (pre-pertubation), damping is zero. Then, damping
begins to increase and oscilate about some value as a result of the (oscilating) pertubation.
