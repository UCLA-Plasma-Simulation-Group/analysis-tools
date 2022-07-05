Here we attempt to solve the dispersion relation in order to obtain roots (values of 
complex omega) as a function of k, for a maxwellian distribution function. Run the script
and have a look at the output generated in each section.

First, we get the data, normalize it (st area under the curve is 1) and then plot f and
df/dv, and compare with theory as a sanity check (f_dfdv.png)

Next we plot the value of the dielectric for a particular k0 in the chose window of complex
omega (dielectric.png). You can see there are clearly some issues with some particular values of omega,
but the solver does still seem to consistently do a good job of finding the correct solution, 
as we'll see in the next section. To mitigate some numerical difficulty, I found it good
to limit the window in complex space that you allow the solver to look in, but just be careful
not to limit it too much so that you don't rule out any unexpected but valid solutions.

Finally, we plot the solutions as a function of k (roots_vs_k0.png). For the real roots,
we compare with theory. At the cost of more time, we could have also compared the imaginary
roots with Landau's analytical solution. We leave that as an exercise to the reader.
Notice that the agreement is good for low kld, and than begins to deviate at higher kld
as kinetic effects become important.
