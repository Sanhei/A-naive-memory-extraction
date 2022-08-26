# Markovian Embedding model in Double-well potential.

## Trajectory making
This will in the document named (timestep_changing)
We put a particle in a double well potenital(potential.cpp), 
![alt text](https://github.com/Sanhei/A-naive-memory-extraction/blob/main/potential.png?raw=true)

using the Markovian embedding system, which the discretizing equation is
$$x(t_{j+1})=2x(t_j)-\frac{\tau_D}{\tau_m}F(x(t))t^2-\frac{\tau_D^2}{\tau_m\gamma}(y(t_j)-x(t_j))\Delta t^2.$$
$$y(t_{j+1})=y(t_j)-\tau_D/\gamma(y(t_j)-x(t_j))\Delta t+\sqrt{2\Delta t}\xi$$

which in the document parameter.h saved all the parameter we need, and the integrator function is stored in the v_verlet.cpp. And in the traj.cpp we just do the iteration.

And it will generate a .txt file which record each step of time position. And a single example:
![alt text](https://github.com/Sanhei/A-naive-memory-extraction/blob/main/0.005_traj.svg?raw=true)
## Trajectory analyse.
We need to extract the memory kernel, which is $$\gamma = \int^\infty_0\Gamma(t)$$.
To calculate this part, we need to get the correlation function of $$C_{vv}$$ and $$C_{\triangledown U x}$$.
### $$\triangledown U$$
For calculating $$\triangledown U$$, we use free energy, which get from the distribution of the position. Therefore, potential.h gives three function.
1. Histogram calculation, which will record the distribution of position in the trajectory. And $$F = -k_BTln(P(x))$$, so this can directly transfer to potential.
2. Gradient of potential.
3. Calculate a small area around boundary, because the behavior is not like stochastic behavior. So will print out how much influnce on the whole calculation. And a corretion, which may a little bit shift due numerical error.



### Correlation calculation
This part, we do correlation function, first we need to calculate the velocity. After that we get everything to calculate the correlation function. For calculating the correlation, we will do three steps for that:
1. Do the Fast Fourier Tranfer(FFT) to the variables. Here in our case, they are velocity, position and derivertive of potential.
2. Get the convolution product of the variables(which got from FFT).
3. Reverse FFT to the convolution product, then we will get correlation function.

## Memory kernel Calculation.
Here we use the equation from paper(Non-markov protien folding).
$$G_n = \frac{2}{\Delta t C^{vv}_0} (C^{\triangledown U q}_n -\frac{C^{\triangledown U q}}{C^{vv}_0} - \Delta t\sum^{n-1}_{i=1}G_{n-i}C^{vv}_i)$$.
Then all the calculation done.

## Program explain.
In our program, we can call matplotlib function in python to directly plot the Potential and Correlaiton. And in the end, after the calculation, program can automatically call the fitting.py, which uses lmfit package to do the fitting. Here we may need the python package path to include this environment.
## Runing the program.
```bash
export PYTHONPATH="~/python_path/lib/python3.x/site-packages"
```
This is to include the path of python, but this step may not necessary, since we can comment all the code of python, the simulation is all done by C++ code.

```bash
cd timestep_changing
g++ -o traj traj.cpp
./traj
```
This part is for the simulation, will go the file, only to compile the code.
```bash
cd ..
make
./memoryextraction
```
Then the data and all figures will in the "figure" director.
## Space and Memory
For now, we use 2GB trajectory, for several parameter settings, it works. However it is not enough, so maybe we want to try 8GB trajectory. As for memory occupied, FFT may use three times bigger than the trajectory size. So for each time correlation function, we record it in a text file, and read it again.
The thing we need to care about is $$C_{\triangledown U x}$$, We need to record an additional trajectory size. So in total, it will $$\mathbf{8}$$ times bigger memory. For the optimization, we just use one core, and didn't set the threads. For FFTW, it should be defined as 4.
