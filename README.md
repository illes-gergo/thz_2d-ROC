# THz 2D project
This project is made for simulating the generation THz pulses by solving a set off coupled differential equations in the Julia programming language.
## Requirements
In order to use this software one must have a ROCm supported graphics card.
## Usage
The simulation can be started by including the ```main.jl``` script, then running the ```runcalc()``` function. Free paramteres can be changed in ```valtozok.jl```, new materials can be implemented in ```fuggvenyek.jl```. In order to implement new materials one has to add functions for refractive index, nonlinear recraftive index, and the cross sections for THz and second-harmonic generation.
