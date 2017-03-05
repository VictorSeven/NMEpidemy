NMEpidemy
=========

Simulation code for a Non-Markovian epidemic spread using the SIS model. This code was done originally for my work on Stochastic Simulation Methods, for the Master in Physics of Complex Systems.

Most common epidemic simulations are Markovian (i.e. they have constant infection/recover rates). In this case, every infected node can track the time it has been infected, and its rate depends on this individual time. The algorithm implemented is a Non-Markovian Gillespie, as in the reference Phys. Rev. E 90, 042108 (Arxiv avaiable [here](https://arxiv.org/abs/1310.0926)). 

The simulation runs over a network (created with [CNetwork](https://github.com/VictorSeven/CNetwork)) and computes the percentage of infected people for several values of the effective-infection rate. It can be seen that the Non-Markovian case really differs from the Markovian one. Simulations were done for an Erdos-Renyi network with average degree 5.

I used this program to see the change of the epidemic threshold when the rates follow a Weibull distribution. Results are plotted below. X axis of the graph is the effective infection rate, while Y is the percentage of infected individuals. For alpha = 0, this is the classical Markovian SIS model.

![alt text](https://github.com/VictorSeven/NMEpidemy/blob/master/source/weib.png "Results")

You can use this code at you own. Feel free to contact for more information about it.


