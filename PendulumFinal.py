#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 25 20:07:27 2019

@author: ginapantano
"""

"""
Final Project: The Double Pendulum
By: Gina Pantano, Miliani Hernandez, Kelli Shar

The double pendulum is a system where two masses are connected by rigid massless
rods. The double pendulum is a great example of a simple physical system which
can become quickly chaotic.   
"""
#Imported Packages
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
import matplotlib.animation as animation

def Pendulum(Y,t,m1,m2,l1,l2,g):
    '''
    The function Pendulum solves for the change in theta1, theta2, omega1, and
    omega 2 at every time step starting with the intiial conditions at t = 0
    for each pendulum mass m1 and m2.
    '''
    
    theta1, dot_theta1, theta2, dot_theta2 = Y
    
    cos = np.cos(theta1-theta2)
    sin = np.sin(theta1-theta2)
    
    dtheta1 = dot_theta1
    ddot_theta1 = ((m2*g*np.sin(theta2)*cos) - (m2*sin*(l1*dot_theta1**2*cos + l2*dot_theta2**2)) 
    -((m1+m2)*g*np.sin(theta1))) / (l1*(m1 + m2*sin**2))
    dtheta2 = dot_theta2
    ddot_theta2 = (((m1+m2)*(l1*dot_theta1**2*sin - g*np.sin(theta2) + g*np.sin(theta1)*cos)) 
    + (m2*l2*dot_theta2**2*sin*cos)) / (l2*(m1 + m2*sin**2))
    
    return [dtheta1, ddot_theta1, dtheta2, ddot_theta2]

def Energy_Function(th1, dot_th1, th2, dot_th2):
    '''
    The function Energy_Function takes in theta1, theta2, omega1, and omega2 at 
    every time step starting with our initial conditions at t = 0. This function
    calculates the kinetic and potential energy, and returns the total energy
    E = T + U.
    '''
    
    cos = np.cos(th1-th2)
    
    T = (.5*(m1+m2)*l1**2*dot_th1**2) + (.5*m2*l2**2*dot_th2**2) + (m2*l1*l2*dot_th1*dot_th2*cos)
    U = -((m1+m2)*g*l1*np.cos(th1)) - (m2*g*l2*np.cos(th2))
    
    return abs(T + U)

#Constants
g = 9.81
m1 = 5     #m1 and m2 in kilograms
m2 = 5      
l1 = 2     #l1 and l2 in meters
l2 = 2

#Time and Initial Conditions
time,dt = np.linspace(0,20,1000, retstep = 'True')   #Time range and time steps 
init = [.2, 3, 0, 2]    #[theta1, dtheta1, theta2, dtheta2]

#Solution Array
soln = odeint(Pendulum,init,time,args=(m1,m2,l1,l2,g,))

#Variable Solution
theta1, dot_theta1, theta2, dot_theta2 = soln[:,0], soln[:,1], soln[:,2], soln[:,3]

#Energy Calculations - Conservation of Energy
Total_Energy = Energy_Function(theta1, dot_theta1, theta2, dot_theta2)
Initial_Energy = Energy_Function(init[0],init[1],init[2],init[3])

#Checks to see if the Total_Energy - Initial_Energy never exceeds 1e-5
#In other words, checks to see if the energy is conserved over all time.
for i in range(len(Total_Energy)):
    '''
    Will print "Energy is not conserved" at a given time step t if the
    change in energy exceeds 1e-5
    '''
    if abs(Total_Energy[i]-Initial_Energy)>1e-5:
        print("Energy is not conserved at time step {}".format(i))
        break
    else:
        pass

#Conversion to Cartesian coordinates of the two pendulum masses.
x1 = l1 * np.sin(theta1)
y1 = -l1 * np.cos(theta1)
x2 = x1 + l2 * np.sin(theta2)
y2 = y1 - l2 * np.cos(theta2)

delta_theta = abs(theta1-theta2)

#Plotting Commands
fig = plt.figure(figsize = (10,10))

#-----------------------------------------------------------------------------
"""
Source: https://matplotlib.org/3.1.1/gallery/animation/double_pendulum_sgskip.html

    This source helped us visualize our data through an animation of our
    two pendulum masses. We wanted to be able to present our data in the 
    clearest way possible. This animation tracks the value of our variables at 
    every time step from 0 to t. 
"""
ax = fig.add_subplot(311, autoscale_on= False, xlim=(-6, 6), ylim=(-6, 0.2))
ax.set_aspect('equal')
ax.grid()

line, = ax.plot([], [], 'co-', lw=1)
time_template = 'time = %.1fs'
time_text = ax.text(0.5, 0.9, '', transform=ax.transAxes)

def init():
    line.set_data([], [])
    time_text.set_text('')
    return line, time_text


def animate(i):
    thisx = [0, x1[i], x2[i]]
    thisy = [0, y1[i], y2[i]]

    line.set_data(thisx, thisy)
    time_text.set_text(time_template % (i*dt))
    return line, time_text

ani = animation.FuncAnimation(fig, animate, range(1, len(soln)),
                              interval=dt*1000, blit=True, init_func=init)

#-----------------------------------------------------------------------------
#Our Plotting Commands

#Cartesian Coordinates
plt.subplot(3,1,2) 
plt.plot(x1,y1,'c', label = 'Top Mass')
plt.plot(x2,y2, label = 'Bottom Mass')
plt.ylabel('y (m)')
plt.xlabel('x (m)')
plt.title('Motion of the Double Pendulum')
plt.legend()

#Polar Coordinates
plt.subplot(3,1,3)
plt.plot(theta1,'c', label = 'Top Angle')
plt.plot(theta2, label = 'Bottom Angle')
plt.ylabel('$\\theta_1$ and $\\theta_2$(radians)')
plt.xlabel('$Time (s)$')
plt.legend()

#Change in Theta over Time
'''
plt.plot(time,delta_theta)
plt.ylabel('$\Delta \\theta$')
plt.xlabel('Time (s)')
plt.title('Angular Frequency')
plt.xlim(0,8)
'''

'''
Note: We quoted this graph out since the most intersting result is when the 
system is chaotic. We have included our images within the LaTeX file portion.
'''
plt.show()



















