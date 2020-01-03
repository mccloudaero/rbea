#!/usr/bin/env python
import os
import sys
import math
import numpy as np
import scipy.interpolate
import matplotlib.pyplot as plt
# Rotor Blade Element Analysis
# Python Script for determining Rotor Performance using Blade Element Theory

# Definitions
# theta - The physical angle the airfoil is rotated from the plane of rotation
# phi	- The change in flow angle due to the induced airflow.  Makes the angle of attack smaller
# alpha_rad	- The effective angle of attack the airfoil experiences
# alpha_rad + phi = theta

print('Rotor Blade Element Analysis\nMcCloud Aero Corp.')

pi = math.pi
pitch_angles = [0]

# Get dir where script lives
script_path = os.path.dirname(os.path.abspath(__file__))
airfoils_path = os.path.join(script_path,'airfoils')

# Read inputs file
exec(open('rotor.inputs').read())

# Load airfoil data
print('Loading Airfoil Data')
cd_splines = []
cl_splines = []
cm_splines = []
spline_dims = []
num_airfoils = len(airfoil)
for airfoil_name in airfoil:
  print(airfoil_name)
  airfoil_path = os.path.join(airfoils_path,airfoil_name+'.dat')
  airfoil_data = np.loadtxt(airfoil_path,skiprows=1,delimiter=',')
  # Create splines
  mach_data = airfoil_data[:,0]
  alpha_data = airfoil_data[:,1]
  cd_data = airfoil_data[:,2]
  cl_data = airfoil_data[:,3]
  cm_data = airfoil_data[:,4]
  # Determine Mach Range
  if min(mach_data) == max(mach_data):
    # Single Mach Number present, Create 1D splines
    spline_dims.append(1) 
    cd_splines.append(scipy.interpolate.splrep(alpha_data,cd_data))
    cl_splines.append(scipy.interpolate.splrep(alpha_data,cl_data))
    cm_splines.append(scipy.interpolate.splrep(alpha_data,cm_data))
  else:
    # Multiple Mach Numbers present, Create 2D splines
    spline_dims.append(2) 
    cd_splines.append(scipy.interpolate.bisplrep(mach_data,alpha_data,cd_data))
    cl_splines.append(scipy.interpolate.bisplrep(mach_data,alpha_data,cl_data))
    cm_splines.append(scipy.interpolate.bisplrep(mach_data,alpha_data,cm_data))


# Initial Calculations
tip_radius = radii[-1]
diameter = tip_radius*2
area = pi*tip_radius**2
avg_chord = (chord[0] + chord[-1])/2
sigma=2*avg_chord/2/pi/tip_radius			# solidity

# Blade Element Calculation
#r_inc = (radii[-1]-radii[0])/(num_elements-1)
#theta_inc = (tip_theta-root_theta)/(num_elements-1)
#chord_inc = (tip_chord-root_chord)/(num_elements-1)


# Initialize total_output file
total_output = open('./output.txt','w')
# Print Data Header
total_output.write(' RPM,  Pitch, T(N),   P(kW), M(N-m), M(kg-cm), T(lbf), P(Hp), M(ft-lb), FM,    M_tip, v_ind, cd\n')
  
# Initialize total_output file
radial_output = open('./radial_output.txt','w')
# Print Data Header
radial_output.write('i\\RPM\t')
for i in range(num_elements):
  radial_output.write(str(i)+'\t')
radial_output.write('\n')

# Specified Theta
#theta = [18.00, 17.34, 16.68, 16.03, 15.37, 14.71, 14.05, 13.39, 12.74, 12.08, 11.42, 10.76, 10.10, 9.44, 8.79, 8.13, 7.47, 6.82, 6.16, 5.50]
#for i in range(num_elements):
# theta_rad[i] = (theta[i]/180*pi)
#print theta

# Setup Rotor model
f_chord = scipy.interpolate.interp1d(radii, chord)
f_theta = scipy.interpolate.interp1d(radii, theta)

radius_interp = np.linspace(radii[0], radii[-1], num=num_elements, endpoint=True)
chord_interp = f_chord(radius_interp)
theta_interp = f_theta(radius_interp)
theta_rad = theta_interp/180.0*pi

# Plot shape and save
plt.plot(radii, chord, 'o')
plt.plot(radius_interp, chord_interp,'o-',fillstyle='none')
plt.xlabel('Radius')
plt.xlim(left=0)
plt.ylabel('Chord')
plt.ylim(bottom=0)
plt.savefig('./rotor_shape.png')

for pitch_angle in pitch_angles:
  for RPM in RPMs:

    # Initial Calcs at RPM
    n = RPM/60               # RPS
    omega = n*2*pi           # angular velocity
    V_tip = omega*tip_radius # tip velocity

    # Initialize lists
    alpha_rad = []
    for i in range(num_elements):
      alpha_rad.append(0)

    # Initialize/Clear Lists
    thrust = []
    drag = []
    drag_coef = []
    moment = []
    torque = []

    inboard_airfoil_index = 0
    outboard_airfoil_index = 1 

    for i in range(num_elements):
      # Guess at initial values of inflow and swirl factor
      # Note: swirl currently isn't used
      v_inflow = 10
      #v_swirl = 0.1
      current_r = radius_interp[i]
      # increment
      r_inc = radius_interp[1]-radius_interp[0] 
      if ( i == 0 or i == num_elements-1):
        r_inc = r_inc/2.0
      # Check inboard/outboard airfoils
      if(current_r > radii[outboard_airfoil_index]):
        # Advance to next set
        inboard_airfoil_index += 1 
        outboard_airfoil_index += 1 

      # Iterate at Each Blade Element
      iter = 0
      finished = False
      while( finished == False):
        v_axial = v_inflow				# axial velocity
        v_radial = omega*radius_interp[i]			# disk plane velocity
        #v_radial = omega*radius[i] - v_swirl	# disk plane velocity
        v_mag = pow((v_axial*v_axial+v_radial*v_radial),0.5)	# local velocity at blade (Resultant of Induced & Radial)
        mach = v_mag/a
        phi_rad = math.atan2(v_axial,v_radial)	# flow angle (radians)
        alpha_rad[i] = theta_rad[i]- phi_rad		# blade angle of attack
        alpha = alpha_rad[i]*180/pi
        #alpha = min(alpha,12)

        # Find section coefficients
        if i == 0 or num_airfoils == 1:
          # Root or Constant
          dim = spline_dims[0]
          if dim == 1:
            cd = scipy.interpolate.splev(alpha,cd_splines[0])
            cl = scipy.interpolate.splev(alpha,cl_splines[0])
            cm = scipy.interpolate.splev(alpha,cm_splines[0]) 
          elif dim == 2:
            cd = scipy.interpolate.bisplev(mach,alpha,cd_splines[0])
            cl = scipy.interpolate.bisplev(mach,alpha,cl_splines[0])
            cm = scipy.interpolate.bisplev(mach,alpha,cm_splines[0]) 
        elif i == num_elements:
          # Tip 
          dim = spline_dims[-1]
          if dim == 1:
            cd = scipy.interpolate.splev(alpha,cd_splines[-1])
            cl = scipy.interpolate.splev(alpha,cl_splines[-1])
            cm = scipy.interpolate.splev(alpha,cm_splines[-1]) 
          elif dim == 2:
            cd = scipy.interpolate.bisplev(mach,alpha,cd_splines[-1])
            cl = scipy.interpolate.bisplev(mach,alpha,cl_splines[-1])
            cm = scipy.interpolate.bisplev(mach,alpha,cm_splines[-1]) 
        else:
          # Interpolate between two airfoils
          # Find inboard and outboard airfoils
          inboard_dim = spline_dims[inboard_airfoil_index]
          if inboard_dim == 1:
            inboard_cd = scipy.interpolate.splev(alpha,cd_splines[inboard_airfoil_index])
            inboard_cl = scipy.interpolate.splev(alpha,cl_splines[inboard_airfoil_index])
            inboard_cm = scipy.interpolate.splev(alpha,cm_splines[inboard_airfoil_index]) 
          elif inboard_dim == 2:
            inboard_cd = scipy.interpolate.bisplev(mach,alpha,cd_splines[inboard_airfoil_index])
            inboard_cl = scipy.interpolate.bisplev(mach,alpha,cl_splines[inboard_airfoil_index])
            inboard_cm = scipy.interpolate.bisplev(mach,alpha,cm_splines[inboard_airfoil_index])
          outboard_dim = spline_dims[outboard_airfoil_index]
          if outboard_dim == 1:
            outboard_cd = scipy.interpolate.splev(alpha,cd_splines[outboard_airfoil_index])
            outboard_cl = scipy.interpolate.splev(alpha,cl_splines[outboard_airfoil_index])
            outboard_cm = scipy.interpolate.splev(alpha,cm_splines[outboard_airfoil_index]) 
          elif inboard_dim == 2:
            outboard_cd = scipy.interpolate.bisplev(mach,alpha,cd_splines[outboard_airfoil_index])
            outboard_cl = scipy.interpolate.bisplev(mach,alpha,cl_splines[outboard_airfoil_index])
            outboard_cm = scipy.interpolate.bisplev(mach,alpha,cm_splines[outboard_airfoil_index])

          # Interpolate cd, cl, cm
          delta_r = radii[outboard_airfoil_index]-radii[inboard_airfoil_index]
          outboard_portion = (current_r-radii[inboard_airfoil_index])/delta_r 
          inboard_portion = 1.0-outboard_portion 
          cd = outboard_portion*outboard_cd + inboard_portion*inboard_cd
          cl = outboard_portion*outboard_cl + inboard_portion*inboard_cl
          cm = outboard_portion*outboard_cm + inboard_portion*inboard_cm
    
        # Apply drag factor
        cd = cd*drag_factor   
   
        q = 0.5*rho*pow(v_mag,2)                                                # local dynamic pressure
        DtDr = q*num_blades*chord_interp[i]*(cl*math.cos(phi_rad)-cd*math.sin(phi_rad))    # thrust grading
        DdDr = q*num_blades*chord_interp[i]*(cd*math.cos(phi_rad)+cl*math.sin(phi_rad))    # drag grading
        DmDr = q*num_blades*chord_interp[i]*cm                                             # moment grading
        DqDr = q*num_blades*chord_interp[i]*radius_interp[i]*(cd*math.cos(phi_rad)+cl*math.sin(phi_rad)) # torque grading
   
        # momentum check on inflow and swirl factors
        v_inflow_new = DtDr/(4*pi*radius_interp[i]*rho*v_axial)
        #v_swirl_new = DqDr/(4*pi*pow(radius[i],3)*rho*v_axial*omega)
           
        # increment iteration count
        iter += 1
  
        # check for convergence
        if ( math.fabs(v_inflow_new-v_inflow)< 0.001):
          finished = True
        # check to see if iterations are stuck
        elif(iter>maximum_iterations):
          finished=True
          print('RPM:',RPM,'element:',i,'exceed maximum number of iterations (',maximum_iterations,')')
  
        # Updates Values
        v_inflow = v_inflow + relaxation_factor*(v_inflow_new-v_inflow)
        #v_swirl = v_inflow + 0.5*(v_swirl_new-v_swirl)    
  
      phi = phi_rad*180/pi
      alpha = alpha_rad[i]*180/pi 
      thrust.append(DtDr*r_inc)
      drag.append(DdDr*r_inc)
      drag_coef.append(cd)
      moment.append(DmDr*r_inc)
      torque.append(DqDr*r_inc)

      #print RPM,i,v_inflow

    # Convert radians to degrees
    radial_output.write('{:d}\t'.format(RPM))
    for i in range(num_elements):
      #theta[i] = theta[i]*180/pi
      alpha_rad[i] = alpha_rad[i]*180/pi
      radial_output.write('{:5.3f}\t'.format(alpha_rad[i]))
    radial_output.write('\n')
 
    # Totals
    T = sum(thrust)                           # total thrust (N)
    P = sum(torque)*omega/1000                # total power (kW)
    M = sum(moment)                           # total moment (N-m) 
    M_kg = M*10.1971621298                    # total moment (kg-cm)
    drag_coef_avg = sum(drag_coef)/num_elements  # total thrust (N)

    # Ideals
    v_ideal = pow((T/(2*rho*pi*pow(tip_radius,2))),0.5) # Ideal Induced Velocity (m/s)
    P_ideal = T*v_ideal/1000 # Ideal Power (kW)
 
    # Compute Coefficients
    ct = T/(rho*n*n*pow(diameter,4))
    cq = sum(torque)/(rho*n*n*pow(diameter,5))

    # Adjust for Tip Loss using Empirical Equations
    B = 1 - pow(2*ct,2)/2
    P = P_ideal/B + (P - P_ideal)
 
    T_imp = T*0.224808943			# N to lbf
    P_imp = P*0.00134102209*1000	# kW to Hp
    M_imp = M*0.737562149			# N-m to ft-lbf
 
    FM = P_ideal/P
    M_tip = v_mag/a
    total_output.write('{:5d}, {: 3.2f}, {:5.2f}, {:5.2f}, {: 6.3f}, {: 8.2f}, {:6.2f}, {:5.2f}, {:8.3f}, {:5.3f}, {:5.3f}, {:5.3f}, {:5.4}\n'.format(RPM,pitch_angle,T,P,M,M_kg,T_imp,P_imp,M_imp,FM,M_tip,v_ideal,drag_coef_avg))

# Close total_output
total_output.close()
radial_output.close()
