#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug  4 22:32:53 2023

@author: espoukeh
"""

import numpy as np


def create_plane(normal_direction, l):
    if normal_direction == 'x':
        P0 = [0, -l/2, -l/2]
        P1 = [0, -l/2, l/2]
        P2 = [0, l/2, -l/2]
        
        return P0, P1, P2
    elif normal_direction == 'y':
        P0 = [-l/2, 0, -l/2]
        P1 = [-l/2, 0, l/2]
        P2 = [l/2, 0, -l/2]
        return P0, P1, P2
    elif normal_direction == 'z':
        P0 = [-l/2, -l/2, 0]
        P1 = [l/2, -l/2, 0]
        P2 = [-l/2, l/2, 0]
        return P0, P1, P2
    
    elif normal_direction == '-z':
        P0 = [l/2, -l/2, 0]
        P1 = [-l/2, -l/2, 0]
        P2 = [l/2, l/2, 0]
        return P0, P1, P2
    
    else:
        print("Invalid input for Normal_direction. Please choose 'x', 'y', or 'z'.")
        return None

        
def rotation_matrix(axis, angle):
    axis = np.asarray(axis)
    axis = axis / np.sqrt(np.dot(axis, axis))
    a = np.cos(angle / 2.0)
    b, c, d = axis * np.sin(angle / 2.0)
    return np.array([[a*a + b*b - c*c - d*d, 2*(b*c - a*d), 2*(b*d + a*c)],
                     [2*(b*c + a*d), a*a + c*c - b*b - d*d, 2*(c*d - a*b)],
                     [2*(b*d - a*c), 2*(c*d + a*b), a*a + d*d - b*b - c*c]])


def main():
    l = float(input("Enter magnitude l: "))
    normal_direction = input("Enter normal_direction ('x', 'y', or 'z'): ")

    P0, P1, P2 = create_plane(normal_direction, l)

    if P0 is not None:
       print("P0:", P0)
       print("P1:", P1)
       print("P2:", P2)
       
       
       # Get the axis from user as three coordinates (x, y, z) in a single line
       axis = list(map(float, input("Enter the coordinates of the axis rotation(x, y, z): ").split(",")))
       angle_degrees = float(input("Enter the angle of rotation in degrees: "))
       angle_radians = np.radians(angle_degrees)
       

       # Rotate the vertices P0, P1, and P2 around the specified axis
       rotation_matrix_ = rotation_matrix(axis, angle_radians)
 
       
       P0_rotated = np.dot(rotation_matrix_, P0)
       P1_rotated = np.dot(rotation_matrix_, P1)
       P2_rotated = np.dot(rotation_matrix_, P2)

       # Shift the result to the new origin (x1, x2, x3)1
       x1, x2, x3 = map(float, input("Enter the coordinates of the new origin (x1, x2, x3): ").split(","))

       P0_rotated_shifted = P0_rotated + np.array([x1, x2, x3])
       P1_rotated_shifted = P1_rotated + np.array([x1, x2, x3])
       P2_rotated_shifted = P2_rotated + np.array([x1, x2, x3])


       print("Array: [%5.3f, %5.3f, %5.3f]" % (P0_rotated_shifted[0], P0_rotated_shifted[1], P0_rotated_shifted[2]))
       print("Array: [%5.3f, %5.3f, %5.3f]" % (P1_rotated_shifted[0], P1_rotated_shifted[1], P1_rotated_shifted[2]))
       print("Array: [%5.3f, %5.3f, %5.3f]" % (P2_rotated_shifted[0], P2_rotated_shifted[1], P2_rotated_shifted[2]))
   
    else:
       print("Invalid input for Normal_direction. Please choose 'x', 'y', or 'z'.")
    

if __name__ == "__main__":
    main()