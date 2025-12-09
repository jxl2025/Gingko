[![Review Assignment Due Date](https://classroom.github.com/assets/deadline-readme-button-22041afd0340ce965d47ae6ef1cefeee28c7c493a6346c4f15d667ab976d596c.svg)](https://classroom.github.com/a/UOBObzOT)
# Lab 5: Light

Welcome to Lab 5! The handout for the lab is located [here](https://browncsci1230.github.io/labs/lab5).


# Task 2 
# (a) How many times do you evaluate to get one RGB color
Three times, once per channel R, G, and B (since colors/light/materials are RGB triplets).

# (b) How many distinct scalar material parameters are needed?
Five in total. They are ka, kd, ks, kr (weights) and the shininess exponent s. I note that the material base colors (ambient/diffuse/specular colors) are RGB, not scalar, and that theattenuation numbers c1, c2, c3 are for the global light/scene, not material.

this