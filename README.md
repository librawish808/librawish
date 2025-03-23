Visual Studio 2022, Windows 11 x64
Project Name: konkuk

This project implements a simple ray tracer using C++ along with OpenGL/GLUT. The final output shows white wherever a ray intersects an object and black for the background. The ray tracer uses intersection logic for both spheres and a plane.

x32 ⇒ x64 Transition Issue:
Initially, I encountered an LNK4272 error when using freeglut32.lib and similar libraries, but I resolved the issue by switching to the 64-bit versions.

freeglut: 3.6.0#1 (x64-windows)

glew: 2.2.0#6 (x64-windows)

glfw3: 3.4#1 (x64-windows)
(Installed via vcpkg)
(Also, Git was used for version control.)

main.cpp Overview:

The render window is set to 512×512 pixels.

I created classes: Ray, Camera, Surface, Plane, Sphere, and Scene.

In the render function, I reversed the loop order (iterating rows from Height-1 down to 0) so that the image displays correctly (i.e., not flipped) when first built and debugged.

For each pixel, if a ray intersects any object, the pixel is set to white; otherwise, it remains black.

I did not modify the resize callback function.

In main(), the camera is positioned at the origin with the standard basis vectors for u, v, and w. I added three spheres to the scene, each with defined center coordinates, radius, and color.

I set up a plane with the equation y = -2 as required, but since it made half the output image white, I commented it out.
