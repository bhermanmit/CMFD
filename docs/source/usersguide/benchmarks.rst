.. _uersguide_bencmarks:

=====================
Benchmarks using CMFD
=====================

Explain purpose of doing benchmarks

.. note:: 
   Specifications and descriptions of the benchmarks were taken directly out of 
   the following two locations:

   - MULLER, E.Z. and Weiss Z.J., "Benchmarking with the Multigroup Diffusion 
     High-Order Response Matrix Method," Annals of Nuclear Energy, Vol. 18, 
     No. 9, pp. 534-544, 1991.
   - SMITH, K.S. "TITLE OF MASTERS"

------
BIBLIS
------

The 2-D BIBLIS problem is a realistic and highly nonseparable 2-group problem 
representative of an actual operating pressurized water reactor (PWR), with a 
checker-board-loaded core.   Homogenized fuel assemblies with widths of 23.1226
cm and seven diffrence compositions are present in the core which is surrouded by a 23.1226 cm homogenized reflector (i.e. the baffle is homogenized with the water reflector) with vacuum external boundary conditions.

.. important::
   The following changes were made from the original benchmark, but should not 
   affect the answers significantly

   - reflecting material place outside of the core so that the computational 
     domain is a rectangle. Reflector boundary not treated explicitly as in 
     orginal specficiations. See references for true core layout.
   - since the reflecting region is large, zero-incoming partial current
     boundary conditions were used instead of vacuum 
