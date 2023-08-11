microstructure
==============

Generating the micro-structures is a non-trival problem because there are many 
requirements such as the partiles should not overlap between each other, 
different shapes of partcles should provide, the particles also shoulf obey the periodical 
boundary conditions, and so on.
With regular shapes, such as disk, 
ellipse in 2D scenario, sphere in 3D scenario; we can obtain the analytical 
formula to describe their shapes. Therefore, it provides much more information 
to design an algorithm to generate the micro-structures. 
Moreover, most of those algorithms are based on techniques like Monte Carlo Simulation, 
Melocular dynamics and so forth, so there are no same micro-structures from two different realizations.
What's more, if the shapes are arbitrary then the generation of micro-structure becomes 
very difficult as an NP hard problem.

--selected methods in this version of rve-simulator

1 Melro, A. R., Camanho, P. P., & Pinho, S. T. (2008). Generation of random distribution of fibres in long-fibre reinforced composites. Composites Science and Technology, 68(9), 2092-2102.


2d micro-structure generation
-----------------------------
Here gives an example of 2D micro-structure generation.

.. code-block:: python

    from rvesimulator.microstructure.circle_particles import CircleParticles
    microstructure_generator_2d = CircleParticles(length=1.0,
                                              width=1.0,
                                              radius_mu=0.050,
                                              radius_std=0.001,
                                              vol_req=0.30)
    # generate by specific seed
    microstructure_generator_2d.generate_microstructure(seed=3)
    # plot
    microstructure_generator_2d.plot_microstructure()

.. image:: /images/mircostructure_show.png
   :alt: mircostructure_show
   :width: 50%
   :align: center

.. code-block:: python

    microstructure_generator_2d.crate_rgmsh(num_discrete=600)
    microstructure_generator_2d.rgmsh_plot(save_fig=True, fig_name="rgmsh.png")

.. image:: /images/rgmsh.png
   :alt: rgmsh
   :width: 50%
   :align: center

