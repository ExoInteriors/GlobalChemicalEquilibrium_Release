.. _solver:

Solvers
=======


In this section, we describe first the general solver strategy, how the chemical network equations are solved. Followed by the available optimizers.


General solver strategy
-----------------------

The chemical network equations are solved by using numerical optimization methods. Currently the AdamW and the AdamaxW algorithms are implemented. Different optimizers could be added.



.. _Nwalkers:

Using multiple walkers
^^^^^^^^^^^^^^^^^^^^^^

For most initial conditions the optimizer is able to find the global solution within 100'000 to 200'000 steps. 
But sometimes, the solver can find only a local minimum. Since the solver takes only a fraction of a second to run, it is most
of the time sufficient to use multiple walkers to avoid being stuck in a local minimum. As default we use 100 walkers.



.. _Nbest:

Using multiple iterations
^^^^^^^^^^^^^^^^^^^^^^^^^

To further improve the results, multiple iterations can be used. For that the best :literal:`nBestSolutions` are used to generate a new set of walkers, by add some perturbation to the best results.



Optimizers
----------

In the following all available optimizers are described. More optimizers could be added to the code.


AdamW
^^^^^


An iteration of the AdamW optimizer includes these steps:

.. math::

	\eta = \eta \cdot 0.9999 \\
	m_t = \beta_1 m_{t - 1} + (1 - \beta_1)  g_t \\
	v_t = \beta_2 v_{t - 1} + (1 - \beta_2)  (g_t \cdot g_t) \\
	\hat{m}_t = \frac{m_t}{1 - \beta_1^t} \\
	\hat{v}_t = \frac{v_t}{1 - \beta_2^t} \\
	\theta_{t + 1} = \theta_t - \eta \frac{\hat{m}_t}{ \sqrt{\hat{v}_t} + \epsilon} - \eta \lambda_{reg} \theta_t \\

	\beta_1^t = \beta_1^t \cdot \beta1 \\
	\beta_2^t = \beta_2^t \cdot \beta2

	
with

 - :math:`\eta`: learning rate, user parameter (default = 0.4).
 - :math:`\lambda_{reg}`: Weight decay, user parameter (default = 0.0001).
 - :math:`\beta_1 = 0.9`
 - :math:`\beta_2 = 0.999`
 - :math:`\epsilon = 10^{-6}`

AdamaxW
^^^^^^^


An interation of the AdamW optimizer includes these steps:

.. math::

	\eta = \eta \cdot 0.9999 \\
	m_t = \beta_1 m_{t - 1} + (1 - \beta_1)  g_t \\
	v_t = max(\beta_2 v_{t - 1}, |g_t|)  \\

	\theta_{t + 1} = \theta_t - \eta \frac{m_t}{ v_t + \epsilon} - \eta \lambda_{reg} \theta_t \\

	\beta_1^t = \beta_1^t \cdot \beta1 \\
	\beta_2^t = \beta_2^t \cdot \beta2

	
with

 - :math:`\eta`: learning rate, user parameter (default = 0.4).
 - :math:`\lambda_{reg}`: Weight decay, user parameter (default = 0.0001).
 - :math:`\beta_1 = 0.9`
 - :math:`\beta_2 = 0.999`
 - :math:`\epsilon = 10^{-6}`


.. https://optimization.cbe.cornell.edu/index.php?title=AdamW
.. https://arxiv.org/pdf/1711.05101v3
