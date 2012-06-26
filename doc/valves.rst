Valve Model
===========
.. plot::

    import matplotlib.pyplot as plt
    import numpy as np
    
    fig=plt.figure()
    ax1=fig.add_subplot(111)
    ax1.fill(np.r_[0.2,0.4,0.4,0.2,0.2],np.r_[0,0,1,1,0],'grey')
    ax1.text(0.41,0.66,r'$\leftarrow p_{low}A_{valve}$',size=20,ha='left',va='center')
    ax1.text(0.41,0.33,r'$\leftarrow k_{valve}x$',size=20,ha='left',va='center')
    ax1.text(0.19,0.66,r'$p_{high}A_{valve}\rightarrow$',size=20,ha='right',va='center')
    ax1.text(0.19,0.33,r'$\frac{1}{2}C_D\rho (V-\dot x)^2 A_{valve}\rightarrow$',size=20,ha='right',va='center')
    ax1.set_xlim(0.1-1,0.1+1)
    ax1.axis('equal')
    ax1.axis('off')
    ax1.set_title('Pressure-dominant Free-body-diagram')
    plt.show()

Pressure-dominant Regime

.. math::

    M_{valve}\ddot x_{valve}+k_{valve}x_{valve} = (p_{high}-p_{low}) A_{valve}+\frac{1}{2}C_D\rho (V-\dot x_{valve})^2A_{valve}
    
Two variables are :math:`x_2=\dot x_{valve}` and :math:`x_1=x_{valve}` where :math:`\ddot x_{valve}=\frac{d}{dt}[\dot x_{valve}]` or :math:`\ddot x_{valve}=\dot x_2`.  Thus the system of derivatives is

.. math::

    \mathbf f_{valves}=\frac{d}{dt}\left[ \begin{array}{c} \dot x_{valve} \\ x_{valve} \end{array} \right]=\left[ \begin{array}{c} \frac{d}{dt}[\dot x_{valve}] \\ \frac{d}{dt}[x_{valve}] \end{array} \right]

Thus the system of equations is given by

.. math::

    \dot x_1 = x_2
    
    M_{valve}\dot x_2+k_{valve}x_1 = (p_{high}-p_{low}) A_{valve}+\frac{1}{2}C_D\rho (V-x_2)^2A_{valve}

which yields the solution for the derivatives of :math:`x_1` and :math:`x_2` of

.. math::

    \dot x_1 = x_2
    
    \dot x_2 = \dfrac{(p_{high}-p_{low}) A_{valve}+\frac{1}{2}C_D\rho (V-x_2)^2A_{valve}-k_{valve}x_1}{M_{valve}}
    
.. plot::

    import matplotlib.pyplot as plt
    import numpy as np
    
    fig=plt.figure()
    ax1=fig.add_subplot(111)
    ax1.fill(np.r_[0.2,0.4,0.4,0.2,0.2],np.r_[0,0,1,1,0],'grey')
    ax1.text(0.41,0.66,r'$\leftarrow p_{low}A_{valve}$',size=20,ha='left',va='center')
    ax1.text(0.41,0.33,r'$\leftarrow k_{valve}x$',size=20,ha='left',va='center')
    ax1.text(0.19,0.8,r'$p_{low}A_{valve} \rightarrow$',size=20,ha='right',va='center')
    ax1.text(0.19,0.5,r'$\frac{1}{2}C_D\rho (V-\dot x)^2 A_{valve} \rightarrow$',size=20,ha='right',va='center')
    ax1.text(0.19,0.2,r'$\rho (V-\dot x)^2 A_{port} \rightarrow$',size=20,ha='right',va='center')
    ax1.set_xlim(0.1-1,0.1+1)
    ax1.axis('equal')
    ax1.axis('off')
    ax1.set_title('Mass-flux-dominant Free-body-diagram')
    plt.show()

And if mass-flux-dominated, force balance given by

.. math::

    M_{valve}\ddot x_{valve}+k_{valve}x_{valve} = \frac{1}{2}C_D\rho (\mathbf V-\dot x_{valve})^2 A_{valve}+\rho (\mathbf V-\dot x_{valve})^2 A_{port}

Which yields the solution for the system of derivatives of 

.. math::

    \dot x_1 = x_2
    
    \dot x_2= \dfrac{\frac{1}{2}C_D\rho (\mathbf V-x_2)^2 A_{valve}+\rho (\mathbf V-x_2)^2 A_{port}-k_{valve}x_1}{M_{valve}}
