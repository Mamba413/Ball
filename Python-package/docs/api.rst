API List
========

Ball functions
-----------------
A complete list of all ball functions provided by Ball.

Ball divergence, ball covariance, and ball correlation statistics
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
These functions compute the estimators for the ball divergence, ball covariance, and ball correlation.

.. autosummary::
   :toctree: functions
   
   Ball.bd
   Ball.bcov
   Ball.bcor
   
Ball Divergence based Equality of Distributions Test
^^^^^^^^^^^^^^^^^^^^^
Performs the nonparametric two-sample or K-sample Ball Divergence test for equality of multivariate distributions.

.. autosummary::
   :toctree: functions
   
   Ball.bd_test

Ball Covariance based (joint) independence Test
^^^^^^^^^^^^^^^^^^^^^^^^^
Ball Covariance test of independence. Ball covariance are generic dependence measures in Banach spaces.

.. autosummary::
   :toctree: functions
   
   Ball.bcov_test

Ball Correlation based Sure Independence Screening (BCor-SIS)
^^^^^^^^^^^^^^^^^^^^^^^^^
Generic non-parametric sure independence screening (SIS) procedure based on Ball Correlation. Ball correlation is a generic measure of dependence in Banach spaces.

.. autosummary::
   :toctree: functions
   
   Ball.bcorsis