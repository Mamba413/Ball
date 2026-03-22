# -*- coding: utf-8 -*-
"""
Empirical Metric Distribution Function (EMDF) and Metric Cramér-von Mises statistic.

References
----------
.. [1] Zhu et al. (2021), Journal of Statistical Software, 97(6), 1-31.
"""

import math
import numpy as np
from scipy.spatial.distance import cdist, pdist, squareform


class MDF:
    """Empirical Metric Distribution Function (EMDF).

    Computes the empirical metric distribution function using rank-based
    acceleration. Supports subsampling of ball centers for scalability.

    Parameters
    ----------
    metric : str, default='euclidean'
        Distance metric passed to ``scipy.spatial.distance.cdist``.
    subsample : str, int, or None, default=None
        Controls how many ball centers to use.

        - None: use all n training points as centers.
        - ``'sqrt'``: use ``floor(sqrt(n))`` randomly selected centers.
        - ``'log2'``: use ``floor(log2(n))`` randomly selected centers.
        - int: use exactly this many centers (clamped to n).

    Examples
    --------
    >>> import numpy as np
    >>> np.random.seed(42)
    >>> X = np.random.randn(100, 2)
    >>> mdf = MDF()
    >>> mdf.fit(X)
    MDF(metric='euclidean', subsample=None)
    >>> emdf = mdf.predict()
    >>> emdf.shape
    (100, 100)
    """

    def __init__(self, metric='euclidean', subsample=None):
        self.metric = metric
        self.subsample = subsample

    def __repr__(self):
        return f"MDF(metric={self.metric!r}, subsample={self.subsample!r})"

    def fit(self, X):
        """Fit the EMDF with training samples.

        Stores all n training points for counting and selects m centers
        for evaluation.

        Precomputation (rank-based acceleration from BD.c):
          1. Compute distances from m centers to all n training points.
          2. Sort distances per center.
          3. Compute center-to-center distances.
          4. Compute rank matrix via searchsorted.

        Parameters
        ----------
        X : array-like, shape (n_samples, n_features)
            Training data.

        Returns
        -------
        self
        """
        X = np.asarray(X, dtype=float)
        if X.ndim == 1:
            X = X.reshape(-1, 1)

        self.X_train_ = X
        self.n_train_ = X.shape[0]
        n = self.n_train_

        # Determine centers
        m = self._resolve_subsample(n)
        if m >= n:
            self.center_indices_ = np.arange(n)
        else:
            rng = np.random.default_rng()
            self.center_indices_ = np.sort(rng.choice(n, size=m, replace=False))

        self.centers_ = X[self.center_indices_]
        self.n_centers_ = len(self.center_indices_)
        m = self.n_centers_

        # Distances from each center to all training points: (m, n)
        self.dist_centers_to_train_ = cdist(self.centers_, X, metric=self.metric)

        # Sorted distances per center: (m, n)
        self.sorted_centers_to_train_ = np.sort(self.dist_centers_to_train_, axis=1)

        # Pairwise distances among centers: (m, m)
        self.dist_centers_ = cdist(self.centers_, self.centers_, metric=self.metric)

        # Rank matrix: rank_matrix_[a, b] = #{l : d(X_l, center_a) <= d(center_a, center_b)}
        self.rank_matrix_ = np.empty((m, m), dtype=int)
        for a in range(m):
            self.rank_matrix_[a] = np.searchsorted(
                self.sorted_centers_to_train_[a],
                self.dist_centers_[a],
                side='right'
            )

        return self

    def predict(self, X=None):
        """Compute the EMDF matrix.

        Parameters
        ----------
        X : array-like, shape (p, n_features), or None
            Query points. If None, computes EMDF on the fitted centers
            using precomputed ranks.

        Returns
        -------
        emdf : ndarray
            - ``X=None``: shape ``(m, m)``, where ``emdf[a, b] = F(center_a, center_b)``.
            - ``X`` given: shape ``(p, p)``, where ``emdf[i, j] = F(X_i, X_j)``.

        EMDF(u, v) = (1 / n_train) * #{l : d(X_train_l, u) <= d(u, v)}
        """
        n = self.n_train_

        if X is None:
            # Use precomputed rank matrix
            return self.rank_matrix_ / n
        else:
            X = np.asarray(X, dtype=float)
            if X.ndim == 1:
                X = X.reshape(-1, 1)

            p = X.shape[0]

            # Distances from each query point to all training points: (p, n)
            dist_to_train = cdist(X, self.X_train_, metric=self.metric)

            # Sort per query point: (p, n)
            sorted_dist = np.sort(dist_to_train, axis=1)

            # Pairwise distances among query points: (p, p)
            dist_query = cdist(X, X, metric=self.metric)

            # Compute EMDF via searchsorted
            emdf = np.empty((p, p), dtype=float)
            for i in range(p):
                ranks = np.searchsorted(sorted_dist[i], dist_query[i], side='right')
                emdf[i] = ranks / n

            return emdf

    def _resolve_subsample(self, n):
        """Resolve the subsample parameter to an integer count of centers."""
        s = self.subsample
        if s is None:
            return n
        if isinstance(s, str):
            s = s.lower()
            if s == 'sqrt':
                return max(1, math.floor(math.sqrt(n)))
            elif s == 'log2':
                return max(1, math.floor(math.log2(n)))
            else:
                raise ValueError(
                    f"subsample must be None, 'sqrt', 'log2', or int, got {s!r}"
                )
        if isinstance(s, (int, np.integer)):
            return max(1, min(int(s), n))
        raise TypeError(
            f"subsample must be None, str, or int, got {type(s).__name__}"
        )


def mcvm(mdf_list, groups, sigma2=None):
    """Compute the Metric Cramer-von Mises statistic.

    Measures distributional divergence across K groups using the EMDF.

    Parameters
    ----------
    mdf_list : list of MDF
        K fitted MDF objects, one per group.
    groups : list of array-like
        K arrays of sample data corresponding to ``mdf_list``.
    sigma2 : float or None, default=None
        Bandwidth for the Gaussian weight ``w(u,v) = exp(-d(u,v)^2 / (2*sigma2))``.
        If None, uses the median heuristic on the pooled sample.

    Returns
    -------
    statistic : float
        The MCVM statistic.

    Examples
    --------
    >>> import numpy as np
    >>> np.random.seed(42)
    >>> X1 = np.random.randn(50, 2)
    >>> X2 = np.random.randn(50, 2) + 2
    >>> mdf1, mdf2 = MDF().fit(X1), MDF().fit(X2)
    >>> stat = mcvm([mdf1, mdf2], [X1, X2])
    >>> stat > 0
    True
    """
    K = len(mdf_list)
    groups = [np.asarray(g, dtype=float) for g in groups]
    for i, g in enumerate(groups):
        if g.ndim == 1:
            groups[i] = g.reshape(-1, 1)

    n_total = sum(g.shape[0] for g in groups)

    # Pool all groups and fit a pooled MDF
    pooled_X = np.vstack(groups)
    subsample_setting = mdf_list[0].subsample
    metric_setting = mdf_list[0].metric
    mdf_pooled = MDF(metric=metric_setting, subsample=subsample_setting)
    mdf_pooled.fit(pooled_X)

    # Compute sigma2 via median heuristic if not provided
    if sigma2 is None:
        # Use distances among pooled MDF centers for efficiency
        dists = squareform(pdist(mdf_pooled.centers_, metric=metric_setting))
        # Extract upper triangle (excluding diagonal)
        upper = dists[np.triu_indices_from(dists, k=1)]
        sigma2 = float(np.median(upper)) ** 2
        if sigma2 == 0.0:
            sigma2 = 1.0  # fallback to avoid division by zero

    statistic = 0.0

    for k in range(K):
        n_k = groups[k].shape[0]
        p_k = n_k / n_total
        m_k = mdf_list[k].n_centers_

        # Group EMDF on its own centers
        F_k = mdf_list[k].predict()  # (m_k, m_k)

        # Pooled EMDF evaluated at group k's centers
        F_pooled = mdf_pooled.predict(mdf_list[k].centers_)  # (m_k, m_k)

        # Gaussian weight matrix
        dist_sq = mdf_list[k].dist_centers_ ** 2
        W = np.exp(-dist_sq / (2.0 * sigma2))

        # Squared difference
        diff_sq = (F_k - F_pooled) ** 2

        # Accumulate: p_k^4 * sum(W * diff_sq) / m_k^2
        # Note: p_k^2 from the formula, and another p_k^2 from the sum normalization
        # Per the plan: p_k^4 * sum / m_k^2
        statistic += (p_k ** 4) * np.sum(W * diff_sq) / (m_k ** 2)

    return statistic
