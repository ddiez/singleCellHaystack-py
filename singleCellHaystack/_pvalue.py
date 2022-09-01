import numpy as np

def estimate_spline_param(x, y, method="bs"):
  from sklearn.preprocessing import SplineTransformer
  from sklearn.preprocessing import FunctionTransformer
  from sklearn.linear_model import LinearRegression
  from sklearn.pipeline import Pipeline
  from sklearn.model_selection import GridSearchCV
  from patsy import cr

  if method == "bs":
    pip = Pipeline([
      ["transformer", SplineTransformer()],
      ["estimator", LinearRegression()]
    ])

    min_knots = 2
    max_knots = 10
    min_degree = 1
    max_degree = 5
    param = {
      "transformer__n_knots": list(range(min_knots, max_knots)),
      "transformer__degree": list(range(min_degree, max_degree))
    }

    cv = GridSearchCV(pip, param, cv=10)
    cv_res = cv.fit(x, y)

    info = {
      "method": "bs",
      "n_knots": cv_res.best_params_["transformer__n_knots"],
      "degree": cv_res.best_params_["transformer__degree"]
    }

  if method == "ns":
    NaturalSplineTransformer = FunctionTransformer(cr)
    pip = Pipeline([
      ["transformer", NaturalSplineTransformer],
      ["estimator", LinearRegression()]
    ])

    # Set dict of degrees of freedom for CV.
    min_df = 3
    max_df = 10
    df = []
    for k in range(min_df, max_df):
      df.append({"df": k})

    param = {
      "transformer__kw_args": df
    }

    cv = GridSearchCV(pip, param, cv=10)
    cv_res = cv.fit(x, y)

    info = {
      "method": "ns",
      "df": cv_res.best_params_["transformer__kw_args"]["df"]
    }

  return info

def calculate_KLD_fit(x, y, method="bs"):
  from sklearn.preprocessing import SplineTransformer
  from sklearn.preprocessing import FunctionTransformer
  from sklearn.linear_model import LinearRegression
  from sklearn.pipeline import Pipeline
  from patsy import cr

  info = estimate_spline_param(x, y, method=method)
  if method == "bs":
    n_knots = info["n_knots"]
    degree = info["degree"]

    pip = Pipeline([
      ["transformer", SplineTransformer(n_knots=n_knots, degree=degree)],
      ["estimator", LinearRegression()]
    ])
    model = pip.fit(x, y)
    y_hat = model.predict(x)

  if method == "ns":
    df = info["df"]

    NaturalSplineTransformer = FunctionTransformer(cr, kw_args={"df": df})
    pip = Pipeline([
      ["transformer", NaturalSplineTransformer],
      ["estimator", LinearRegression()]
    ])
    model = pip.fit(x, y)
    y_hat = model.predict(x)

  return {
    "model": model,
    "spline": info,
    "y_hat": y_hat
  }

def calculate_Pval(KLD, KLD_rand, cv, cv_rand, method="bs", verbose=False):
  from sklearn.preprocessing import SplineTransformer
  from sklearn.linear_model import LinearRegression
  from sklearn.pipeline import Pipeline
  from patsy import cr
  from scipy.stats import norm

  if (verbose):
    print("> calculating P values ...")

  KLD_log = np.log(KLD)
  KLD_rand_log = np.log(KLD_rand)
  KLD_rand_mean = np.mean(KLD_rand_log, axis=1)
  KLD_rand_sd = np.std(KLD_rand_log, axis=1)
  cv_log = np.log(cv).reshape(-1, 1)
  cv_rand_log = np.log(cv_rand).reshape(-1, 1)

  KLD_rand_mean_fit = calculate_KLD_fit(cv_rand_log, KLD_rand_mean, method=method)
  KLD_rand_sd_fit = calculate_KLD_fit(cv_rand_log, KLD_rand_sd, method=method)

  KLD_rand_mean_model = KLD_rand_mean_fit["model"]
  KLD_rand_sd_model = KLD_rand_sd_fit["model"]

  KLD_mean = KLD_rand_mean_model.predict(cv_log)
  KLD_sd = KLD_rand_sd_model.predict(cv_log)

  logpval = norm.logsf(KLD_log, loc=KLD_mean, scale=KLD_sd)/np.log(10)
  pval = 10 ** logpval

  return {
    "pval": pval,
    "logpval": logpval,
    "method": method,
    "CV": cv_rand_log,
    "rand_mean": KLD_rand_mean,
    "rand_sd": KLD_rand_sd,
    "rand_mean_model": KLD_rand_mean_model,
    "rand_sd_model": KLD_rand_sd_model,
    "rand_mean_spline": KLD_rand_mean_fit["spline"],
    "rand_sd_spline": KLD_rand_sd_fit["spline"],
    "rand_mean_hat": KLD_rand_mean_fit["y_hat"],
    "rand_sd_hat": KLD_rand_sd_fit["y_hat"]
  }
