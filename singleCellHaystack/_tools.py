def scale(x, mean=None, std=None):
  import numpy as np

  if mean is None:
    mean = np.mean(x, axis=0)
  
  if std is None:
    std = np.std(x, axis=0)

  x = (x - mean) / std
  return x, mean, std


def scale_inv(x, mean, std):
  x = x * std + mean
  return(x)