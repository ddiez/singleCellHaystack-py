import singleCellHaystack as hs
import numpy as np

m = np.array([[0, 1], [2, 3]])

def test_scale():
  m_scaled, m_mean, m_std = hs.scale(m)
  assert np.all(m_scaled == np.array([[-1, -1], [1, 1]]))

  m_scaled, m_mean, m_std = hs.scale(m, mean=10, std=2)
  assert np.all(m_scaled == np.array([[-5, -4.5], [-4, -3.5]]))

def test_scale_inv():
  m_scaled, m_mean, m_std = hs.scale(m)
  m_back = hs.scale_inv(m_scaled, m_mean, m_std)
  assert np.all(m == m_back)