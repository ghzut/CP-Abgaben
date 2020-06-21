#include<Eigen/Dense>
using namespace std;

Vector2d verlet(function<Vector2d(Vector2d)> F, Vector4d y, Vector4d y_m_1, double &c)
{
  double m = 1.;
  double h = 0.01;
  Vector4d ret;
  Vector2d r_n_m_1 << y_m_1(0),y_m_1(1);
  Vector2d r_n << y(0), y(1);
  Vector2d v_n_m_1 << y_m_1(2),y_m_1(3);
  Vector2d v_n << y(2), y(3);
  Vector2d a_n = F(r_n)/m;
  Vector2d r_n_p_1 = 2 * r_n - r_n_m_1 + a_n * pow(h,2.);
  Vector2d v_n_p_1 = r_n_p_1 - r_n_m_1;
  v_n_p_1 /= 2*h;
  ret << r_n_p_1(0), r_n_p_1(1), v_n_p_1(0), v_n_p_1(1); 
  return ret;
}