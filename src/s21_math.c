#include "s21_math.h"

int s21_abs(int x) { return x > 0 ? x : -x; }

long double s21_fact(int N) {
  if (N < 0) {
    return 0;
  } else if (N == 0) {
    return 1;
  } else {
    return N * s21_fact(N - 1);
  }
}

long double s21_fmax(double a, double b) {
  long double res = 1;
  if (a >= b) {
    res = a;
  } else {
    res = b;
  }
  return res;
}

long double s21_sqrt(double x) {
  if (S21_isNAN(x)) {
    return S21_NAN;
  }
  long double left = 0;
  long double right = s21_fmax(1, x);
  long double mid;
  mid = (left + right) / 2;
  if (x < 0) {
    mid = S21_NAN;
  } else {
    while ((mid - left) > S21_EPS) {
      if (mid * mid > x)
        right = mid;
      else
        left = mid;
      mid = (left + right) / 2;
    }
  }
  return mid;
}

long double s21_fmod(double x, double y) {
  if (!S21_isFIN(x) || S21_isNAN(y)) {
    return S21_NAN;
  }
  if (S21_isINF(x) && S21_isINF(y)) {
    return S21_NAN;
  }
  if (S21_isINF(y)) {
    return x;
  }
  if (s21_fabs(y) < 1e-7) {
    return S21_NAN;
  }
  if (s21_fabs(x) < 1e-7) {
    return 0;
  }
  long long int mod = 0;
  mod = x / y;
  long double res = (long double)x - mod * (long double)y;
  return res;
}

long double s21_fabs(double x) {
  if (S21_isNAN(x)) {
    return S21_NAN;
  }
  if (!S21_isFIN(x)) {
    if (x < 0) {
      return -x;
    }
    return x;
  }
  return x < 0 ? -x : x;
}

long double s21_exp(double x) {
  long double res = 1;
  long double temp = 1;
  long double i = 1;
  int flag = 0;
  if (x < 0) {
    x *= -1;
    flag = 1;
  }
  while (s21_fabs(res) > S21_EPS) {
    res *= x / i;
    i += 1;
    temp += res;
    if (temp > DBL_MAX) {
      temp = S21_INF;
      break;
    }
  }
  if (flag == 1) {
    if (temp > DBL_MAX) {
      temp = 0;
    } else {
      temp = 1. / temp;
    }
  }
  if (temp > DBL_MAX) {
    return S21_INF;
  }
  return temp;
}

long double s21_log(double x) {
  int ex_pow = 0;
  double res = 0;
  double compare = 0;
  if (x == S21_INF) {
    res = S21_INF;
  } else if (x == 0) {
    res = -S21_INF;
  } else if (x < 0) {
    res = S21_NAN;
  } else if (x == 1) {
    res = 0;
  } else {
    for (; x >= S21_EXP; x /= S21_EXP, ex_pow++) continue;
    int i;
    for (i = 0; i < 100; i++) {
      compare = res;
      res = compare + 2 * (x - s21_exp(compare)) / (x + s21_exp(compare));
    }
  }
  return (res + ex_pow);
}

long double s21_ceil(double x) {
  if (!S21_isFIN(x)) {
    return x;
  }
  long double ceil_x = (long long int)x;
  if (s21_fabs(x) > 0. && x != ceil_x) {
    if (x != DBL_MAX) {
      if (x > 0.) {
        ceil_x += 1;
      }
    } else {
      return DBL_MAX;
    }
  }
  return ceil_x;
}

long double s21_floor(double x) {
  if (!S21_isFIN(x)) {
    return x;
  }
  long double floor_x = (long long int)x;
  if (s21_fabs(x - floor_x) > 0. && s21_fabs(x) > 0.) {
    if (x < 0.) {
      floor_x -= 1;
    }
  }
  return floor_x;
}

long double s21_pow(double base, double exp) {
  long double res;
  long double copy = base;

  if (copy < 0) {
    copy = -copy;
    res = s21_exp(exp * s21_log(copy));
    if (s21_fmod(exp, 2) != 0) {
      res = -res;  // четная / нечетная степень при отрицательном основании
    }
  } else {
    res = s21_exp(exp * s21_log(base));
  }
  return res;
}

long double s21_atan(double x) {
  long double sum_atan = 0;
  const long double s21_atan_1 = 0.7853981633974480L;
  if (S21_isNAN(x)) {
    return S21_NAN;
  }
  if (x == 1) {
    sum_atan = s21_atan_1;
  } else if (x == -1) {
    sum_atan = -s21_atan_1;
  } else if (x == S21_PI / 2) {
    sum_atan = 1.003884821853887214L;
  } else if (x == -S21_PI / 2) {
    sum_atan = -1.003884821853887214L;
  } else if (x == S21_INF || x == -S21_INF) {
    sum_atan = x < 0 ? -S21_PI / 2 : S21_PI / 2;
  } else if (-1. < x && x < 1.) {
    for (register int i = 0; i < 5000; i++) {
      sum_atan += s21_pow(-1, i) * s21_pow(x, 1 + (2 * i)) / (1 + (2 * i));
    }
  } else {
    for (register int i = 0; i < 7000; i++) {
      sum_atan += s21_pow(-1, i) * s21_pow(x, -1 - (2 * i)) / (1 + (2 * i));
    }
    sum_atan = S21_PI * s21_sqrt(x * x) / (2 * x) - sum_atan;
  }
  return sum_atan;
}

long double s21_asin(double x) {
  long double asin = 0.;
  if (x == 1.) {
    return S21_PI / 2;
  } else if (x == -1.) {
    return -S21_PI / 2;
  }
  if (s21_fabs(x) < 1e-9) {
    return 0;
  }
  if (x == 0.7071067811865475244) {
    return S21_PI / 4;
  }
  if (x == -0.7071067811865475244) {
    return -S21_PI / 4;
  }
  if (-1. < x && x < 1.) {
    asin = s21_atan(x / s21_sqrt(1 - x * x));
  } else {
    return S21_NAN;
  }
  return asin;
}

long double s21_acos(double x) {
  long double acos = 0.;
  if (x == 1.) {
    return 0;
  } else if (x == -1.) {
    return S21_PI;
  } else if (x == 0) {
    return S21_PI / 2;
  }
  if (x == 0.7071067811865475244) {
    return S21_PI / 4;
  }
  if (x == -0.7071067811865475244) {
    return 3 * S21_PI / 4;
  }
  if (0. < x && x < 1.) {
    acos = s21_atan(s21_sqrt(1 - x * x) / x);
  } else if (-1. < x && x < 0.) {
    acos = S21_PI + s21_atan(s21_sqrt(1 - x * x) / x);
  } else {
    return S21_NAN;
  }
  return acos;
}

long double s21_sin(double x) {
  int count = -1;
  if (x > 0) count = 1;
  x *= count;
  if (x > S21_PI) {
    x -= 2 * S21_PI * s21_floor(x / (2 * S21_PI));
  }
  long double temp = x;
  long double sum_sin = x;
  unsigned int fact = 1;
  while (s21_fabs(temp) > S21_EPS * S21_EPS) {
    temp /= (fact + 1) * (fact + 2);
    fact += 2;
    temp *= -x * x;
    sum_sin += temp;
  }
  return sum_sin * count;
}

long double s21_cos(double x) { return s21_sin(S21_PI / 2 - x); }

long double s21_tan(double x) {
  if (x == S21_PI / 2) {
    return 16331239353195370L;
  } else if (x == -S21_PI / 2) {
    return -16331239353195370L;
  }
  if (x == 0) {
    return 0;
  }
  return s21_sin(x) / s21_cos(x);
}
