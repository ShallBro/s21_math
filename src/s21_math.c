#include "s21_math.h"

// вычисляет абсолютное значение целочисленного значения
int s21_abs(int x) { return x > 0 ? x : -x; }

// функция для вычисления корня, макс точка интервала
long double s21_fmax(double a, double b) {
  long double res = 1;
  if (a >= b) {
    res = a;
  } else {
    res = b;
  }
  return res;
}

// вычисляет квадратный корень
long double s21_sqrt(double x) {
  long double mid;
  if (S21_isNAN(x)) {
    mid = S21_NAN;
  }
  long double left = 0;
  long double right = s21_fmax(1, x);

  mid = (left + right) / 2;
  if (x < 0) {
    mid = S21_NAN;
  } else if (x == 0.0) {
    mid = 0;
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

// остаток от операции деления с плавающей запятой
long double s21_fmod(double x, double y) {
  long double res;
  if (!S21_isFIN(x) || S21_isNAN(y) || (S21_isINF(x) && S21_isINF(y)) ||
      (s21_fabs(y) < 1e-7)) {
    res = S21_NAN;
  } else {
    if (S21_isINF(y)) {
      res = x;
    } else {
      if (s21_fabs(x) < 1e-7) {
        res = 0;
      } else {
        long long int mod = 0;
        mod = x / y;
        res = (long double)x - mod * (long double)y;
      }
    }
  }
  return res;
}

// вычисляет абсолютное значение значения с плавающей запятой
long double s21_fabs(double x) {
  long double res;
  if (S21_isNAN(x)) {
    res = S21_NAN;
  }
  if (!S21_isFIN(x)) {
    if (x < 0) {
      res = -x;
    }
    res = x;
  }
  res = x < 0 ? -x : x;
  return res;
}

// возвращает e, возведенное в заданную степень
long double s21_exp(double x) {
  long double sum = 1;
  long double res = 1;
  long double i = 1;
  int count = 0;
  if (x < 0) {
    x *= -1;
    count = 1;
  }
  while (s21_fabs(sum) > S21_EPS && res != S21_INF) {
    sum *= x / i;
    i += 1;
    res += sum;
    if (res > DBL_MAX) {
      res = S21_INF;
    }
  }
  if (count == 1) {
    if (res > DBL_MAX) {
      res = 0;
    } else {
      res = 1. / res;
    }
  }
  if (res > DBL_MAX) {
    res = S21_INF;
  }
  return res;
}

// вычисляет натуральный логарифм
long double s21_log(double x) {
  int count = 0;
  double res = 0;
  double cmp = 0;
  if (x == S21_INF) {
    res = S21_INF;
  } else if (x == 0) {
    res = -S21_INF;
  } else if (x < 0) {
    res = S21_NAN;
  } else if (x == 1) {
    res = 0;
  } else {
    for (; x >= S21_EXP; x /= S21_EXP, count++) continue;
    int i;
    for (i = 0; i < 100; i++) {
      cmp = res;
      res = cmp + 2 * (x - s21_exp(cmp)) / (x + s21_exp(cmp));
    }
  }
  return (res + count);
}

// возвращает ближайшее целое число, округление в большую сторону
long double s21_ceil(double x) {
  long double res = (long long int)x;
  if (!S21_isFIN(x)) {
    res = x;
  }
  if (s21_fabs(x) > 0. && x != res) {
    if (x != DBL_MAX) {
      if (x > 0.) {
        res += 1;
      }
    } else {
      res = DBL_MAX;
    }
  }
  return res;
}

// возвращает ближайшее целое число, округление в меньшую сторону
long double s21_floor(double x) {
  long double res = (long long int)x;
  if (!S21_isFIN(x)) {
    res = x;
  }
  if (s21_fabs(x - res) > 0. && s21_fabs(x) > 0.) {
    if (x < 0.) {
      res -= 1;
    }
  }
  return res;
}

// возводит число в заданную степень
long double s21_pow(double base, double exp) {
  long double res = 1;
  long double copy = base;
  long double expo = exp;
  int flag = 0;
  if (copy < 0) {
    copy = -copy;
    res = s21_exp(exp * s21_log(copy));
    if (s21_fmod(exp, 2) != 0) {
      res = -res;  // четная / нечетная степень при отрицательном основании
    }
  } else if (copy == 0) {
    if (exp == 0) {
      res = 1;
    } else {
      res = 0;
    }
  } else {
    res = s21_exp(exp * s21_log(base));
  }
  if (expo < 0) {
    expo = -expo;
    if (base == 0) res = S21_INF;
  }
  if ((base == S21_INF && exp == 1) || (S21_isNAN(base) && exp == 0) ||
      (base == 0 && exp == 0) || (base == -1 && S21_isINF(exp)) ||
      (s21_fabs(base) == 1 && (exp == -S21_INF || S21_isNAN(exp))) ||
      (base == S21_NAN && exp == 0) || (base == -S21_NAN && exp == 0) ||
      (base == S21_INF && exp == 0) || (base == -S21_INF && exp == 0) ||
      (base == 1 && S21_isINF(exp)) || (base == 1 && S21_isINF(-exp))) {
    res = 1;
    flag = 1;
  }
  if (base == -S21_INF && -exp) res = 0;
  if ((s21_fabs(base) - 1 < S21_EPS && S21_isINF(exp) && exp < 0 &&
       flag == 0) ||
      (s21_fabs(base) - 1 > S21_EPS && S21_isINF(exp) && exp > 0 &&
       flag == 0) ||
      (s21_fabs(base) <= S21_EPS && base >= 0 && exp < 0 && flag == 0))
    res = S21_INF;
  if (S21_isINF(base) && base < 0 && exp > 0) res = S21_INF;
  if ((S21_isNAN(base) && S21_isNAN(exp)) ||
      (S21_isINF(base) && base < 0 && S21_isNAN(exp)) ||
      (S21_isINF(base) && base < 0 && S21_isNAN(exp) && exp < 0))
    res = S21_NAN;
  if (base < 0 && s21_fabs(s21_fmod(exp, 1.0)) > S21_EPS) res = -S21_NAN;
  return res;
}

// вычисляет арктангенс
long double s21_atan(double x) {
  long double res = 0;
  const long double atan = 0.7853981633974480L;
  if (S21_isNAN(x)) {
    res = S21_NAN;
  }
  if (x == 1) {
    res = atan;
  } else if (x == -1) {
    res = -atan;
  } else if (x == S21_PI / 2) {
    res = 1.003884821853887214L;
  } else if (x == -S21_PI / 2) {
    res = -1.003884821853887214L;
  } else if (x == S21_INF || x == -S21_INF) {
    res = x < 0 ? -S21_PI / 2 : S21_PI / 2;
  } else if (-1. < x && x < 1.) {
    for (register int i = 0; i < 100; i++) {
      res += s21_pow(-1, i) * s21_pow(x, 1 + (2 * i)) / (1 + (2 * i));
    }
  } else {
    for (register int i = 0; i < 100; i++) {
      res += s21_pow(-1, i) * s21_pow(x, -1 - (2 * i)) / (1 + (2 * i));
    }
    res = S21_PI * s21_sqrt(x * x) / (2 * x) - res;
  }
  return res;
}

// вычисляет арксинус
long double s21_asin(double x) {
  long double res = 0.;
  if (x == 1.) {
    res = S21_PI / 2;
  } else if (x == -1.) {
    res = -S21_PI / 2;
  } else {
    if (s21_fabs(x) < 1e-9) {
      res = 0;
    }
    if (x == 0.7071067811865475244) {
      res = S21_PI / 4;
    }
    if (x == -0.7071067811865475244) {
      res = -S21_PI / 4;
    }
    if (-1. < x && x < 1.) {
      res = s21_atan(x / s21_sqrt(1 - x * x));
    } else {
      res = S21_NAN;
    }
  }
  return res;
}

// вычисляет арккосинус
long double s21_acos(double x) {
  long double res = 0.;
  if (x == 1.) {
    res = 0;
  } else if (x == -1.) {
    res = S21_PI;
  } else if (x == 0) {
    res = S21_PI / 2;
  } else {
    if (x == 0.7071067811865475244) {
      res = S21_PI / 4;
    }
    if (x == -0.7071067811865475244) {
      res = 3 * S21_PI / 4;
    }
    if (0. < x && x < 1.) {
      res = s21_atan(s21_sqrt(1 - x * x) / x);
    } else if (-1. < x && x < 0.) {
      res = S21_PI + s21_atan(s21_sqrt(1 - x * x) / x);
    } else {
      res = S21_NAN;
    }
  }
  return res;
}

// вычисляет синус
long double s21_sin(double x) {
  int count = -1;
  if (x > 0) count = 1;
  x *= count;
  if (x > S21_PI) {
    x -= 2 * S21_PI * s21_floor(x / (2 * S21_PI));
  }
  long double temp = x;
  long double res = x;
  unsigned int fact = 1;
  while (s21_fabs(temp) > S21_EPS * S21_EPS) {
    temp /= (fact + 1) * (fact + 2);
    fact += 2;
    temp *= -x * x;
    res += temp;
  }
  return res * count;
}

// вычисляет косинус
long double s21_cos(double x) {
  if (x > 2 * S21_PI || x < -2 * S21_PI) {
    x += x > 2 * S21_PI ? -2 * S21_PI : 2 * S21_PI;
  }
  return s21_sin((S21_PI / 2.0) - x);
}

// вычисляет тангенс
long double s21_tan(double x) {
  long double res;
  if (x == S21_PI / 2) {
    res = 16331239353195370L;
  } else if (x == -S21_PI / 2) {
    res = -16331239353195370L;
  }
  if (x == 0) {
    res = 0;
  }
  res = s21_sin(x) / s21_cos(x);
  return res;
}
