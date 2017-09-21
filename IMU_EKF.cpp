#include "IMU_EKF.h"

static void mrdivide(const float A[12], const float B[9], float y[12]);
static float norm(const float x[3]);
static void mrdivide(const float A[12], const float B[9], float y[12])
{
  int rtemp;
  int r1;
  float b_A[9];
  int r2;
  int r3;
  float maxval;
  float a21;
  for (rtemp = 0; rtemp < 9; rtemp++) {
    b_A[rtemp] = B[rtemp];
  }

  r1 = 0;
  r2 = 1;
  r3 = 2;
  maxval = (float)fabs(B[0]);
  a21 = (float)fabs(B[1]);
  if (a21 > maxval) {
    maxval = a21;
    r1 = 1;
    r2 = 0;
  }

  if ((float)fabs(B[2]) > maxval) {
    r1 = 2;
    r2 = 1;
    r3 = 0;
  }

  b_A[r2] = B[r2] / B[r1];
  b_A[r3] /= b_A[r1];
  b_A[3 + r2] -= b_A[r2] * b_A[3 + r1];
  b_A[3 + r3] -= b_A[r3] * b_A[3 + r1];
  b_A[6 + r2] -= b_A[r2] * b_A[6 + r1];
  b_A[6 + r3] -= b_A[r3] * b_A[6 + r1];
  if ((float)fabs(b_A[3 + r3]) > (float)fabs(b_A[3 + r2])) {
    rtemp = r2;
    r2 = r3;
    r3 = rtemp;
  }

  b_A[3 + r3] /= b_A[3 + r2];
  b_A[6 + r3] -= b_A[3 + r3] * b_A[6 + r2];
  for (rtemp = 0; rtemp < 4; rtemp++) {
    y[rtemp + (r1 << 2)] = A[rtemp] / b_A[r1];
    y[rtemp + (r2 << 2)] = A[4 + rtemp] - y[rtemp + (r1 << 2)] * b_A[3 + r1];
    y[rtemp + (r3 << 2)] = A[8 + rtemp] - y[rtemp + (r1 << 2)] * b_A[6 + r1];
    y[rtemp + (r2 << 2)] /= b_A[3 + r2];
    y[rtemp + (r3 << 2)] -= y[rtemp + (r2 << 2)] * b_A[6 + r2];
    y[rtemp + (r3 << 2)] /= b_A[6 + r3];
    y[rtemp + (r2 << 2)] -= y[rtemp + (r3 << 2)] * b_A[3 + r3];
    y[rtemp + (r1 << 2)] -= y[rtemp + (r3 << 2)] * b_A[r3];
    y[rtemp + (r1 << 2)] -= y[rtemp + (r2 << 2)] * b_A[r2];
  }
}

static float norm(const float x[3])
{
  float y;
  float scale;
  int k;
  float absxk;
  float t;
  y = 0.0F;
  scale = 1.17549435E-38F;
  for (k = 0; k < 3; k++) {
    absxk = (float)fabs(x[k]);
    if (absxk > scale) {
      t = scale / absxk;
      y = 1.0F + y * t * t;
      scale = absxk;
    } else {
      t = absxk / scale;
      y += t * t;
    }
  }

  return scale * (float)sqrt(y);
}

void IMU_EKF(float P[16], float q[4], const float Cov_info[3], const float
             omega[3], float accel[3], float mag[3], float dt, signed char ini,
             signed char use_mag)
{
  int i;
  float y;
  static const signed char iv0[16] = { 10, 0, 0, 0, 0, 10, 0, 0, 0, 0, 10, 0, 0,
    0, 0, 10 };

  signed char I[16];
  float b_y[16];
  int i0;
  float Ak[4];
  float b_Ak[16];
  int i1;
  float Hk[12];
  float S[9];
  float b_Hk[12];
  float Kk[12];
  float fv4[3];
  float b_accel[3];
  float qe[4];
  float scale;
  float absxk;
  float t;
  if (ini == 1) {
    for (i = 0; i < 16; i++) {
      P[i] = iv0[i];
    }

    for (i = 0; i < 4; i++) {
      q[i] = 0.0F;
    }

    q[0] = 1.0F;
  }

  y = dt / 2.0F;
  for (i = 0; i < 16; i++) {
    I[i] = 0;
  }

  for (i = 0; i < 4; i++) {
    I[i + (i << 2)] = 1;
  }

  b_y[0] = 0.0F;
  b_y[4] = y * -omega[0];
  b_y[8] = y * -omega[1];
  b_y[12] = y * -omega[2];
  b_y[1] = y * omega[0];
  b_y[5] = 0.0F;
  b_y[9] = y * omega[2];
  b_y[13] = y * -omega[1];
  b_y[2] = y * omega[1];
  b_y[6] = y * -omega[2];
  b_y[10] = 0.0F;
  b_y[14] = y * omega[0];
  b_y[3] = y * omega[2];
  b_y[7] = y * omega[1];
  b_y[11] = y * -omega[0];
  b_y[15] = 0.0F;
  for (i = 0; i < 4; i++) {
    for (i0 = 0; i0 < 4; i0++) {
      b_Ak[i0 + (i << 2)] = b_y[i0 + (i << 2)] + (float)I[i0 + (i << 2)];
    }
  }

  for (i = 0; i < 4; i++) {
    Ak[i] = 0.0F;
    for (i0 = 0; i0 < 4; i0++) {
      Ak[i] += b_Ak[i + (i0 << 2)] * q[i0];
    }
  }

  for (i = 0; i < 4; i++) {
    q[i] = Ak[i];
    for (i0 = 0; i0 < 4; i0++) {
      b_y[i + (i0 << 2)] = 0.0F;
      for (i1 = 0; i1 < 4; i1++) {
        b_y[i + (i0 << 2)] += b_Ak[i + (i1 << 2)] * P[i1 + (i0 << 2)];
      }
    }
  }

  for (i = 0; i < 4; i++) {
    for (i0 = 0; i0 < 4; i0++) {
      P[i + (i0 << 2)] = 0.0F;
      for (i1 = 0; i1 < 4; i1++) {
        P[i + (i0 << 2)] += b_y[i + (i1 << 2)] * b_Ak[i0 + (i1 << 2)];
      }
    }
  }

  for (i = 0; i < 4; i++) {
    P[i + (i << 2)] += Cov_info[0];
  }

  y = norm(accel);
  Hk[0] = -2.0F * q[2];
  Hk[3] = 2.0F * q[3];
  Hk[6] = -2.0F * q[0];
  Hk[9] = 2.0F * q[1];
  Hk[1] = 2.0F * q[1];
  Hk[4] = 2.0F * q[0];
  Hk[7] = 2.0F * q[3];
  Hk[10] = 2.0F * q[2];
  Hk[2] = 2.0F * q[0];
  Hk[5] = -2.0F * q[1];
  Hk[8] = -2.0F * q[2];
  Hk[11] = 2.0F * q[3];
  for (i = 0; i < 3; i++) {
    for (i0 = 0; i0 < 4; i0++) {
      b_Hk[i + 3 * i0] = 0.0F;
      for (i1 = 0; i1 < 4; i1++) {
        b_Hk[i + 3 * i0] += Hk[i + 3 * i1] * P[i1 + (i0 << 2)];
      }
    }

    for (i0 = 0; i0 < 3; i0++) {
      S[i + 3 * i0] = 0.0F;
      for (i1 = 0; i1 < 4; i1++) {
        S[i + 3 * i0] += b_Hk[i + 3 * i1] * Hk[i0 + 3 * i1];
      }
    }

    accel[i] /= y;
  }

  for (i = 0; i < 3; i++) {
    S[i + 3 * i] += Cov_info[1];
  }

  for (i = 0; i < 4; i++) {
    for (i0 = 0; i0 < 3; i0++) {
      b_Hk[i + (i0 << 2)] = 0.0F;
      for (i1 = 0; i1 < 4; i1++) {
        b_Hk[i + (i0 << 2)] += P[i + (i1 << 2)] * Hk[i0 + 3 * i1];
      }
    }
  }

  mrdivide(b_Hk, S, Kk);
  fv4[0] = 2.0F * q[1] * q[3] - 2.0F * q[0] * q[2];
  fv4[1] = 2.0F * q[0] * q[1] + 2.0F * q[2] * q[3];
  fv4[2] = ((q[0] * q[0] - q[1] * q[1]) - q[2] * q[2]) + q[3] * q[3];
  for (i = 0; i < 3; i++) {
    b_accel[i] = accel[i] - fv4[i];
  }

  for (i = 0; i < 4; i++) {
    qe[i] = 0.0F;
    for (i0 = 0; i0 < 3; i0++) {
      qe[i] += Kk[i + (i0 << 2)] * b_accel[i0];
    }
  }

  qe[3] = 0.0F;
  for (i = 0; i < 4; i++) {
    q[i] += qe[i];
    for (i0 = 0; i0 < 4; i0++) {
      b_y[i + (i0 << 2)] = 0.0F;
      for (i1 = 0; i1 < 3; i1++) {
        b_y[i + (i0 << 2)] += Kk[i + (i1 << 2)] * Hk[i1 + 3 * i0];
      }
    }

    for (i0 = 0; i0 < 4; i0++) {
      y = 0.0F;
      for (i1 = 0; i1 < 4; i1++) {
        y += b_y[i + (i1 << 2)] * P[i1 + (i0 << 2)];
      }

      b_Ak[i + (i0 << 2)] = P[i + (i0 << 2)] - y;
    }
  }

  for (i = 0; i < 4; i++) {
    for (i0 = 0; i0 < 4; i0++) {
      P[i0 + (i << 2)] = b_Ak[i0 + (i << 2)];
    }
  }

  if (use_mag == 1) {
    y = norm(mag);
    Hk[0] = 2.0F * q[3];
    Hk[3] = 2.0F * q[2];
    Hk[6] = 2.0F * q[1];
    Hk[9] = 2.0F * q[0];
    Hk[1] = 2.0F * q[0];
    Hk[4] = -2.0F * q[1];
    Hk[7] = -2.0F * q[2];
    Hk[10] = -2.0F * q[3];
    Hk[2] = -2.0F * q[1];
    Hk[5] = -2.0F * q[0];
    Hk[8] = 2.0F * q[3];
    Hk[11] = 2.0F * q[2];
    for (i = 0; i < 3; i++) {
      for (i0 = 0; i0 < 4; i0++) {
        b_Hk[i + 3 * i0] = 0.0F;
        for (i1 = 0; i1 < 4; i1++) {
          b_Hk[i + 3 * i0] += Hk[i + 3 * i1] * P[i1 + (i0 << 2)];
        }
      }

      for (i0 = 0; i0 < 3; i0++) {
        S[i + 3 * i0] = 0.0F;
        for (i1 = 0; i1 < 4; i1++) {
          S[i + 3 * i0] += b_Hk[i + 3 * i1] * Hk[i0 + 3 * i1];
        }
      }

      mag[i] /= y;
    }

    for (i = 0; i < 3; i++) {
      S[i + 3 * i] += Cov_info[2];
    }

    for (i = 0; i < 4; i++) {
      for (i0 = 0; i0 < 3; i0++) {
        b_Hk[i + (i0 << 2)] = 0.0F;
        for (i1 = 0; i1 < 4; i1++) {
          b_Hk[i + (i0 << 2)] += P[i + (i1 << 2)] * Hk[i0 + 3 * i1];
        }
      }
    }

    mrdivide(b_Hk, S, Kk);
    fv4[0] = 2.0F * q[1] * q[2] + 2.0F * q[0] * q[3];
    fv4[1] = ((q[0] * q[0] - q[1] * q[1]) - q[2] * q[2]) - q[3] * q[3];
    fv4[2] = 2.0F * q[2] * q[3] - 2.0F * q[0] * q[1];
    for (i = 0; i < 3; i++) {
      b_accel[i] = mag[i] - fv4[i];
    }

    for (i = 0; i < 4; i++) {
      qe[i] = 0.0F;
      for (i0 = 0; i0 < 3; i0++) {
        qe[i] += Kk[i + (i0 << 2)] * b_accel[i0];
      }
    }

    qe[1] = 0.0F;
    qe[2] = 0.0F;
    for (i = 0; i < 4; i++) {
      q[i] += qe[i];
      for (i0 = 0; i0 < 4; i0++) {
        b_y[i + (i0 << 2)] = 0.0F;
        for (i1 = 0; i1 < 3; i1++) {
          b_y[i + (i0 << 2)] += Kk[i + (i1 << 2)] * Hk[i1 + 3 * i0];
        }
      }

      for (i0 = 0; i0 < 4; i0++) {
        y = 0.0F;
        for (i1 = 0; i1 < 4; i1++) {
          y += b_y[i + (i1 << 2)] * P[i1 + (i0 << 2)];
        }

        b_Ak[i + (i0 << 2)] = P[i + (i0 << 2)] - y;
      }
    }

    for (i = 0; i < 4; i++) {
      for (i0 = 0; i0 < 4; i0++) {
        P[i0 + (i << 2)] = b_Ak[i0 + (i << 2)];
      }
    }
  }

  y = 0.0F;
  scale = 1.17549435E-38F;
  for (i = 0; i < 4; i++) {
    absxk = (float)fabs(q[i]);
    if (absxk > scale) {
      t = scale / absxk;
      y = 1.0F + y * t * t;
      scale = absxk;
    } else {
      t = absxk / scale;
      y += t * t;
    }
  }

  y = scale * (float)sqrt(y);
  for (i = 0; i < 4; i++) {
    q[i] /= y;
  }
}
