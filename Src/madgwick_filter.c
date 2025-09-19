// madgwick_filter.c
#include "madgwick_filter.h"

static double beta_g = 0.1;

void madgwick_init(Quaternion* q) {
    q->w = 1.0f;
    q->x = 0.0f;
    q->y = 0.0f;
    q->z = 0.0f;
}

static double inv_sqrt(double x) {
    return 1.0f / sqrtf(x);
}

void madgwick_update(Quaternion* q, double ax, double ay, double az, 
                     double gx, double gy, double gz, double dt) {
    double recipNorm;
    double s0, s1, s2, s3;
    double qDot1, qDot2, qDot3, qDot4;
    double _2q0, _2q1, _2q2, _2q3, _4q0, _4q1, _4q2, _8q1, _8q2, q0q0, q1q1, q2q2, q3q3;

    qDot1 = 0.5f * (-q->x * gx - q->y * gy - q->z * gz);
    qDot2 = 0.5f * (q->w * gx + q->y * gz - q->z * gy);
    qDot3 = 0.5f * (q->w * gy - q->x * gz + q->z * gx);
    qDot4 = 0.5f * (q->w * gz + q->x * gy - q->y * gx);

    if(!((ax == 0.0f) && (ay == 0.0f) && (az == 0.0f))) {


        recipNorm = inv_sqrt(ax * ax + ay * ay + az * az);
        ax *= recipNorm;
        ay *= recipNorm;
        az *= recipNorm;


        _2q0 = 2.0f * q->w;
        _2q1 = 2.0f * q->x;
        _2q2 = 2.0f * q->y;
        _2q3 = 2.0f * q->z;
        _4q0 = 4.0f * q->w;
        _4q1 = 4.0f * q->x;
        _4q2 = 4.0f * q->y;
        _8q1 = 8.0f * q->x;
        _8q2 = 8.0f * q->y;
        q0q0 = q->w * q->w;
        q1q1 = q->x * q->x;
        q2q2 = q->y * q->y;
        q3q3 = q->z * q->z;


        s0 = _4q0 * q2q2 + _2q2 * ax + _4q0 * q1q1 - _2q1 * ay;
        s1 = _4q1 * q3q3 - _2q3 * ax + 4.0f * q0q0 * q->x - _2q0 * ay - _4q1 + _8q1 * q1q1 + _8q1 * q2q2 + _4q1 * az;
        s2 = 4.0f * q0q0 * q->y + _2q0 * ax + _4q2 * q3q3 - _2q3 * ay - _4q2 + _8q2 * q1q1 + _8q2 * q2q2 + _4q2 * az;
        s3 = 4.0f * q1q1 * q->z - _2q1 * ax + 4.0f * q2q2 * q->z - _2q2 * ay;
        recipNorm = inv_sqrt(s0 * s0 + s1 * s1 + s2 * s2 + s3 * s3); 
        s0 *= recipNorm;
        s1 *= recipNorm;
        s2 *= recipNorm;
        s3 *= recipNorm;


        qDot1 -= beta_g * s0;
        qDot2 -= beta_g * s1;
        qDot3 -= beta_g * s2;
        qDot4 -= beta_g * s3;
    }

    q->w += qDot1 * dt;
    q->x += qDot2 * dt;
    q->y += qDot3 * dt;
    q->z += qDot4 * dt;


    recipNorm = inv_sqrt(q->w * q->w + q->x * q->x + q->y * q->y + q->z * q->z);
    q->w *= recipNorm;
    q->x *= recipNorm;
    q->y *= recipNorm;
    q->z *= recipNorm;
}
