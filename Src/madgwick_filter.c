// madgwick_filter.c
#include "madgwick_filter.h"

static double beta_g = 0.1;
static double q0, q1, q2, q3;


void madgwick_init(QuaternionDouble* q) {
    q->w = 1.0f;
    q->x = 0.0f;
    q->y = 0.0f;
    q->z = 0.0f;
}

static double inv_sqrt(double x) {
	float halfx = 0.5f * x;
	float y = x;
	long i = *(long*)&y;
	i = 0x5f3759df - (i>>1);
	y = *(float*)&i;
	y = y * (1.5f - (halfx * y * y));
	return y;
}

void madgwick_update(QuaternionDouble* q, double ax, double ay, double az, 
                     double gx, double gy, double gz, double dt) {
    if (dt <= 0.0) return;

    q0 = q->w;
    q1 = q->x ;
    q2 = q->y ;
    q3 = q->z ;


    double recipNorm;
    double s0, s1, s2, s3;
    double qDot1, qDot2, qDot3, qDot4;
    double _2q0, _2q1, _2q2, _2q3, _4q0, _4q1, _4q2, _8q1, _8q2;
    double q0q0, q1q1, q2q2, q3q3;

    qDot1 = 0.5 * (-q1 * gx - q2 * gy - q3 * gz);
    qDot2 = 0.5 * (q0 * gx + q2 * gz - q3 * gy);
    qDot3 = 0.5 * (q0 * gy - q1 * gz + q3 * gx);
    qDot4 = 0.5 * (q0 * gz + q1 * gy - q2 * gx);


    if ((ax * ax + ay * ay + az * az) > 1e-12) {

        recipNorm = inv_sqrt(ax * ax + ay * ay + az * az);
        ax *= recipNorm;
        ay *= recipNorm;
        az *= recipNorm;


        _2q0 = 2.0 * q0;
        _2q1 = 2.0 * q1;
        _2q2 = 2.0 * q2;
        _2q3 = 2.0 * q3;
        _4q0 = 4.0 * q0;
        _4q1 = 4.0 * q1;
        _4q2 = 4.0 * q2;
        _8q1 = 8.0 * q1;
        _8q2 = 8.0 * q2;
        q0q0 = q0 * q0;
        q1q1 = q1 * q1;
        q2q2 = q2 * q2;
        q3q3 = q3 * q3;

        s0 = _4q0 * q2q2 + _2q2 * ax + _4q0 * q1q1 - _2q1 * ay;
        s1 = _4q1 * q3q3 - _2q3 * ax + 4.0 * q0q0 * q1 - _2q0 * ay - _4q1 + _8q1 * q1q1 + _8q1 * q2q2 + _4q1 * az;
        s2 = 4.0 * q0q0 * q2 + _2q0 * ax + _4q2 * q3q3 - _2q3 * ay - _4q2 + _8q2 * q1q1 + _8q2 * q2q2 + _4q2 * az;
        s3 = 4.0 * q1q1 * q3 - _2q1 * ax + 4.0 * q2q2 * q3 - _2q2 * ay;


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


    q0 += qDot1 * dt;
    q1 += qDot2 * dt;
    q2 += qDot3 * dt;
    q3 += qDot4 * dt;

    recipNorm = inv_sqrt(q0 * q0 + q1 * q1 + q2 * q2 + q3 * q3);
    q0 *= recipNorm;
    q1 *= recipNorm;
    q2 *= recipNorm;
    q3 *= recipNorm;


    q->w = q0;
    q->x = q1;
    q->y = q2;
    q->z = q3;
}
