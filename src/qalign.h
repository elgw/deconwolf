#pragma once
#include <stdint.h>

typedef struct
{
    union{
        float delta[3];
        struct{
            float dx;
            float dy;
            float dz;
        };
    };
    int32_t idx1;
    int32_t idx2;
} qalign_pair;


typedef struct {
    /* Input data */
    const float * A; /* 3*nA elements, x, y, z for each point */
    const float * B;
    int64_t nA;
    int64_t nB;

    /* Settings */
    int verbose;
    /* Largest expected delative difference in magnification. For example 1/1000  */
    float rel_error;
    /* Standard deviation of the localization accuracy
     * given in pixels. */
    float localication_sigma;
    /* Number of points to use at most */
    int64_t npoint;

    /* Results */
    qalign_pair * H; /* Hypotheses */
    float * S; /* Scaling */
    float * SZ;
    int64_t nH; /* Number of hypotheses */
} qalign_config;

qalign_config * qalign_config_new();

int qalign(qalign_config *);

void qalign_config_free(qalign_config *);
