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
    const float * A;
    const float * B;
    int64_t nA;
    int64_t nB;
    float similarity_threshold;
    float deviation_x;
    int verbose;
    int64_t npoint;

    /* Results */
    qalign_pair * H; /* Hypotheses */
    int64_t nH; /* Number of hypotheses */
} qalign_config;

qalign_config * qalign_config_new();

int qalign(qalign_config *);

void qalign_config_free(qalign_config *);
