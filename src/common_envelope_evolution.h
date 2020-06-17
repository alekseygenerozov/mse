#include "types.h"
extern "C"
{

int common_envelope_evolution(ParticlesMap *particlesMap, int binary_index, int index1, int index2, double t, int *integration_flag);
int binary_common_envelope_evolution(ParticlesMap *particlesMap, int binary_index, int index1, int index2, double t, int *integration_flag);
int triple_common_envelope_evolution(ParticlesMap *particlesMap, int binary_index, int index1, int index2, double t, int *integration_flag);

}
