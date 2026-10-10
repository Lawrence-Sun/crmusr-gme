// Compatibility entry point for the earlier direct-sampling experiment.
// The main generator now implements this method; do not duplicate the model.
#include "generate_michel.C"

void generate_michel_direct(const char *input = "pmu_output/polarizationData.root",
                            const char *outputDirectory = "michel_direct_output",
                            long long eventsPerComponent = 100000000LL,
                            unsigned int seed = 20241025, bool overwrite = false) {
    generate_michel(input, outputDirectory, eventsPerComponent, seed, overwrite);
}
