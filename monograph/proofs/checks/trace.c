/* Replay one common control sequence. Reference states are row vectors.
 * The integer transition table has been checked exhaustively against the
 * packed C implementation and the independent arbitrary-precision update.
 * The real-valued reference error reported here is a binary64 diagnostic. */
#include <math.h>
#include <stddef.h>
#include <stdint.h>
#include <string.h>

double trace_maximum(size_t length, const uint8_t *tokens, size_t dimension,
                     const double *initial, const double *matrices,
                     const double *readout, const int32_t *transition,
                     const double *machine_output) {
    double state[dimension], next[dimension];
    memcpy(state, initial, dimension * sizeof(double));
    int32_t index = 0;
    double maximum = 0.0;
    for (size_t t = 0; t < length; ++t) {
        size_t a = tokens[t];
        const double *matrix = matrices + a * dimension * dimension;
        double answer = 0.0;
        for (size_t j = 0; j < dimension; ++j) {
            double sum = 0.0;
            for (size_t i = 0; i < dimension; ++i)
                sum += state[i] * matrix[i * dimension + j];
            next[j] = sum;
            answer += sum * readout[j];
        }
        memcpy(state, next, dimension * sizeof(double));
        index = transition[2 * index + a];
        double error = fabs(answer - machine_output[index]);
        if (!isfinite(error)) return NAN;
        if (error > maximum) maximum = error;
    }
    return maximum;
}
