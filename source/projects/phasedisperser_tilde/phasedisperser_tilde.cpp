#include "c74_min.h"

#include <array>
#include <cmath>

using namespace c74::min;

struct AllPassCoefficients {
    double c0 = 0.0f;
    double c1 = 0.0f;
    double c2 = 0.0f;
    double c3 = 0.0f;
    double c4 = 0.0f;
};

class AllPassFilter {
public:
    AllPassFilter() = default;

    AllPassFilter(double frequency, double sr, double q) {
        double w0 = 2 * 3.14159265358979323846 * frequency / sr;
        double cosw0 = cos(w0);
        double alpha = sin(w0) * (2 * q);

        double a1 = 0;
        double a2 = 0;
        double b0 = 1 - alpha;
        double b1 = -2 * cosw0;
        double b2 = 1 + alpha;
        double a0 = 1 + alpha;
        a1 = -2 * cosw0;
        a2 = 1 - alpha;
        
        co.c0 = b0 / a0;
        co.c1 = b1 / a0;
        co.c2 = b2 / a0;
        co.c3 = a1 / a0;
        co.c4 = a2 / a0;
    }

    void copy_coefficients_from(const AllPassFilter& filter) {
        co.c0 = filter.co.c0;
        co.c1 = filter.co.c1;
        co.c2 = filter.co.c2;
        co.c3 = filter.co.c3;
        co.c4 = filter.co.c4;
    }
    
    void reset() {
        xm1 = 0;
        xm2 = 0;
        ym1 = 0;
        ym2 = 0;
    }

    void processBlock(const double* in, double* out, long numSamples) {
        double x = 0;
        double y = 0;

        for (auto i = 0; i < numSamples; i++) {
            x = in[i];
            y = co.c0 * x + 
                co.c1 * xm1 +
                co.c2 * xm2 -
                co.c3 * ym1 -
                co.c4 * ym2;
            
            ym2 = ym1;
            ym1 = y;
            xm2 = xm1;
            xm1 = x;

            out[i] = y;
        }
    }

private:
    AllPassCoefficients co{};

    double xm1 = 0.0f;
    double xm2 = 0.0f;
    double ym1 = 0.0f;
    double ym2 = 0.0f;
};

class phasedisperser : public object<phasedisperser>, public vector_operator<> {
public:
    MIN_DESCRIPTION {"A massive allpass filter designed for phase dispersion"};
    MIN_TAGS {"filter"};
    MIN_AUTHOR {"Isabel Kaspriskie"};

    inlet<> in1 {this, "(signal) Input1"};
    inlet<> in2 {this, "(signal) Input2"};
    outlet<> out1 {this, "(signal) Output1", "signal"};
    outlet<> out2 {this, "(signal) Output2", "signal"};

    attribute<number, threadsafe::no, limit::clamp> a_frequency {
        this, "frequency", 700, range {20, 20000}
    };

    // "number of iterations"
    attribute<number, threadsafe::no, limit::clamp> a_intensity {
        this, "intensity", 25, range {1, 50}
    };

    attribute<number, threadsafe::no, limit::clamp> a_q {
        this, "q", 0.7, range {0.0, 1.41}
    };

    attribute<number, threadsafe::no, limit::clamp> a_mix {
        this, "mix", 0.8, range {0.0, 1.0}
    };

    message<> dspsetup {this, "dspsetup",
        MIN_FUNCTION {
            return {};
        }
    };

    void operator()(audio_bundle input, audio_bundle output) override {
        double* in1 = input.samples(0);
        double* in2 = input.samples(1);
        double* out1 = output.samples(0);
        double* out2 = output.samples(1);

        long sampleFrames = input.frame_count();

        auto iterations = static_cast<long>(a_intensity);
        if (iterations > cur_iterations) {
            if (cur_iterations > 0) {
                for (int i = cur_iterations - 1; i < iterations; i++) {
                    filterL[i].copy_coefficients_from(filterL[cur_iterations - 1]);
                    filterR[i].copy_coefficients_from(filterR[cur_iterations - 1]);
                }
            } else {
                setupFilters();
            }
            cur_iterations = iterations;
        } else if (iterations < cur_iterations) {
            cur_iterations = iterations;
        }

        // copy input to temp array
        for (int i = 0; i < sampleFrames; i++) {
            temp1[i] = in1[i];
            temp2[i] = in2[i];
            if (abs(temp1[i]) >= noiseFloor || abs(temp2[i]) >= noiseFloor) {
                samplesSinceSilence = 0;
            }
        }

        // filter the audio
        if (samplesSinceSilence < deactivateAfterSamples && cur_iterations != 0 && a_mix > 0) {
            for (int i = 0; i < cur_iterations; i++) {
                filterL[i].processBlock(temp1.data(), left.data(), sampleFrames);
                filterR[i].processBlock(temp2.data(), right.data(), sampleFrames);

                for (int j = 0; j < sampleFrames; j++) {
                    temp1[j] = left[j];
                    temp2[j] = right[j];
                }
            }
        }

        for (int i = 0; i < sampleFrames; i++) {
            if (abs(temp1[i]) >= noiseFloor || abs(temp2[i]) >= noiseFloor) {
                samplesSinceSilence = 0;
            } else if (samplesSinceSilence < 32768) {
                // int overflow protection
                samplesSinceSilence++;
            }

            out1[i] = temp1[i] * a_mix + in1[i] * (1 - a_mix);
            out2[i] = temp2[i] * a_mix + in2[i] * (1 - a_mix);
        }
    }

private:
    static constexpr long maxFilters = 50;
    static constexpr long maxBufferSize = 4096;
    static constexpr long deactivateAfterSamples = 16384;
    static constexpr double noiseFloor = 0.000007;

    std::array<AllPassFilter, maxFilters> filterL;
    std::array<AllPassFilter, maxFilters> filterR;
    std::array<double, maxBufferSize> temp1;
    std::array<double, maxBufferSize> temp2;
    std::array<double, maxBufferSize> left;
    std::array<double, maxBufferSize> right;

    int freq = 0;
    long cur_iterations = 0;
    double q = 0;
    long samplesSinceSilence = 1;
    double last_freq = 0;

    void setupFilters() {
        // https://www.musicdsp.org/en/latest/Other/260-exponential-curve-for.html
        freq = floor(exp((16 + a_frequency * 100 * 1.20103)*log(1.059)) * 8.17742);
        q = a_q * sqrt(2);

        if (q <= 0.005) {
            q = 0.005;
        }

        filterL[0] = AllPassFilter(freq, 44100.0f, q);
        filterR[0] = AllPassFilter(freq, 44100.0f, q);

        for (int i = 1; i < a_intensity; i++) {
            filterL[i].copy_coefficients_from(filterL[0]);
            filterR[i].copy_coefficients_from(filterR[0]);
            // attempt to prevent the filters from generating noise that could damage audio equipment
            if (abs(a_frequency - last_freq) > a_frequency / 10 && freq < 500) {
                filterL[i].reset();
                filterR[i].reset();
            }
        }
    }
};


MIN_EXTERNAL(phasedisperser);
