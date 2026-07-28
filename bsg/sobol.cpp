// Sobol low-discrepancy sequence generator (see sobol.hpp for the sizing rationale).
#include "sobol.hpp"
#include <stdexcept>
#include <string>

namespace StochasticSimulator {
namespace {

/**
 * Direction-number table.
 *
 * Each row is a primitive polynomial x^s + a_1 x^(s-1) + ... + a_(s-1) x + 1, with `a` packing
 * the inner coefficients a_1..a_(s-1) (a_1 = MSB), plus s initial direction numbers m_1..m_s.
 * Rows are ordered by polynomial degree then by increasing `a`, the conventional Sobol ordering.
 *
 * Dimensions 1-8 carry the classic tuned initial numbers. Dimensions 9-19 use m_i = 1, which is
 * a *valid* Sobol construction (every m_i is odd and m_i < 2^i) but is not uniformity-tuned, so
 * its low-discrepancy constant is worse than a tabulated row's. If you need many decorrelated
 * operands and care about that, pass the matching line of the official Joe-Kuo table to
 * SobolRNG::set_direction_numbers(). Dimension 1 is the van der Corput / bit-reversal sequence
 * and needs no polynomial at all.
 */
struct DirectionRow {
    unsigned s;
    uint32_t a;
    uint32_t m[6];
};

// Indexed by (dimension - 2); dimension 1 is special-cased.
constexpr DirectionRow DIRECTION_TABLE[] = {
    {1,  0,  {1}},                  // dim 2   x + 1
    {2,  1,  {1, 3}},               // dim 3   x^2 + x + 1
    {3,  1,  {1, 3, 1}},            // dim 4   x^3 + x + 1
    {3,  2,  {1, 1, 1}},            // dim 5   x^3 + x^2 + 1
    {4,  1,  {1, 1, 3, 3}},         // dim 6   x^4 + x + 1
    {4,  4,  {1, 3, 5, 13}},        // dim 7   x^4 + x^3 + 1
    {5,  2,  {1, 1, 5, 5, 17}},     // dim 8   x^5 + x^2 + 1
    {5,  4,  {1, 1, 1, 1, 1}},      // dim 9   x^5 + x^3 + 1            (untuned m)
    {5,  7,  {1, 1, 1, 1, 1}},      // dim 10  x^5 + x^3 + x^2 + x + 1  (untuned m)
    {5, 11,  {1, 1, 1, 1, 1}},      // dim 11  x^5 + x^4 + x^2 + x + 1  (untuned m)
    {5, 13,  {1, 1, 1, 1, 1}},      // dim 12  x^5 + x^4 + x^3 + x + 1  (untuned m)
    {5, 14,  {1, 1, 1, 1, 1}},      // dim 13  x^5 + x^4 + x^3 + x^2+ 1 (untuned m)
    {6,  1,  {1, 1, 1, 1, 1, 1}},   // dim 14  x^6 + x + 1              (untuned m)
    {6, 13,  {1, 1, 1, 1, 1, 1}},   // dim 15  x^6 + x^4 + x^3 + x + 1  (untuned m)
    {6, 16,  {1, 1, 1, 1, 1, 1}},   // dim 16  x^6 + x^5 + 1            (untuned m)
    {6, 19,  {1, 1, 1, 1, 1, 1}},   // dim 17  x^6 + x^5 + x^2 + x + 1  (untuned m)
    {6, 22,  {1, 1, 1, 1, 1, 1}},   // dim 18  x^6 + x^5 + x^3 + x^2+ 1 (untuned m)
    {6, 25,  {1, 1, 1, 1, 1, 1}},   // dim 19  x^6 + x^5 + x^4 + x + 1  (untuned m)
};

// Index of the least significant ZERO bit -- the bit the Gray-code recurrence flips next.
inline unsigned least_significant_zero(uint64_t value) {
    unsigned index = 0;
    while (value & 1ull) {
        value >>= 1;
        ++index;
    }
    return index;
}

unsigned width_of(StreamLength lengthMode) {
    switch (lengthMode) {
        case StreamLength::Length_128:   return 7;
        case StreamLength::Length_256:   return 8;
        case StreamLength::Length_512:   return 9;
        case StreamLength::Length_1024:  return 10;
        case StreamLength::Length_4096:  return 12;
        case StreamLength::Length_16384: return 14;
        case StreamLength::Length_65536: return 16;
    }
    throw std::invalid_argument("Unknown StreamLength mode.");
}

// Builds v_1..v_width from one primitive-polynomial row. Any odd m_i < 2^i is a legal seed.
std::vector<uint32_t> build_vectors(unsigned width, unsigned s, uint32_t a,
                                    const uint32_t* m_init, std::size_t m_len) {
    if (s == 0 || m_len != s) {
        throw std::invalid_argument("Expected exactly s initial direction numbers.");
    }
    for (unsigned i = 1; i <= s; ++i) {
        uint32_t seed = m_init[i - 1];
        if ((seed & 1u) == 0) {
            throw std::invalid_argument("Initial direction numbers must be odd.");
        }
        if (i < 32 && seed >= (1u << i)) {
            throw std::invalid_argument("Initial direction number m_i must be less than 2^i.");
        }
    }

    // m_i, 1-based; extend past s with the primitive-polynomial recurrence:
    //   m_i = 2 a_1 m_(i-1) ^ 4 a_2 m_(i-2) ^ ... ^ 2^(s-1) a_(s-1) m_(i-s+1) ^ 2^s m_(i-s) ^ m_(i-s)
    std::vector<uint32_t> m(width + 1, 0);
    for (unsigned i = 1; i <= s && i <= width; ++i) {
        m[i] = m_init[i - 1];
    }
    for (unsigned i = s + 1; i <= width; ++i) {
        m[i] = m[i - s] ^ (m[i - s] << s);
        for (unsigned k = 1; k <= s - 1; ++k) {
            m[i] ^= ((a >> (s - 1 - k)) & 1u) * (m[i - k] << k);
        }
    }

    // v_i = m_i scaled to the top of the word: m_i * 2^(width - i).
    std::vector<uint32_t> v(width, 0);
    for (unsigned i = 1; i <= width; ++i) {
        v[i - 1] = m[i] << (width - i);
    }
    return v;
}

// Direction vectors for a tabulated dimension at a given width.
std::vector<uint32_t> vectors_for_dimension(unsigned width, unsigned dim) {
    if (dim == 1) {
        // van der Corput: v_i = 2^(width - i), i.e. plain bit reversal of the counter.
        std::vector<uint32_t> v(width, 0);
        for (unsigned i = 1; i <= width; ++i) {
            v[i - 1] = 1u << (width - i);
        }
        return v;
    }
    const DirectionRow& row = DIRECTION_TABLE[dim - 2];
    return build_vectors(width, row.s, row.a, row.m, row.s);
}

/**
 * Below a certain width the polynomial recurrence never runs -- only the s seed values reach the
 * direction vectors -- so several tabulated rows collapse onto the same sequence. Two "different"
 * dimensions that are secretly identical are perfectly correlated and would silently destroy any
 * multiply, so a dimension is only usable once it differs from every LOWER dimension. Lower
 * dimensions are canonical: dimension 1 is always valid.
 */
bool collides_with_lower_dimension(unsigned width, unsigned dim, unsigned& collides_with) {
    std::vector<uint32_t> mine = vectors_for_dimension(width, dim);
    for (unsigned other = 1; other < dim; ++other) {
        if (vectors_for_dimension(width, other) == mine) {
            collides_with = other;
            return true;
        }
    }
    return false;
}

}  // namespace

SobolRNG::SobolRNG(unsigned width_bits, unsigned dim) {
    if (width_bits == 0 || width_bits > MAX_WIDTH) {
        throw std::invalid_argument("Sobol width must be between 1 and 32 bits.");
    }
    width = width_bits;
    period = 1ull << width;

    load_tabulated_dimension(dim);
    reset();
}

SobolRNG::SobolRNG(StreamLength lengthMode, unsigned dim)
    : SobolRNG(width_of(lengthMode), dim) {}

unsigned SobolRNG::width_for_length(uint64_t stream_length) {
    if (stream_length == 0) {
        throw std::invalid_argument("Stream length must be at least 1.");
    }
    if (stream_length > (1ull << MAX_WIDTH)) {
        throw std::invalid_argument("Stream length exceeds the 2^32 maximum Sobol period.");
    }

    // Smallest width whose period holds the whole stream: ceil(log2(stream_length)).
    unsigned bits = 1;
    while ((1ull << bits) < stream_length) {
        ++bits;
    }
    return bits;
}

SobolRNG SobolRNG::for_length(uint64_t stream_length, unsigned dim) {
    return SobolRNG(width_for_length(stream_length), dim);
}

void SobolRNG::load_tabulated_dimension(unsigned dim) {
    if (dim == 0 || dim > MAX_DIMENSION) {
        throw std::invalid_argument(
            "Sobol dimension must be between 1 and " + std::to_string(MAX_DIMENSION) +
            "; use set_direction_numbers() to supply your own row beyond that.");
    }
    dimension = dim;
    direction = vectors_for_dimension(width, dim);

    // Refuse a dimension that is a secret duplicate of a lower one at this width -- see
    // collides_with_lower_dimension(). Silently returning a perfectly correlated stream would be
    // far worse than failing here.
    unsigned duplicate_of = 0;
    if (collides_with_lower_dimension(width, dim, duplicate_of)) {
        throw std::invalid_argument(
            "Sobol dimension " + std::to_string(dim) + " is identical to dimension " +
            std::to_string(duplicate_of) + " at width " + std::to_string(width) +
            " (streams would be perfectly correlated). It needs width >= " +
            std::to_string(minimum_width_for_dimension(dim)) +
            "; widen the stream or pick a lower dimension.");
    }
}

unsigned SobolRNG::minimum_width_for_dimension(unsigned dim) {
    if (dim == 0 || dim > MAX_DIMENSION) {
        throw std::invalid_argument("Sobol dimension must be between 1 and " +
                                    std::to_string(MAX_DIMENSION) + ".");
    }
    for (unsigned w = 1; w <= MAX_WIDTH; ++w) {
        unsigned duplicate_of = 0;
        if (!collides_with_lower_dimension(w, dim, duplicate_of)) {
            return w;
        }
    }
    throw std::logic_error("Sobol dimension never becomes distinct.");  // unreachable
}

void SobolRNG::set_direction_numbers(unsigned s, uint32_t a, const std::vector<uint32_t>& m_init) {
    build_direction_vectors(s, a, m_init);
    reset();
}

void SobolRNG::build_direction_vectors(unsigned s, uint32_t a, const std::vector<uint32_t>& m_init) {
    direction = build_vectors(width, s, a, m_init.data(), m_init.size());
}

uint32_t SobolRNG::next() {
    uint32_t point = state;

    ++counter;
    if (counter >= period) {
        // One full period consumed: the sequence restarts cleanly at 0.
        counter = 0;
        state = 0;
    } else {
        // Gray-code recurrence: flip the direction vector for the bit that just changed.
        state ^= direction[least_significant_zero(counter - 1)];
    }

    return point;
}

double SobolRNG::next_uniform() {
    return static_cast<double>(next()) / static_cast<double>(period);
}

bool SobolRNG::next_bit(double probability) {
    if (probability < 0.0 || probability > 1.0) {
        throw std::invalid_argument("Probability must be between 0.0 and 1.0");
    }

    // Same comparator semantics as BitstreamGenerator: truncating threshold, strict ">".
    // p = 1.0 gives threshold == period, which beats every point, so the stream is all ones.
    uint64_t threshold = static_cast<uint64_t>(probability * static_cast<double>(period));
    return threshold > static_cast<uint64_t>(next());
}

uint32_t SobolRNG::value_at(uint64_t index) const {
    uint64_t gray = index ^ (index >> 1);
    uint32_t point = 0;
    for (unsigned i = 0; i < width && gray != 0; ++i, gray >>= 1) {
        if (gray & 1ull) {
            point ^= direction[i];
        }
    }
    return point;
}

void SobolRNG::jump_to(uint64_t index) {
    counter = index % period;
    state = value_at(counter);
}

uint64_t SobolRNG::get_max_cycles() const { return period; }
unsigned SobolRNG::get_width() const { return width; }
unsigned SobolRNG::get_dimension() const { return dimension; }
uint64_t SobolRNG::get_position() const { return counter; }

void SobolRNG::reset() {
    counter = 0;
    state = 0;  // the Sobol sequence starts at the origin
}

std::vector<bool> generate_sobol_stream(double probability, std::size_t length, unsigned dimension) {
    SobolRNG rng = SobolRNG::for_length(static_cast<uint64_t>(length), dimension);

    std::vector<bool> stream;
    stream.reserve(length);
    for (std::size_t i = 0; i < length; ++i) {
        stream.push_back(rng.next_bit(probability));
    }
    return stream;
}

}  // namespace StochasticSimulator
